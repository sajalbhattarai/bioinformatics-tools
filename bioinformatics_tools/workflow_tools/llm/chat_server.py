#!/usr/bin/env python3
"""chat_server.py -- a minimal, long-lived inference endpoint for the genome chat.

Runs INSIDE llm.sif on a GPU node and holds the model in memory so the web app
can ask questions interactively. This is deliberately not vLLM or TGI: llm.sif
has neither, the GPUs here are AMD (the phase15 rule exports
PYTORCH_HIP_ALLOC_CONF, i.e. ROCm), and ROCm support in those servers is patchy.
The container does have a working torch + transformers, so this uses those plus
the standard library's own HTTP server -- no fastapi, no uvicorn, nothing to
install.

Distinct from score-genes-llm.py, which runs the SAME model as a batch job to
produce scores. This one only answers questions about results that already
exist; it never writes to a run.

Service discovery: a SLURM allocation is not a stable address -- the node
changes every time and the job dies at walltime. On startup this writes

    {"host": ..., "port": ..., "pid": ..., "model": ..., "started": ...}

to --advertise (an atomic replace), and removes it on clean shutdown. The API
reads that file to find the current endpoint, and treats "file missing" as
"chat is offline" rather than an error.

Endpoints:
    GET  /health  -> {"ok": true, "model": ...}
    POST /chat    -> {"system": str, "prompt": str, "max_tokens": int}
                     => {"text": str}

There is no auth here: it binds inside the cluster and the public-facing
authenticated surface is the API's /v1/llm/chat, which proxies to it.
"""
from __future__ import annotations

import argparse
import json
import os
import socket
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

MODEL = None
TOKENIZER = None
MODEL_NAME = ""
# Time of the last request; the idle watchdog frees the GPU when the page is gone.
LAST_SEEN = time.time()
# Serialises generation: concurrent .generate() calls would interleave KV cache state.
_GEN_LOCK = threading.Lock()


def log(msg: str) -> None:
    print(f"[chat_server] {msg}", flush=True)


def load_model(model_path: str, dtype: str) -> None:
    global MODEL, TOKENIZER, MODEL_NAME
    import torch
    from transformers import AutoModelForCausalLM, AutoTokenizer

    MODEL_NAME = model_path
    log(f"loading {model_path} (dtype={dtype}) …")
    t0 = time.time()
    TOKENIZER = AutoTokenizer.from_pretrained(model_path)
    if TOKENIZER.pad_token_id is None:
        TOKENIZER.pad_token = TOKENIZER.eos_token
    torch_dtype = {"bfloat16": torch.bfloat16, "float16": torch.float16,
                   "float32": torch.float32}[dtype]
    MODEL = AutoModelForCausalLM.from_pretrained(
        model_path, torch_dtype=torch_dtype, device_map="auto",
    )
    MODEL.eval()
    dev = next(MODEL.parameters()).device
    log(f"loaded on {dev} in {time.time() - t0:.0f}s")


def _build_inputs(system: str, prompt: str):
    """Returns tokenised chat-template inputs for the batch and streaming paths."""
    messages = [{"role": "system", "content": system},
                {"role": "user", "content": prompt}]
    try:
        text = TOKENIZER.apply_chat_template(
            messages, tokenize=False, add_generation_prompt=True)
    except Exception:
        text = f"{system}\n\n{prompt}\n\n"
    return TOKENIZER(text, return_tensors="pt").to(MODEL.device)


def generate_stream(system: str, prompt: str, max_tokens: int):
    """Yields text as it is generated, via a worker thread and TextIteratorStreamer."""
    import torch
    from transformers import TextIteratorStreamer

    inputs = _build_inputs(system, prompt)
    streamer = TextIteratorStreamer(TOKENIZER, skip_prompt=True,
                                    skip_special_tokens=True)
    kwargs = dict(**inputs, max_new_tokens=max_tokens, do_sample=False,
                  pad_token_id=TOKENIZER.pad_token_id, streamer=streamer)

    def run():
        with _GEN_LOCK, torch.no_grad():
            MODEL.generate(**kwargs)

    t = threading.Thread(target=run, daemon=True)
    t.start()
    for piece in streamer:
        if piece:
            yield piece
    t.join()


def generate(system: str, prompt: str, max_tokens: int) -> str:
    import torch

    inputs = _build_inputs(system, prompt)
    with _GEN_LOCK, torch.no_grad():
        out = MODEL.generate(
            **inputs,
            max_new_tokens=max_tokens,
            do_sample=False,              # deterministic output
            pad_token_id=TOKENIZER.pad_token_id,
        )
    # Slice off the prompt so only the completion is returned.
    return TOKENIZER.decode(out[0][inputs["input_ids"].shape[1]:],
                            skip_special_tokens=True).strip()


class Handler(BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"

    def _send(self, code: int, payload: dict) -> None:
        body = json.dumps(payload).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):  # noqa: N802
        global LAST_SEEN
        p = self.path.rstrip("/")
        if p == "/health":
            # Health polls from an open page count as activity.
            LAST_SEEN = time.time()
            self._send(200, {"ok": MODEL is not None, "model": MODEL_NAME,
                             "idle_s": 0})
        elif p == "/idle":
            self._send(200, {"idle_s": round(time.time() - LAST_SEEN)})
        else:
            self._send(404, {"error": "not found"})

    def do_POST(self):  # noqa: N802
        global LAST_SEEN
        p = self.path.rstrip("/")
        if p == "/shutdown":
            # Replies before stopping so the caller (often a sendBeacon) sees success.
            self._send(200, {"stopping": True})
            log("shutdown requested — exiting")
            threading.Thread(target=self.server.shutdown, daemon=True).start()
            return
        if p != "/chat":
            self._send(404, {"error": "not found"})
            return
        LAST_SEEN = time.time()
        try:
            n = int(self.headers.get("Content-Length") or 0)
            body = json.loads(self.rfile.read(n) or b"{}")
        except Exception as exc:
            self._send(400, {"error": f"bad request: {exc}"})
            return
        prompt = (body.get("prompt") or "").strip()
        if not prompt:
            self._send(400, {"error": "prompt is required"})
            return
        system = body.get("system") or ""
        max_tokens = int(body.get("max_tokens") or 600)

        if not body.get("stream"):
            try:
                self._send(200, {"text": generate(system, prompt, max_tokens)})
            except Exception as exc:
                log(f"generation failed: {type(exc).__name__}: {exc}")
                self._send(500, {"error": f"{type(exc).__name__}: {exc}"})
            return

        # Streams chunked UTF-8 text; headers go out before the first token.
        try:
            self.send_response(200)
            self.send_header("Content-Type", "text/plain; charset=utf-8")
            self.send_header("Transfer-Encoding", "chunked")
            self.send_header("Cache-Control", "no-cache")
            self.end_headers()
            for piece in generate_stream(system, prompt, max_tokens):
                data = piece.encode("utf-8")
                self.wfile.write(b"%X\r\n" % len(data) + data + b"\r\n")
                self.wfile.flush()
            self.wfile.write(b"0\r\n\r\n")
            self.wfile.flush()
        except (BrokenPipeError, ConnectionResetError):
            log("client disconnected during stream")
        except Exception as exc:
            log(f"stream failed: {type(exc).__name__}: {exc}")
            # Headers are already sent, so the only signal left is the body.
            try:
                msg = f"\n[stream error: {type(exc).__name__}: {exc}]".encode()
                self.wfile.write(b"%X\r\n" % len(msg) + msg + b"\r\n0\r\n\r\n")
            except Exception:
                pass

    def log_message(self, fmt, *args):
        log(fmt % args)


def advertise(path: Path, host: str, port: int, model: str) -> None:
    """Writes the endpoint file atomically (temp file + replace)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps({
        "host": host, "port": port, "pid": os.getpid(),
        "model": model, "started": time.time(),
    }) + "\n")
    tmp.replace(path)
    log(f"advertised {host}:{port} -> {path}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--model", required=True, help="path to the model directory")
    ap.add_argument("--advertise", required=True,
                    help="where to publish host/port for the API to discover")
    ap.add_argument("--port", type=int, default=0, help="0 = pick a free port")
    ap.add_argument("--dtype", default="bfloat16",
                    choices=["bfloat16", "float16", "float32"])
    # The page polls /health every 60 s, so 5 min of silence means it is closed.
    ap.add_argument("--idle-timeout", type=int, default=300,
                    help="exit after this many seconds with no /chat or /health "
                         "(0 disables). Releases the GPU when the page is gone.")
    args = ap.parse_args()

    if not Path(args.model).is_dir():
        sys.exit(f"model directory not found: {args.model}")

    load_model(args.model, args.dtype)

    srv = ThreadingHTTPServer(("0.0.0.0", args.port), Handler)
    host = socket.getfqdn()
    port = srv.server_address[1]
    adv = Path(args.advertise)
    advertise(adv, host, port, args.model)

    if args.idle_timeout > 0:
        def watchdog():
            while True:
                time.sleep(15)
                idle = time.time() - LAST_SEEN
                if idle >= args.idle_timeout:
                    log(f"idle {idle:.0f}s >= {args.idle_timeout}s — releasing the GPU")
                    srv.shutdown()
                    return
        threading.Thread(target=watchdog, daemon=True).start()
        log(f"idle timeout: {args.idle_timeout}s")

    log("ready — POST /chat, POST /shutdown, GET /health")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        # Removes the advert only if a newer server has not replaced it.
        try:
            if adv.is_file() and json.loads(adv.read_text()).get("pid") == os.getpid():
                adv.unlink()
                log("advert removed")
        except Exception:
            pass


if __name__ == "__main__":
    main()
