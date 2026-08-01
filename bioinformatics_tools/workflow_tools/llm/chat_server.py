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
# Last time anyone touched this server. The page sends a stop on exit, but that
# never arrives if the browser crashes, the laptop sleeps, or the VPN drops --
# exactly the cases where a GPU would otherwise sit idle until walltime. This
# is the backstop that makes those cases self-correcting.
LAST_SEEN = time.time()
# Generation is serialised: one set of weights, and concurrent .generate() calls
# on the same model would contend for the GPU and can interleave KV cache state.
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


def generate(system: str, prompt: str, max_tokens: int) -> str:
    import torch

    # Llama 3.1 instruct chat template. Using the tokenizer's own template
    # rather than hand-rolling the special tokens keeps this correct if the
    # fused model ships a different one.
    messages = [{"role": "system", "content": system},
                {"role": "user", "content": prompt}]
    try:
        text = TOKENIZER.apply_chat_template(
            messages, tokenize=False, add_generation_prompt=True)
    except Exception:
        text = f"{system}\n\n{prompt}\n\n"

    inputs = TOKENIZER(text, return_tensors="pt").to(MODEL.device)
    with _GEN_LOCK, torch.no_grad():
        out = MODEL.generate(
            **inputs,
            max_new_tokens=max_tokens,
            do_sample=False,              # deterministic: same evidence -> same answer
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
            # A health poll counts as activity: while the page is open it polls,
            # so the idle clock only advances once nobody is watching.
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
            # Fast path: the page was closed. Ack first, then stop — otherwise
            # the caller (often a sendBeacon during unload) sees a dropped
            # connection and cannot tell success from failure.
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
        try:
            text = generate(
                body.get("system") or "",
                prompt,
                int(body.get("max_tokens") or 600),
            )
            self._send(200, {"text": text})
        except Exception as exc:
            log(f"generation failed: {type(exc).__name__}: {exc}")
            self._send(500, {"error": f"{type(exc).__name__}: {exc}"})

    def log_message(self, fmt, *args):
        log(fmt % args)


def advertise(path: Path, host: str, port: int, model: str) -> None:
    """Publish the endpoint atomically so a reader never sees a half-written file."""
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
    # 5 minutes, not 30. The panel polls /health every 60s while it is open, so
    # silence for five minutes means nobody is watching -- and that signal is
    # reliable in a way page-exit events are not. pagehide does NOT fire on
    # SvelteKit client-side navigation, so "user clicked back" never reached the
    # stop endpoint and a GPU sat held until the old 30-minute timeout. This is
    # the mechanism that actually frees the GPU; the stop endpoint is only a
    # fast path for when it happens to work.
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
        # Only remove the advert if it is still ours; a newer server may have
        # replaced it, and deleting that would take chat offline for no reason.
        try:
            if adv.is_file() and json.loads(adv.read_text()).get("pid") == os.getpid():
                adv.unlink()
                log("advert removed")
        except Exception:
            pass


if __name__ == "__main__":
    main()
