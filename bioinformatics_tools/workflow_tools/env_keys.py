"""Creates the backend secret keys before the licensing gate imports the API layer.

Missing BSP_SECRET_KEY / BSP_ENCRYPTION_KEY are generated and appended to the backend
.env that dane-api reads; existing keys are never overwritten. Uses only the stdlib
(cryptography lazily) so importing it does not load the API package.
"""
from __future__ import annotations

import os
import secrets
from pathlib import Path

_KEYS = ("BSP_SECRET_KEY", "BSP_ENCRYPTION_KEY")


def _repo_root() -> Path:
    # Same root as api/main.py's _PROJECT_ROOT, so the .env written is the one it reads.
    return Path(__file__).resolve().parents[2]


def _read_env(path: Path) -> dict:
    values: dict[str, str] = {}
    try:
        for line in path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                key, _, val = line.partition("=")
                values[key.strip()] = val.strip()
    except OSError:
        pass
    return values


def _generate(key: str) -> str:
    if key == "BSP_ENCRYPTION_KEY":
        from cryptography.fernet import Fernet
        return Fernet.generate_key().decode()
    return secrets.token_urlsafe(32)


def ensure_api_keys() -> None:
    """Sets both keys in os.environ, generating and appending any missing ones to .env (best effort)."""
    env_path = _repo_root() / ".env"
    existing = _read_env(env_path)
    added: dict[str, str] = {}
    for key in _KEYS:
        val = os.environ.get(key) or existing.get(key)
        if not val:
            val = _generate(key)
            added[key] = val
        os.environ.setdefault(key, val)

    if not added:
        return
    try:
        prefix = ""
        if env_path.exists() and env_path.stat().st_size > 0:
            if not env_path.read_text(encoding="utf-8").endswith("\n"):
                prefix = "\n"
        with env_path.open("a", encoding="utf-8") as f:
            f.write(prefix + "".join(f"{k}={v}\n" for k, v in added.items()))
        env_path.chmod(0o600)
    except OSError:
        pass  # the keys are still set in os.environ for this run
