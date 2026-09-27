"""
Centralized SSH connection configuration.

Provides a single place to manage host, username, and connection setup
instead of hardcoding paramiko boilerplate in every function.

All SSH/SFTP operations are API-layer only. The CLI runs directly on the
cluster and has no need to establish outbound SSH connections.

API usage:
    Call make_user_connection(host, username, private_key_str), which reads
    the user's cluster credentials from the database record, loads the
    decrypted private key into memory (never written to disk), and returns
    a ready SSHConnection.
"""
import io
import logging
import os
import shlex
import threading
import time

import paramiko

LOGGER = logging.getLogger(__name__)

# The repository margie-frontend's hpc-connect.sh clones for the API
# (BACKEND_REPO_URL), so the workflow checkout and the API share one codebase.
_DANE_WF_REPO_URL = 'https://github.com/sajalbhattarai/margie-backend.git'
# Branch that margie_sb remote checkouts track (override with BSP_MARGIE_SB_REF).
# The margie workflow is never auto-synced, so this applies only to margie_sb.
_MARGIE_SB_REF = os.getenv('BSP_MARGIE_SB_REF', 'for-website-deployment')

# ---- Connection pool ----
# Clients are reused per (host, username, key) to avoid a full SSH handshake per
# request; a dead client is discarded and replaced on the next request.
_POOL: dict = {}
_POOL_LOCK = threading.Lock()
# A pooled client is kept for as long as it works (an unclosed paramiko
# transport keeps its thread alive); keepalives let a dropped one be replaced.
_KEEPALIVE_SECONDS = 30
# Connections per user, used in turn: sshd allows 10 sessions per connection
# (MaxSessions), which one busy page can exceed.
_POOL_SIZE = 4


def _alive(client) -> bool:
    try:
        transport = client.get_transport()
        return transport is not None and transport.is_active() and transport.is_authenticated()
    except Exception:
        return False


def _pool_key(host, username, pkey, key_filename):
    # Keys are identified by fingerprint, so the key material is not held in the dict key.
    fp = None
    if pkey is not None:
        try:
            fp = pkey.get_fingerprint().hex()
        except Exception:
            fp = id(pkey)
    return (host, username, fp, key_filename)


def close_pooled_connections():
    """Closes every pooled client, for shutdown and tests."""
    with _POOL_LOCK:
        for entry in _POOL.values():
            for client in entry['slots']:
                if client is None:
                    continue
                try:
                    client.close()
                except Exception:
                    pass
        _POOL.clear()

_KEY_CLASSES = (
    paramiko.RSAKey,
    paramiko.Ed25519Key,
    paramiko.ECDSAKey,
)


def load_private_key(key_str: str) -> paramiko.PKey:
    """
    Auto-detect SSH key type and return a paramiko PKey object.
    Tries RSA, Ed25519, ECDSA, and DSS in order.
    Raises ValueError if none succeed.
    """
    for key_class in _KEY_CLASSES:
        try:
            return key_class.from_private_key(io.StringIO(key_str.strip()))
        except (paramiko.SSHException, Exception):
            continue
    raise ValueError('Unsupported or invalid SSH private key format')


class SSHConnection:
    """Manages paramiko SSH connections with configurable host/user/key."""

    def __init__(
        self,
        host: str | None = None,
        username: str | None = None,
        pkey: paramiko.PKey | None = None,
        key_filename: str | None = None,
    ):
        self.host = host
        self.username = username
        self.pkey = pkey               # in-memory key object (API usage)
        self.key_filename = key_filename   # file path (CLI fallback)

    def _handshake(self) -> paramiko.SSHClient:
        """Builds and authenticates a new client, outside the pool."""
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        connect_kwargs: dict = {'username': self.username}
        if self.pkey:
            connect_kwargs['pkey'] = self.pkey
        elif self.key_filename:
            connect_kwargs['key_filename'] = self.key_filename
        # With neither set, paramiko falls back to the system SSH agent.
        ssh.connect(self.host, **connect_kwargs)
        transport = ssh.get_transport()
        if transport is not None:
            transport.set_keepalive(_KEEPALIVE_SECONDS)
        LOGGER.debug('Connected to %s as %s', self.host, self.username)
        return ssh

    def connect(self, pooled: bool = True) -> paramiko.SSHClient:
        """Returns a live SSH connection, reusing a pooled one when possible.

        pooled=False returns a private client outside the pool, for callers that
        hold it beyond one request (e.g. ssh_slurm.submit_ssh_job for a whole run).
        """
        if not self.host or not self.username:
            raise ValueError(
                'SSHConnection requires host and username. '
                'Use make_user_connection() in API context, or set host/username '
                'from the user config for CLI usage.'
            )
        if not pooled:
            return self._handshake()

        key = _pool_key(self.host, self.username, self.pkey, self.key_filename)
        with _POOL_LOCK:
            entry = _POOL.setdefault(key, {'slots': [None] * _POOL_SIZE, 'next': 0,
                                           'making': [threading.Lock() for _ in range(_POOL_SIZE)]})
            slot = entry['next']
            entry['next'] = (slot + 1) % _POOL_SIZE
            client = entry['slots'][slot]
        # is_active() alone is not enough: a closed client has _transport None,
        # so anything not verifiably usable is replaced.
        if client is not None and _alive(client):
            LOGGER.debug('Reusing pooled SSH connection %d to %s', slot, self.host)
            return client
        # One handshake per slot at a time; concurrent requests wait for it.
        with entry['making'][slot]:
            with _POOL_LOCK:
                current = entry['slots'][slot]
            if current is not None and _alive(current):
                return current
            if current is not None:
                # Only dead clients are closed; a live one may still be in use.
                try:
                    current.close()
                except Exception:
                    pass
            ssh = self._handshake()
            with _POOL_LOCK:
                entry['slots'][slot] = ssh
            return ssh


def runs_here(connection: 'SSHConnection') -> bool:
    """Returns True when this process runs on the connection's cluster as the same
    user (the desktop app's API on a login node), so files can be read directly."""
    import getpass
    import socket
    if not connection.username or connection.username != getpass.getuser():
        return False
    host = (connection.host or '').lower()
    fqdn = socket.getfqdn().lower()
    return bool(host) and (fqdn == host or fqdn.endswith('.' + host))


def make_user_connection(
    cluster_host: str,
    cluster_username: str,
    private_key_str: str,
) -> SSHConnection:
    """
    Build an SSHConnection for a specific user's cluster account.

    Accepts the plaintext (already-decrypted) private key string, loads it
    into a paramiko PKey object in memory, and returns an SSHConnection.
    The key is never written to disk.
    """
    pkey = load_private_key(private_key_str)
    return SSHConnection(host=cluster_host, username=cluster_username, pkey=pkey)


def ensure_remote_dane_wf(conn: SSHConnection, *, repo_url: str = _DANE_WF_REPO_URL, timeout: float = 300.0) -> None:
    """
    Ensures the user's cluster account has a working `dane_wf` before the first job.

    Clones into ~/bioinformatics-tools and runs `uv sync` only when no built
    dane_wf exists; an existing checkout or symlink is never touched. Raises
    RuntimeError with the remote output on failure.
    """
    command = f'''set -e
if [ ! -x "$HOME/bioinformatics-tools/.venv/bin/dane_wf" ]; then
    if [ ! -d "$HOME/bioinformatics-tools" ]; then
        git clone {shlex.quote(repo_url)} "$HOME/bioinformatics-tools"
    fi
    cd "$HOME/bioinformatics-tools"
    [ -f pyproject.toml ]
    export PATH="$HOME/.local/bin:$PATH"
    if ! command -v uv >/dev/null 2>&1; then
        curl -LsSf https://astral.sh/uv/install.sh | sh >/dev/null 2>&1 || true
    fi
    command -v uv >/dev/null 2>&1
    uv sync
fi
test -x "$HOME/bioinformatics-tools/.venv/bin/dane_wf"
'''
    ssh = conn.connect()
    _, stdout, stderr = ssh.exec_command(command, timeout=timeout)
    stdout.channel.settimeout(timeout)
    exit_code = stdout.channel.recv_exit_status()
    output = (stdout.read().decode(errors='replace') + stderr.read().decode(errors='replace')).strip()
    if exit_code != 0:
        raise RuntimeError(
            f'Could not provision dane_wf on the remote account (exit {exit_code}): {output or "(no output)"}'
        )


def sync_remote_dane_wf(conn: SSHConnection, *, ref: str = _MARGIE_SB_REF, timeout: float = 30.0) -> str:
    """
    Fast-forwards ~/bioinformatics-tools to `ref` on origin when it is behind and clean.

    Cheap in the common case (a `git fetch` plus a SHA comparison); the
    `git reset --hard` + `uv sync` runs only when a new deployment has landed.
    Never raises. Returns 'unprovisioned', 'not-a-git-checkout', 'up-to-date',
    'dirty-skipped', 'updated', 'disabled' or 'error: <detail>'.

    BSP_SKIP_DANE_WF_SYNC=1 disables it, for a developer whose API and cluster
    are one machine (scripts/dev-local/start.sh sets it).
    """
    if os.environ.get('BSP_SKIP_DANE_WF_SYNC'):
        LOGGER.info('dane_wf version-sync check disabled by BSP_SKIP_DANE_WF_SYNC')
        return 'disabled'
    command = f'''cd "$HOME/bioinformatics-tools" 2>/dev/null || {{ echo NO_CHECKOUT; exit 0; }}
[ -d .git ] || {{ echo NOT_GIT; exit 0; }}
git fetch origin {shlex.quote(ref)} --quiet 2>&1
local_sha=$(git rev-parse HEAD 2>/dev/null)
remote_sha=$(git rev-parse origin/{shlex.quote(ref)} 2>/dev/null)
if [ -z "$remote_sha" ] || [ "$local_sha" = "$remote_sha" ]; then
    echo UP_TO_DATE
    exit 0
fi
if [ -n "$(git status --porcelain)" ]; then
    echo DIRTY_SKIPPED
    exit 0
fi
git checkout {shlex.quote(ref)} --quiet 2>/dev/null || true
git reset --hard "$remote_sha" --quiet
export PATH="$HOME/.local/bin:$PATH"
uv sync --quiet
echo "UPDATED $remote_sha"
'''
    try:
        ssh = conn.connect()
        _, stdout, stderr = ssh.exec_command(command, timeout=timeout)
        stdout.channel.settimeout(timeout)
        exit_code = stdout.channel.recv_exit_status()
        output = (stdout.read().decode(errors='replace') + stderr.read().decode(errors='replace')).strip()
    except Exception as exc:
        LOGGER.warning('dane_wf version-sync check raised, skipping: %s', exc)
        return f'error: {exc}'

    if exit_code != 0:
        LOGGER.warning('dane_wf version-sync check failed (exit %s): %s', exit_code, output)
        return f'error: exit {exit_code}: {output}'

    if 'NO_CHECKOUT' in output:
        return 'unprovisioned'
    if 'NOT_GIT' in output:
        return 'not-a-git-checkout'
    if 'DIRTY_SKIPPED' in output:
        LOGGER.info('Remote ~/bioinformatics-tools has local changes; skipped auto-update.')
        return 'dirty-skipped'
    if 'UPDATED' in output:
        LOGGER.info('Updated remote ~/bioinformatics-tools: %s', output)
        return 'updated'
    return 'up-to-date'


