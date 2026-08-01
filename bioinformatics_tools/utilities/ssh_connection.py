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
import threading
import time

import paramiko

LOGGER = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Connection pool.
#
# connect() used to build a fresh paramiko.SSHClient and complete a full TCP +
# SSH + public-key handshake on EVERY call. Every file listing, history load
# and status poll in the GUI paid that -- typically several hundred ms to well
# over a second against an HPC login node -- before doing any actual work. That
# is why the file and history lists felt slow: almost all of the wait was
# reconnecting, not listing.
#
# Clients are now reused per (host, username, key). A pooled client is handed
# back only if its transport is still active; a dead one is discarded and
# replaced, so a dropped VPN or a bounced login node self-heals on the next
# request rather than raising.
_POOL: dict = {}
_POOL_LOCK = threading.Lock()
# Long enough to cover a browsing session, short enough that a stale client is
# not kept around indefinitely.
_POOL_TTL = 600.0


def _pool_key(host, username, pkey, key_filename):
    # Keys are objects; their fingerprint identifies the credential without
    # holding the material in the dict key.
    fp = None
    if pkey is not None:
        try:
            fp = pkey.get_fingerprint().hex()
        except Exception:
            fp = id(pkey)
    return (host, username, fp, key_filename)


def close_pooled_connections():
    """Drop every pooled client. For shutdown and tests."""
    with _POOL_LOCK:
        for client, _ in _POOL.values():
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

    def connect(self) -> paramiko.SSHClient:
        """Return a live SSH connection, reusing a pooled one when possible."""
        if not self.host or not self.username:
            raise ValueError(
                'SSHConnection requires host and username. '
                'Use make_user_connection() in API context, or set host/username '
                'from the user config for CLI usage.'
            )
        key = _pool_key(self.host, self.username, self.pkey, self.key_filename)
        now = time.time()
        with _POOL_LOCK:
            hit = _POOL.get(key)
            if hit:
                client, created = hit
                transport = client.get_transport()
                if transport is not None and transport.is_active() and now - created < _POOL_TTL:
                    LOGGER.debug('Reusing pooled SSH connection to %s', self.host)
                    return client
                # Dead or expired: bin it and fall through to reconnect, so a
                # dropped connection heals silently instead of erroring.
                _POOL.pop(key, None)
                try:
                    client.close()
                except Exception:
                    pass

        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        connect_kwargs: dict = {'username': self.username}
        if self.pkey:
            connect_kwargs['pkey'] = self.pkey
        elif self.key_filename:
            connect_kwargs['key_filename'] = self.key_filename
        # If neither is set, paramiko falls back to the system SSH agent (CLI default)
        ssh.connect(self.host, **connect_kwargs)
        LOGGER.debug('Connected to %s as %s', self.host, self.username)
        with _POOL_LOCK:
            _POOL[key] = (ssh, time.time())
        return ssh


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


