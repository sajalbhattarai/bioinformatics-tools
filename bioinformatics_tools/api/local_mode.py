"""
Local mode: the API runs on the user's own computer instead of an HPC cluster.

The GUI's local option starts the API with BSP_LOCAL_MODE=1 on 127.0.0.1;
a hosted deployment never sets it and keeps requiring cluster credentials.
"""
import getpass
import os

# Stored as cluster_host for local accounts (the users table's cluster columns are NOT NULL).
LOCAL_CLUSTER_HOST = 'localhost'


def is_local_mode() -> bool:
    return os.getenv('BSP_LOCAL_MODE') == '1'


def local_account() -> dict:
    """Returns the cluster-side fields for a local account: this machine and user."""
    return {
        'cluster_host': LOCAL_CLUSTER_HOST,
        'cluster_username': getpass.getuser(),
        'home_dir': os.path.expanduser('~'),
    }
