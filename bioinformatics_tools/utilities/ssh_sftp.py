"""
SFTP file operations over SSH.

Provides remote directory listing, file streaming, and YAML config
read/write via paramiko SFTP, plus a paginated line-range reader (via
SSH exec, not SFTP) for the file viewer.

All functions are API-layer only. Pass a per-user SSHConnection built
with make_user_connection() for every call.
"""
import logging
import re
import shlex
import stat
import time

import yaml

from bioinformatics_tools.utilities.ssh_connection import SSHConnection

LOGGER = logging.getLogger(__name__)

_PAGE_SENTINEL = "___MARGIE_PAGE_SENTINEL___"


def list_remote_dir(
    remote_path: str,
    connection: SSHConnection,
) -> list[dict]:
    """List files and directories in a remote path via SFTP.

    Returns a list of dicts: {name, type, size, mtime}; mtime is seconds
    since the epoch, or None when the server does not report it.
    """
    ssh = connection.connect()
    sftp = ssh.open_sftp()
    entries = []
    for attr in sftp.listdir_attr(remote_path):
        entry_type = 'directory' if stat.S_ISDIR(attr.st_mode) else 'file'
        entries.append({
            'name': attr.filename,
            'type': entry_type,
            'size': attr.st_size,
            'mtime': attr.st_mtime,
        })
    sftp.close()
    pass  # pooled client: closing it would defeat SSHConnection's pool
    return entries


def list_remote_dir_checked(
    remote_path: str,
    connection: SSHConnection,
) -> list[dict]:
    """Lists a remote directory, checking it is a directory in the same SFTP session.

    Raises FileNotFoundError when the path is missing and NotADirectoryError when it is a file."""
    ssh = connection.connect()
    sftp = ssh.open_sftp()
    try:
        try:
            attr = sftp.stat(remote_path)
        except FileNotFoundError:
            raise FileNotFoundError(f'Path not found on cluster: {remote_path}')
        if not stat.S_ISDIR(attr.st_mode):
            raise NotADirectoryError(remote_path)
        return [{
            'name': a.filename,
            'type': 'directory' if stat.S_ISDIR(a.st_mode) else 'file',
            'size': a.st_size,
            'mtime': a.st_mtime,
        } for a in sftp.listdir_attr(remote_path)]
    finally:
        sftp.close()


def stream_remote_file(
    remote_path: str,
    connection: SSHConnection,
):
    """Opens a remote file eagerly and returns a generator that streams it in 8KB chunks.

    The open happens at call time, so callers can catch FileNotFoundError /
    IOError before the StreamingResponse sends its headers.
    """
    ssh = connection.connect()
    sftp = ssh.open_sftp()
    f = sftp.open(remote_path, 'rb')   # raises FileNotFoundError/IOError here if absent

    def _chunks():
        try:
            while True:
                chunk = f.read(8192)
                if not chunk:
                    break
                yield chunk
        finally:
            f.close()
            sftp.close()
            pass  # pooled client: closing it would defeat SSHConnection's pool

    return _chunks()


def read_remote_yaml(
    remote_path: str,
    connection: SSHConnection,
) -> dict:
    """Read and parse a YAML file from the remote cluster.

    Returns the parsed dict, or an empty dict if the file does not exist.
    """
    ssh = connection.connect()
    sftp = ssh.open_sftp()
    try:
        with sftp.open(remote_path, 'r') as f:
            content = f.read().decode('utf-8')
        return yaml.safe_load(content) or {}
    except FileNotFoundError:
        LOGGER.info('Remote config not found at %s — returning empty dict', remote_path)
        return {}
    finally:
        sftp.close()
        pass  # pooled client: closing it would defeat SSHConnection's pool


def check_remote_file(
    path: str,
    connection: SSHConnection,
) -> None:
    """Verify a remote file exists and is a regular file via SFTP.

    Raises FileNotFoundError if the path does not exist on the cluster,
    or IsADirectoryError if it resolves to a directory rather than a file.
    """
    ssh = connection.connect()
    sftp = ssh.open_sftp()
    try:
        attr = sftp.stat(path)
        if stat.S_ISDIR(attr.st_mode):
            raise IsADirectoryError(f'Path is a directory, not a file: {path}')
    except FileNotFoundError:
        raise FileNotFoundError(f'File not found on cluster: {path}')
    finally:
        sftp.close()
        pass  # pooled client: closing it would defeat SSHConnection's pool


def stat_remote_file(path: str, connection: SSHConnection) -> tuple[float, int]:
    """Returns (mtime, size) for a remote file from one SFTP stat, used as a cache key.

    Raises FileNotFoundError if the path does not exist on the cluster.
    """
    ssh = connection.connect()
    sftp = ssh.open_sftp()
    try:
        attr = sftp.stat(path)
    except FileNotFoundError:
        raise FileNotFoundError(f'File not found on cluster: {path}')
    finally:
        sftp.close()
        pass  # pooled client: closing it would defeat SSHConnection's pool
    return (attr.st_mtime, attr.st_size)


def check_remote_path_kind(path: str, connection: SSHConnection) -> str:
    """Returns 'file' or 'directory' for a remote path.

    Raises FileNotFoundError if the path does not exist on the cluster.
    """
    ssh = connection.connect()
    sftp = ssh.open_sftp()
    try:
        attr = sftp.stat(path)
    except FileNotFoundError:
        raise FileNotFoundError(f'Path not found on cluster: {path}')
    finally:
        sftp.close()
        pass  # pooled client: closing it would defeat SSHConnection's pool
    return 'directory' if stat.S_ISDIR(attr.st_mode) else 'file'


def write_remote_yaml(
    remote_path: str,
    data: dict,
    connection: SSHConnection,
) -> None:
    """Write a dict as YAML to a remote path via SFTP.

    Creates parent directories on the remote if they do not exist.
    """
    ssh = connection.connect()

    # Ensure parent directory exists
    parent = remote_path.rsplit('/', 1)[0]
    if parent:
        ssh.exec_command(f'mkdir -p {parent}')

    sftp = ssh.open_sftp()
    content = yaml.dump(data, default_flow_style=False, allow_unicode=True)
    with sftp.open(remote_path, 'w') as f:
        f.write(content)
    sftp.close()
    pass  # pooled client: closing it would defeat SSHConnection's pool
    LOGGER.info('Wrote remote config to %s', remote_path)


def write_remote_text_file(
    remote_path: str,
    content: str,
    connection: SSHConnection,
) -> None:
    """Writes raw text to a remote path via SFTP, creating parent directories.
    Used by the file explorer's Save action.
    """
    ssh = connection.connect()

    parent = remote_path.rsplit('/', 1)[0]
    if parent:
        ssh.exec_command(f'mkdir -p {shlex.quote(parent)}')

    sftp = ssh.open_sftp()
    try:
        with sftp.open(remote_path, 'w') as f:
            f.write(content)
    finally:
        sftp.close()
        pass  # pooled client: closing it would defeat SSHConnection's pool
    LOGGER.info('Wrote remote text file to %s', remote_path)


def copy_remote_directory(
    src_path: str,
    dest_path: str,
    connection: SSHConnection,
) -> None:
    """Copies a remote directory tree to a new path with one rsync on the cluster.

    Excludes .snakemake/ (it embeds the old absolute path and is regenerated)
    and original_container_outputs/*/stage/ (per-run scratch). Creates dest_path
    if needed; raises RuntimeError if rsync exits non-zero.
    """
    ssh = connection.connect()
    try:
        quoted_src = shlex.quote(src_path.rstrip('/') + '/')
        quoted_dest = shlex.quote(dest_path)
        cmd = (
            f'mkdir -p {quoted_dest} && '
            f'rsync -a '
            f"--exclude='.snakemake' "
            f"--exclude='original_container_outputs/*/stage' "
            f'{quoted_src} {quoted_dest}'
        )
        LOGGER.info('Copying remote directory %s -> %s', src_path, dest_path)
        _, stdout, stderr = ssh.exec_command(cmd)
        exit_code = stdout.channel.recv_exit_status()
        if exit_code != 0:
            err = stderr.read().decode('utf-8', errors='replace')
            raise RuntimeError(f'rsync failed (exit {exit_code}) copying {src_path} to {dest_path}: {err}')
    finally:
        pass  # pooled client: closing it would defeat SSHConnection's pool


def stage_selected_genomes(
    source_dir: str,
    names: list[str],
    connection: SSHConnection,
    label: str = '',
) -> str:
    """Builds a folder of symlinks to the chosen genomes and returns its absolute path.

    A workflow annotates everything in its input folder, so a subset gets its
    own folder under the user's home; old selections are kept as a record.

    Raises ValueError for a name that is not a plain file name, and
    RuntimeError if the remote command fails.
    """
    if not names:
        raise ValueError('No genomes were selected.')
    for name in names:
        # Names come from the browser; anything with a slash could escape source_dir.
        if not name or '/' in name or name in ('.', '..'):
            raise ValueError(f'Not a genome file name: {name!r}')

    stamp = time.strftime('%Y%m%d-%H%M%S')
    suffix = f"-{re.sub(r'[^A-Za-z0-9]+', '-', label).strip('-')}" if label else ''
    rel = f'.margie/selections/{stamp}{suffix}'

    src = shlex.quote(source_dir.rstrip('/'))
    links = ' && '.join(
        f'ln -sfn {src}/{shlex.quote(n)} "$d"/{shlex.quote(n)}' for n in names
    )
    # printf returns the absolute path without a trailing newline.
    cmd = f'd="$HOME"/{shlex.quote(rel)} && mkdir -p "$d" && {links} && printf %s "$d"'

    ssh = connection.connect()
    LOGGER.info('Staging %d selected genome(s) from %s', len(names), source_dir)
    _, stdout, stderr = ssh.exec_command(cmd)
    out = stdout.read().decode('utf-8', errors='replace').strip()
    exit_code = stdout.channel.recv_exit_status()
    if exit_code != 0 or not out:
        err = stderr.read().decode('utf-8', errors='replace')
        raise RuntimeError(
            f'Could not stage the selected genomes on the cluster (exit {exit_code}): {err or "no path returned"}'
        )
    return out


def _build_path_rewrite_script(directory: str, old_path: str, new_path: str) -> str:
    """Builds the Python source that rewrite_path_references() runs remotely.

    Paths are embedded via repr() and replaced with a literal str.replace(), so
    no escaping is needed. Kept separate so tests can run the script locally.
    """
    return (
        "import os\n"
        f"old = {old_path!r}\n"
        f"new = {new_path!r}\n"
        "modified = 0\n"
        f"for root, _, files in os.walk({directory!r}):\n"
        "    for name in files:\n"
        "        path = os.path.join(root, name)\n"
        "        try:\n"
        "            with open(path, 'r', encoding='utf-8') as f:\n"
        "                content = f.read()\n"
        "        except (UnicodeDecodeError, OSError):\n"
        "            continue\n"
        "        if old not in content:\n"
        "            continue\n"
        "        try:\n"
        "            with open(path, 'w', encoding='utf-8') as f:\n"
        "                f.write(content.replace(old, new))\n"
        "            modified += 1\n"
        "        except OSError:\n"
        "            continue\n"
        "print(modified)\n"
    )


def rewrite_path_references(
    directory: str,
    old_path: str,
    new_path: str,
    connection: SSHConnection,
) -> int:
    """Replaces every literal old_path with new_path in text files under directory.

    Used after copy_remote_directory() on Resume to fix provenance paths some
    tools write into their outputs. Runs as an embedded Python script over SSH
    to avoid sed escaping; binary files and per-file errors are skipped.

    Returns the number of files modified.
    """
    ssh = connection.connect()
    try:
        script = _build_path_rewrite_script(directory, old_path, new_path)
        remote_python = '~/bioinformatics-tools/.venv/bin/python'
        cmd = f'{remote_python} -c {shlex.quote(script)}'
        _, stdout, stderr = ssh.exec_command(cmd)
        out = stdout.read().decode().strip()
        exit_code = stdout.channel.recv_exit_status()
        if exit_code != 0:
            err = stderr.read().decode('utf-8', errors='replace')
            LOGGER.warning('rewrite_path_references failed (exit %d) for %s: %s', exit_code, directory, err)
            return 0
        return int(out) if out.isdigit() else 0
    finally:
        pass  # pooled client: closing it would defeat SSHConnection's pool


def read_remote_file_page(
    remote_path: str,
    start_row: int,
    end_row: int,
    connection: SSHConnection,
    known_total_lines: int | None = None,
) -> dict:
    """Reads lines [start_row, end_row] (1-indexed, inclusive) of a remote text
    file, plus its header line and total line count, in one SSH exec_command.

    sed quits after end_row, so cost scales with page position, not file size.
    known_total_lines skips the wc -l pass.

    Returns {"total_lines": int, "header": str, "lines": list[str]}.
    Raises FileNotFoundError if remote_path does not exist, or
    IsADirectoryError if it resolves to a directory.
    """
    check_remote_file(remote_path, connection)

    quoted = shlex.quote(remote_path)
    parts = []
    if known_total_lines is None:
        parts.append(f'wc -l < {quoted}')
    # `1{p;q}` prints line 1 and quits.
    parts.append(f"sed -n '1{{p;q}}' {quoted}")
    # `p` must precede `q` on end_row, or the page's last row is dropped.
    parts.append(f"sed -n '{start_row},{end_row}p;{end_row}q' {quoted}")
    script = f'; echo {_PAGE_SENTINEL}; '.join(parts)

    ssh = connection.connect()
    try:
        _, stdout, stderr = ssh.exec_command(script)
        output = stdout.read().decode('utf-8', errors='replace')
        err = stderr.read().decode('utf-8', errors='replace').strip()
        if err:
            LOGGER.warning('read_remote_file_page stderr for %s: %s', remote_path, err)
    finally:
        pass  # pooled client: closing it would defeat SSHConnection's pool

    sections = output.split(f'{_PAGE_SENTINEL}\n')
    if known_total_lines is None:
        total_lines = int(sections[0].strip() or 0)
        header = sections[1].rstrip('\n')
        data_block = sections[2]
    else:
        total_lines = known_total_lines
        header = sections[0].rstrip('\n')
        data_block = sections[1]

    lines = data_block.split('\n')
    if lines and lines[-1] == '':
        lines.pop()

    return {"total_lines": total_lines, "header": header, "lines": lines}
