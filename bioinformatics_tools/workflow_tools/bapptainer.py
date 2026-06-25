"""
All code related to fetching containers
"""
import hashlib
import json
import logging
import shutil
import subprocess
import sys
from pathlib import Path

import requests
from tqdm import tqdm

from bioinformatics_tools.workflow_tools.models import ApptainerKey

LOGGER = logging.getLogger(__name__)


class CacheSifError(Exception):
    """Raised when there's an error caching SIF files"""
    pass


class ApptainerRunError(Exception):
    """Raised when a container run (via its baked-in runscript) exits non-zero"""
    pass

CACHE_DIR = Path.home() / ".cache" / "bioinformatics-tools"


def _resolve_cache_dir(local_sif_dir: str | Path | None = None) -> Path:
    """Return the configured cache directory or the default cache location."""
    if local_sif_dir:
        return Path(local_sif_dir).expanduser()
    return CACHE_DIR


def verify_sha256(file_path: Path, expected_sha256: str) -> bool:
    """Verify file SHA256 checksum"""
    sha256 = hashlib.sha256()
    with open(file_path, 'rb') as f:
        for chunk in iter(lambda: f.read(8192), b''):
            sha256.update(chunk)
    return sha256.hexdigest() == expected_sha256


def download_container_with_progress(url: Path, dest: Path, expected_sha256: str | None = None):
    '''Download with progress bar and compare against a SHA256 if available'''
    registry_url = 'https://github.com/wintermutant/biotools-containers/releases/download'
    latest_version ='v0.0.2'
    full_url = f"{registry_url}/{latest_version}/{url}"
    LOGGER.info('Downloading data from the registry...%s', full_url)
    response = requests.get(full_url, stream=True, timeout=None)
    response.raise_for_status()

    total_size = int(response.headers.get('content-length', 0))

    with open(dest, 'wb') as f, tqdm(
        desc=dest.name, total=total_size,
        unit="B", unit_scale=True, unit_divisor=1024
    ) as progress_bar:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
            progress_bar.update(len(chunk))
    
    if expected_sha256:
        print('Verifying SHA256...')
        if verify_sha256(dest, expected_sha256):
            print('Verification successful')
        else:
            dest.unlink()  # Delete the corrupted file
            raise ValueError('SHA256 verification failed')
    
    return dest


def find_apptainer_command(apptainer_path: str | None = None) -> str | None:
    '''Check for system-level apptainer command'''
    commands_to_check = [apptainer_path, 'apptainer.lima', 'apptainer']
    for cmd in commands_to_check:
        if cmd and shutil.which(cmd):
            LOGGER.debug('Running command: %s', cmd)
            return cmd
    return None


def init_cache(local_sif_dir: str | Path | None = None):
    '''standard or custom cache location'''
    _resolve_cache_dir(local_sif_dir).mkdir(parents=True, exist_ok=True)


def get_cached_file(filename: Path, local_sif_dir: str | Path | None = None) -> Path | None:
    '''Search default cache directory for a file'''
    # TODO: Later, put the version in the cached name
    cached_file = _resolve_cache_dir(local_sif_dir) / filename
    if cached_file.exists():
        print('Found file in cache. Using that')
        return cached_file
    return None


def init_apptainer():
    pass


def run_apptainer():
    pass


def pull_container_from_ghcr(container_name: Path, tag: str, output_filename: str | None = None,
                             local_sif_dir: str | Path | None = None) -> Path:
    '''Pull container from GitHub Container Registry using apptainer

    Args:
        container_name: Name of the container (e.g., 'prodigal')
        tag: Version tag (e.g., '2.6.3-v1.0')
        output_filename: Optional output filename (e.g., 'prodigal.sif').
                        If None, defaults to '{container_name}.sif'

    Returns:
        Path to the pulled .sif file in cache directory

    Example:
        pull_container_from_ghcr('prodigal', '2.6.3-v1.0')
        # Pulls docker://ghcr.io/wintermutant/prodigal:2.6.3-v1.0
        # Saves to ~/.cache/bioinformatics-tools/prodigal.sif
    '''
    # Determine output filename
    if output_filename is None:
        output_filename = f'{container_name}.sif'

    # Determine docker URL early so both paths can use it
    str_container_name = str(container_name)
    docker_url = f'docker://ghcr.io/wintermutant/{str_container_name}:{tag}'

    # Check if already cached
    cached = get_cached_file(Path(output_filename), local_sif_dir=local_sif_dir)
    if cached:
        LOGGER.info('Found %s in cache', output_filename)
        _emit_container_metadata(str_container_name, tag, str(cached), "cached", docker_url)
        return cached

    dest = _resolve_cache_dir(local_sif_dir) / output_filename

    LOGGER.info('Pulling %s from GitHub Container Registry...', docker_url)

    apptainer_cmd = find_apptainer_command()
    if not apptainer_cmd:
        raise RuntimeError('Apptainer not found. Cannot pull container.')

    # Ensure destination directory exists
    dest.parent.mkdir(parents=True, exist_ok=True)

    # Build apptainer pull command
    cmd = [apptainer_cmd, 'pull', str(dest), docker_url]
    LOGGER.info('Running: %s', ' '.join(cmd))

    try:
        subprocess.run(cmd, check=True)
        LOGGER.info('Successfully pulled %s to %s', container_name, dest)
        _emit_container_metadata(str_container_name, tag, str(dest), "downloaded", docker_url)
        return dest
    except subprocess.CalledProcessError as e:
        LOGGER.error('Failed to pull container: %s', e)
        raise


def get_sif_files(sif_paths: list[Path]):
    pass


def cache_sif_files(sif_paths: list[tuple[str, str]], local_sif_dir: str | Path | None = None):
    '''ensure all paths are in ~/.cache and accessible

    Raises:
        CacheSifError: If any SIF file cannot be cached
    '''
    for sif_name, sif_version in sif_paths:
        try:
            get_verified_sif_file(sif_name, sif_version, local_sif_dir=local_sif_dir)
        except Exception as e:
            LOGGER.critical('Issue getting a cached file: %s:%s', sif_name, sif_version)
            raise CacheSifError(f'Failed to cache {sif_name}:{sif_version}') from e


def get_verified_sif_file(sif_name: str, sif_version: str,
                          local_sif_dir: str | Path | None = None):
    '''search cache for sif file, if does not exist then download'''
    # sif_path = [(prodigal, 2.6.3), (program, version)]
    sif_path = Path(sif_name)
    LOGGER.info('SIF path: %s Version: %s', sif_path, sif_version)
    docker_url = f'docker://ghcr.io/wintermutant/{sif_path}:{sif_version}'  # TODO: Add a list of container repos

    cached = get_cached_file(sif_path, local_sif_dir=local_sif_dir)
    if cached:
        source = "local" if local_sif_dir else "cached"
        _emit_container_metadata(sif_name, sif_version, str(cached), source, docker_url)
        return cached
    dest = _resolve_cache_dir(local_sif_dir) / sif_path
    LOGGER.info('Destination: %s', dest)
    #TODO: un-hardcode this
    verified_sif_file = pull_container_from_ghcr(sif_path, sif_version, local_sif_dir=local_sif_dir)
    LOGGER.info('Verified sif: %s', verified_sif_file)
    return verified_sif_file


def locate_local_sif_files(sif_paths: list[tuple[str, str]], local_sif_dir: str | Path | None = None):
    '''Resolve each SIF file against the local filesystem only — never contacts
    the container registry. Emits container metadata for every file found so
    the job's container table shows exactly where it was loaded from.'''
    cache_dir = _resolve_cache_dir(local_sif_dir)
    for sif_name, sif_version in sif_paths:
        resolved = cache_dir / sif_name
        if resolved.exists():
            _emit_container_metadata(sif_name, sif_version, str(resolved), "local", "")
        else:
            LOGGER.warning('SIF file not found locally (no registry pull attempted): %s', resolved)


def _emit_container_metadata(name: str, version: str, path: str, source: str, docker_url: str):
    """Log structured container metadata for the task runner to parse."""
    metadata = {
        "name": name,
        "version": version,
        "path": path,
        "resolved_path": path,
        "source": source,
        "registry_url": docker_url,
        "docker_url": docker_url,
    }
    LOGGER.info('__CONTAINER__:%s', json.dumps(metadata))


def run_apptainer_container(app_obj: ApptainerKey, container_args: list[str]) -> int:
    """
    Run an Apptainer container with the specified arguments.
    Returns exit code from the container execution
    """
    sif_path = Path(app_obj.sif_path)
    LOGGER.debug('Sif path: %s', sif_path)

    if not sif_path:
        sys.exit('No bueno, brotha')

    verified_sif_file = get_verified_sif_file(sif_path)

    apptainer_command = find_apptainer_command(app_obj.executable)
    if apptainer_command is None:
        LOGGER.error('Apptainer not found. Please install Apptainer LiMa')
        return 127

    # Build the command
    cmd = [apptainer_command, 'exec', str(verified_sif_file)] + container_args
    cmd_string = ' '.join(cmd)
    LOGGER.info(cmd_string)

    # Execute the container
    try:
        result = subprocess.run(
            cmd,
            capture_output=False,  # Stream output directly to terminal
            text=True,
            check=True
        )
        return result.returncode
    except FileNotFoundError:
        LOGGER.error("Apptainer not found. Please install Apptainer/Singularity.")
        LOGGER.error("See: https://apptainer.org/docs/admin/main/installation.html")
        return 127
    except Exception as e:
        LOGGER.error("Failed to run container: %s", e)
        return 1


def run_apptainer_run(sif_path: str | Path, args: list[str], binds: list[tuple[str, str, str]],
                       log_path: str | Path | None = None, env: dict[str, str] | None = None) -> None:
    """Run a margie_sb container via its baked-in runscript (`apptainer run`, not
    `exec`) — every margie_sb entrypoint.sh is installed as the container's
    runscript and documents its own `apptainer run <sif> ...` usage.

    Args:
        sif_path: resolved, already-existing local .sif path.
        args: CLI args passed straight to the container's entrypoint (e.g. ['-i', '/input', '-o', '/output']).
        binds: list of (host_path, container_path, mode) tuples, mode is 'ro' or 'rw'.
        log_path: if given, stdout+stderr are teed to this file instead of streamed only.
        env: optional environment variables to set inside the container (e.g. a
            few entrypoints, like gtdbtk's, branch on env vars rather than flags).

    Raises:
        ApptainerRunError: apptainer is missing, the sif is missing, a 'ro' bind
            source does not exist, or the container exits non-zero.
    """
    apptainer_command = find_apptainer_command()
    if apptainer_command is None:
        raise ApptainerRunError('Apptainer not found. Please install Apptainer/Singularity.')

    resolved_sif = Path(sif_path)
    if not resolved_sif.exists():
        raise ApptainerRunError(f'SIF file not found: {resolved_sif}')

    bind_args = []
    for host_path, container_path, mode in dict.fromkeys(binds):  # de-dupe, preserve order
        host = Path(host_path)
        if mode == 'ro' and not host.exists():
            raise ApptainerRunError(f'Bind source does not exist: {host} (-> {container_path})')
        host.mkdir(parents=True, exist_ok=True)
        bind_args += ['-B', f'{host}:{container_path}:{mode}']

    env_args = []
    for key, value in (env or {}).items():
        env_args += ['--env', f'{key}={value}']

    cmd = [apptainer_command, 'run'] + bind_args + env_args + [str(resolved_sif)] + args
    LOGGER.info('Running: %s', ' '.join(cmd))

    if log_path:
        resolved_log = Path(log_path)
        resolved_log.parent.mkdir(parents=True, exist_ok=True)
        with open(resolved_log, 'w') as log_file:
            result = subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT, text=True)
    else:
        resolved_log = None
        result = subprocess.run(cmd, text=True)

    if result.returncode != 0:
        detail = f' (see {resolved_log})' if resolved_log else ''
        raise ApptainerRunError(
            f'apptainer run failed for {resolved_sif.name} (exit {result.returncode}){detail}'
        )
