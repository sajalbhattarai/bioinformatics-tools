"""
Genome and protein identity hashes for MARGIE's caches (stdlib hashlib).

genome_hash() covers each record's id and upper-case bases in file order, ignoring
formatting; ids and order count because RASTtk numbering depends on them.
legacy_file_hash() is the file-byte hash still used to find older cache entries.
"""
from __future__ import annotations

import gzip
import hashlib
from pathlib import Path
from typing import IO, Iterator

# Marks a sequence hash, so it is never mistaken for a byte hash of a file.
GENOME_HASH_PREFIX = 'fa1:'


def _open_text(path: str | Path) -> IO[str]:
    p = str(path)
    if p.endswith('.gz'):
        return gzip.open(p, 'rt', encoding='utf-8', errors='replace', newline=None)
    return open(p, 'r', encoding='utf-8', errors='replace', newline=None)


def genome_hash(path: str | Path) -> str:
    """Returns the SHA-256 of the records' ids and bases, independent of file formatting."""
    sha = hashlib.sha256()
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith('>'):
                words = line[1:].split()
                sha.update(b'>' + (words[0] if words else '').encode() + b'\n')
            else:
                seq = ''.join(line.split()).upper()
                if seq:
                    sha.update(seq.encode())
    return GENOME_HASH_PREFIX + sha.hexdigest()


def legacy_file_hash(path: str | Path) -> str:
    """Returns the SHA-256 of the file's bytes (the older genome identity)."""
    sha = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            sha.update(chunk)
    return sha.hexdigest()


def genome_hashes(path: str | Path) -> list[str]:
    """Returns [genome_hash, legacy_file_hash]; lookups try both, new records use the first."""
    return [genome_hash(path), legacy_file_hash(path)]


def normalise_protein(seq: str) -> str:
    return ''.join(seq.split()).upper().rstrip('*')


def protein_hash(seq: str) -> str:
    return hashlib.sha256(normalise_protein(seq).encode()).hexdigest()


def read_fasta(path: str | Path) -> Iterator[tuple[str, str, str]]:
    """Yields (id, header without '>', sequence) for each FASTA record."""
    header = None
    parts: list[str] = []
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if header is not None:
                    words = header.split()
                    yield (words[0] if words else ''), header, ''.join(parts)
                header = line[1:].rstrip('\n')
                parts = []
            else:
                parts.append(''.join(line.split()))
    if header is not None:
        words = header.split()
        yield (words[0] if words else ''), header, ''.join(parts)
