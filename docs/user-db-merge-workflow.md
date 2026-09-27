# User DB Merge Workflow

This document defines the admin workflow for merging per-user writable
SQLite databases into a new shared version while preventing duplicate
organisms.

## Objectives

- Keep runtime writes user-specific to avoid shared-file lock conflicts.
- Produce immutable shared versions (`shared-vN.db` -> `shared-vN+1.db`).
- Deduplicate organism rows by sequence identity hash.
- Keep merge provenance/audit for reproducibility.

## Versioning

- Shared base: `shared-vN.db`
- User branch DBs: `username-shared-vN-vK.db`
- Merge output: `shared-vN+1.db`
- Latest alias: `shared-current.db` (updated only after validation)

## Merge Inputs

- One base shared DB version.
- One or more user DB versions derived from that base.
- A dedupe key (`fasta_hash` preferred).

## Deduplication Rules

1. Organism identity is determined by canonical `fasta_hash`.
2. If hash does not exist in target: insert.
3. If hash exists: merge metadata deterministically.
4. For critical-field conflicts: write to conflict audit, do not overwrite silently.

## Safe Merge Procedure

1. Copy `shared-vN.db` to new candidate `shared-vN+1.db`.
2. Merge each user DB into the candidate in deterministic order.
3. Record audit counts:
   - inserted
   - duplicate-skipped
   - metadata-updated
   - conflicts
4. Run integrity checks (`PRAGMA integrity_check`, unique-key checks).
5. Publish candidate by switching alias `shared-current.db`.

## Tooling Scaffold

A planning/validation script exists at:

- `backend/scripts/merge_user_databases.py`

Current status:

- Validates input files and dedupe key.
- Reports planned merge steps as JSON.
- Execution mode is intentionally not implemented yet.

## Suggested Next Implementation Steps

1. Add table-specific upsert SQL for each organism-bearing table.
2. Add an explicit merge audit table in the output shared DB.
3. Add a dry-run diff summary (how many rows would be inserted/updated/conflicted).
4. Add CI test fixtures with overlapping hashes across source DBs.
