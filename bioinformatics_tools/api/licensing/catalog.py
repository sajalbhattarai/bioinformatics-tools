"""Load licensing terms + catalog, record acceptances, and check acceptance.

Data files (same directory):
  - terms.md               versioned acceptance terms (first line carries the version)
  - licensing_catalog.json per-tool license metadata + provenance

Environment overrides:
  - LICENSE_RECORDS_DIR   base dir for the signed record copies
                          (default: /depot/lindems/data/margie/licensing-records)
  - MARGIE_OPERATOR       operator identifier used in the record path (default: lindems)
"""
from __future__ import annotations

import hashlib
import json
import os
import re
from collections.abc import Iterable
from datetime import datetime, timezone
from pathlib import Path

from bioinformatics_tools.api.database import get_db

_HERE = Path(__file__).resolve().parent
_TERMS_PATH = _HERE / "terms.md"
_CATALOG_PATH = _HERE / "licensing_catalog.json"

_DEFAULT_RECORDS_DIR = "/depot/lindems/data/margie/licensing-records"
_VERSION_RE = re.compile(r"terms_version:\s*([0-9A-Za-z._-]+)")

# The acknowledgments the user must check. Keep ids stable — they are stored in
# each acceptance record. Text mirrors terms.md ("Your acknowledgments").
ACK_ITEMS = [
    {
        "id": "lawful_use",
        "label": "I will not use the MARGIE pipeline, or its outputs, for any purpose "
                 "prohibited by applicable law or by the licenses of the underlying tools.",
    },
    {
        "id": "attribution",
        "label": "I will cite the MARGIE developers and authors, and the authors of the "
                 "underlying tools, whenever this work or its outputs are used or published.",
    },
    {
        "id": "third_party_licenses_obtained",
        "label": "For every tool that requires a separate license or permission (Phobius, "
                 "SignalP 4.x, SignalP 6.0, MEROPS; and, for commercial use, TMbed/ProtT5 "
                 "weights, TCDB, and the full KEGG database), I have already obtained the "
                 "necessary license(s) or permission(s) directly from the provider.",
    },
    {
        "id": "record_keeping_consent",
        "label": "I understand and agree that an exact copy of these accepted terms — with my "
                 "username, the UTC date/time of acceptance, and the IP address I accepted "
                 "from — will be recorded and sent to the MARGIE developers for legal "
                 "record-keeping.",
    },
]
_REQUIRED_ACK_IDS = {item["id"] for item in ACK_ITEMS}

# How the user intends to use MARGIE. Drives whether commercial_restricted tools
# are gated. Keep ids stable — stored in each acceptance record.
USAGE_ACADEMIC = "academic"
USAGE_COMMERCIAL = "commercial"
USAGE_TYPES = [
    {"id": USAGE_ACADEMIC, "label": "Academic or non-profit research"},
    {"id": USAGE_COMMERCIAL,
     "label": "Commercial or other (for-profit, government contract, or any non-academic use)"},
]
_VALID_USAGE = {USAGE_ACADEMIC, USAGE_COMMERCIAL}

# Tiers that require the user to hold their own license/permission for a tool.
_BLOCKED_TIER = "blocked"
_COMMERCIAL_TIER = "commercial_restricted"


def gated_tool_ids() -> dict[str, set[str]]:
    """Tool ids grouped by the two gated tiers ({'blocked': {...}, 'commercial_restricted': {...}})."""
    catalog = load_catalog()
    out: dict[str, set[str]] = {_BLOCKED_TIER: set(), _COMMERCIAL_TIER: set()}
    for t in catalog.get("tools", []):
        tier = t.get("tier")
        if tier in out:
            out[tier].add(t["id"])
    return out


def disabled_tool_ids(usage_type: str | None, licensed_ids: Iterable[str] | None) -> set[str]:
    """Tool ids the user is NOT entitled to run, given their usage type and the
    tools they have licensed themselves.

    - ``blocked`` tools (Phobius, SignalP 4/6, MEROPS): disabled unless the user
      has licensed that specific tool — this holds for everyone, academic or
      not, because "blocked" means you must obtain your own copy.
    - ``commercial_restricted`` tools (TMbed, TCDB, KEGG): free for academic /
      non-profit use; for commercial use, disabled unless the user has licensed
      that specific tool.
    """
    licensed = set(licensed_ids or [])
    gated = gated_tool_ids()
    disabled = {tid for tid in gated[_BLOCKED_TIER] if tid not in licensed}
    if (usage_type or "").lower() == USAGE_COMMERCIAL:
        disabled |= {tid for tid in gated[_COMMERCIAL_TIER] if tid not in licensed}
    return disabled


def get_entitlement(username: str) -> dict:
    """Return {usage_type, licensed_tools} from the user's CURRENT-terms
    acceptance (usage_type=None, licensed_tools=[] if they haven't accepted)."""
    current = load_terms()["version"]
    with get_db() as db:
        row = db.execute(
            "SELECT usage_type, licensed_tools FROM license_acceptances "
            "WHERE username = ? AND terms_version = ? ORDER BY id DESC LIMIT 1",
            (username, current),
        ).fetchone()
    if not row:
        return {"usage_type": None, "licensed_tools": []}
    try:
        licensed = json.loads(row["licensed_tools"]) if row["licensed_tools"] else []
    except (TypeError, ValueError):
        licensed = []
    return {"usage_type": row["usage_type"], "licensed_tools": licensed}


def save_depot_record(username: str, record: dict, terms_text: str,
                      timestamp: str | None = None) -> str | None:
    """Best-effort mirror of an acceptance record to the depot records dir.
    Returns the path, or None on any failure (never raises) — used by the CLI
    gate, which must not fail a run just because the shared record dir is
    unreachable."""
    ts = timestamp or datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    try:
        return _write_record_file(
            username=username, record=record, terms_text=terms_text, timestamp=ts,
        )
    except Exception:  # noqa: BLE001 — best-effort only
        return None


def _records_dir() -> Path:
    return Path(os.getenv("LICENSE_RECORDS_DIR", _DEFAULT_RECORDS_DIR))


def _operator() -> str:
    return os.getenv("MARGIE_OPERATOR", "lindems")


def load_terms() -> dict:
    """Return {version, text, sha256} for the current terms.md."""
    text = _TERMS_PATH.read_text(encoding="utf-8")
    first_line = text.splitlines()[0] if text else ""
    m = _VERSION_RE.search(first_line)
    version = m.group(1) if m else "unknown"
    sha256 = hashlib.sha256(text.encode("utf-8")).hexdigest()
    return {"version": version, "text": text, "sha256": sha256}


def load_catalog() -> dict:
    """Return the parsed licensing_catalog.json."""
    return json.loads(_CATALOG_PATH.read_text(encoding="utf-8"))


def build_terms_payload() -> dict:
    """The full payload the frontend needs to render the gate."""
    terms = load_terms()
    catalog = load_catalog()
    tools = catalog.get("tools", [])
    gated = [t for t in tools if t.get("tier") in ("blocked", "commercial_restricted")]
    return {
        "terms_version": terms["version"],
        "terms_sha256": terms["sha256"],
        "terms_markdown": terms["text"],
        "acknowledgments": ACK_ITEMS,
        "usage_types": USAGE_TYPES,
        "catalog_version": catalog.get("catalog_version"),
        "gated_tools": gated,
        "tools": tools,
        "tier_definitions": catalog.get("tier_definitions", {}),
    }


def has_accepted_current_terms(username: str) -> bool:
    """True if the user has an acceptance row matching the CURRENT terms version."""
    current = load_terms()["version"]
    with get_db() as db:
        row = db.execute(
            "SELECT 1 FROM license_acceptances WHERE username = ? AND terms_version = ? LIMIT 1",
            (username, current),
        ).fetchone()
    return row is not None


def revoke_current_acceptance(username: str) -> int:
    """Delete the current-terms acceptance row(s) for a user so the license
    gate re-prompts them on the analyze page.

    Only the app-side gate state is cleared; the immutable depot record of
    the original acceptance is intentionally left in place as legal history
    (an acceptance genuinely happened, even if the user later revokes going
    forward). Returns the number of rows removed (0 if the user hadn't
    accepted the current terms)."""
    current = load_terms()["version"]
    with get_db() as db:
        cur = db.execute(
            "DELETE FROM license_acceptances WHERE username = ? AND terms_version = ?",
            (username, current),
        )
        return cur.rowcount


def _write_record_file(
    *, username: str, record: dict, terms_text: str, timestamp: str
) -> str:
    """Write the exact terms copy + machine-readable record under the records dir.

    Layout: <records_dir>/<timestamp>/<operator>/<username>/{terms.txt,record.json}
    Returns the directory path as a string. Best-effort: raises on failure so the
    caller can decide whether to still record the DB row.
    """
    operator = _operator()
    # sanitize username for a path segment
    safe_user = re.sub(r"[^A-Za-z0-9._-]", "_", username) or "user"
    dest = _records_dir() / timestamp / operator / safe_user
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "terms.txt").write_text(terms_text, encoding="utf-8")
    (dest / "record.json").write_text(
        json.dumps(record, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return str(dest)


def record_acceptance(
    *,
    username: str,
    accepted_items: list[str],
    ip_address: str | None,
    user_agent: str | None,
    usage_type: str,
    licensed_tools: list[str] | None = None,
) -> dict:
    """Validate + persist an acceptance. Returns {terms_version, accepted_at, depot_record_path}.

    Raises ValueError if not all required acknowledgments were accepted, or if
    usage_type is not one of USAGE_TYPES.
    """
    if not _REQUIRED_ACK_IDS.issubset(set(accepted_items)):
        missing = _REQUIRED_ACK_IDS - set(accepted_items)
        raise ValueError(f"All acknowledgments are required; missing: {sorted(missing)}")
    if usage_type not in _VALID_USAGE:
        raise ValueError(f"usage_type must be one of {sorted(_VALID_USAGE)}; got {usage_type!r}")
    licensed_tools = sorted(set(licensed_tools or []))

    terms = load_terms()
    accepted_at = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    timestamp_dir = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")

    record = {
        "username": username,
        "operator": _operator(),
        "accepted_at_utc": accepted_at,
        "ip_address": ip_address,
        "user_agent": user_agent,
        "terms_version": terms["version"],
        "terms_sha256": terms["sha256"],
        "catalog_version": load_catalog().get("catalog_version"),
        "accepted_acknowledgments": accepted_items,
        "usage_type": usage_type,
        "licensed_tools": licensed_tools,
    }

    # Write the exact copy to the depot records dir. If that fails we still
    # record the DB row (with a null path) so the acceptance is not lost, but
    # we surface the error path so operators can tell a copy is missing.
    depot_path = None
    depot_error = None
    try:
        depot_path = _write_record_file(
            username=username, record=record, terms_text=terms["text"],
            timestamp=timestamp_dir,
        )
    except Exception as exc:  # noqa: BLE001 — record must not be lost on FS error
        depot_error = f"{type(exc).__name__}: {exc}"

    with get_db() as db:
        db.execute(
            """INSERT INTO license_acceptances
               (username, terms_version, terms_sha256, accepted_items,
                ip_address, user_agent, accepted_at, depot_record_path,
                usage_type, licensed_tools)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (
                username, terms["version"], terms["sha256"],
                json.dumps(accepted_items), ip_address, user_agent,
                accepted_at, depot_path,
                usage_type, json.dumps(licensed_tools),
            ),
        )

    return {
        "terms_version": terms["version"],
        "accepted_at": accepted_at,
        "depot_record_path": depot_path,
        "depot_record_error": depot_error,
    }
