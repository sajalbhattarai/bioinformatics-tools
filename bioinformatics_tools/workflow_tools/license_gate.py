"""Command-line licensing gate for the MARGIE pipeline.

A run must not proceed until the operator has accepted the current licensing
terms. There are two ways that acceptance can already be satisfied:

  1. The web app (dane-api) verifies acceptance in its database before it SSHes
     in and runs `dane_wf`. It passes the result down as environment variables
     (see ENV_* below), so this gate never re-prompts a web-initiated run — it
     just honours the entitlement the web side already recorded.

  2. A previous interactive CLI acceptance, saved to
     ~/.config/bioinformatics-tools/license-acceptance.json.

If neither applies, an interactive run prompts the user to read and accept the
terms (and record their usage type + which license-required tools they hold);
a non-interactive run with no prior acceptance is refused with a clear message
rather than hanging on input.

The entitlement (usage type + self-licensed tools) then drives which tools are
disabled for the run — see ``disabled_tool_ids`` in the licensing catalog.
"""
from __future__ import annotations

import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

from bioinformatics_tools.api.licensing import catalog

# Set by dane-api on the remote `dane_wf` command once it has confirmed the
# user's acceptance in its own database.
ENV_ACCEPTED = "MARGIE_LICENSE_ACCEPTED"   # terms version (or any truthy value)
ENV_USAGE = "MARGIE_USAGE_TYPE"            # "academic" | "commercial"
ENV_LICENSED = "MARGIE_LICENSED_TOOLS"     # comma-separated tool ids


class LicenseError(Exception):
    """Raised when licensing has not been accepted and cannot be obtained."""


def _config_path() -> Path:
    base = os.getenv("XDG_CONFIG_HOME") or (Path.home() / ".config")
    return Path(base) / "bioinformatics-tools" / "license-acceptance.json"


def _split_csv(value: str) -> list[str]:
    return [x.strip() for x in value.split(",") if x.strip()]


def _load_local() -> dict | None:
    path = _config_path()
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def ensure_cli_license() -> dict:
    """Return the caller's entitlement ``{usage_type, licensed_tools, source}``.

    Prompts for first-time acceptance when run interactively. Raises
    ``LicenseError`` if the terms have not been accepted and we cannot prompt.
    """
    # 1) Web/API path — acceptance already verified; honour the passed entitlement.
    if os.environ.get(ENV_ACCEPTED):
        return {
            "usage_type": os.environ.get(ENV_USAGE) or None,
            "licensed_tools": _split_csv(os.environ.get(ENV_LICENSED, "")),
            "source": "api",
        }

    current = catalog.load_terms()["version"]

    # 2) A saved CLI acceptance of the current terms.
    saved = _load_local()
    if saved and saved.get("terms_version") == current:
        return {
            "usage_type": saved.get("usage_type"),
            "licensed_tools": saved.get("licensed_tools", []),
            "source": "cli-saved",
        }

    # 3) First CLI run (or terms changed): prompt, or refuse if non-interactive.
    if not (sys.stdin.isatty() and sys.stdout.isatty()):
        raise LicenseError(
            "MARGIE licensing terms have not been accepted on this machine. "
            "Run the pipeline once in an interactive terminal to read and accept "
            "the terms, or use the web app (which records acceptance for you)."
        )
    return _interactive_accept(current)


# --------------------------------------------------------------------------- #
# Interactive acceptance
# --------------------------------------------------------------------------- #
def _ask_required(prompt: str) -> str:
    while True:
        value = input(prompt).strip()
        if value:
            return value
        print("  (required)")


def _ask_yes_no(prompt: str) -> bool:
    while True:
        value = input(f"{prompt} [y/N]: ").strip().lower()
        if value in ("y", "yes"):
            return True
        if value in ("", "n", "no"):
            return False
        print("  Please answer y or n.")


def _print_license_disclosure(cat: dict) -> None:
    """Print the license of every third-party tool / database, verbatim from the
    catalog, so the user can review and accept them (and credit the authors)."""
    tools = sorted(cat.get("tools", []), key=lambda t: (t.get("phase", 0), t.get("name", "")))
    print("\n" + "=" * 74)
    print("PER-TOOL / DATABASE LICENSE DETAILS")
    print("=" * 74)
    print("MARGIE runs these third-party tools and databases, each under its own")
    print("license. By continuing you agree to use each responsibly, within its")
    print("license, and to credit its authors.")
    for t in tools:
        allowed = "; ".join(t.get("allowed") or [])
        not_allowed = "; ".join(t.get("not_allowed") or [])
        print("\n" + "-" * 74)
        print(f"  Tool / database:       {t.get('name', t.get('id', ''))}")
        print(f"  License type:          {t.get('license', '')}")
        print(f"  Academic use:          {t.get('academic_use', '')}")
        print(f"  Research use:          {t.get('research_use', '')}")
        print(f"  Commercial use:        {t.get('commercial_use', '')}")
        print(f"  Permission required:   {t.get('user_action', '')}")
        if allowed:
            print(f"  You may:               {allowed}")
        if not_allowed:
            print(f"  You may not:           {not_allowed}")
        if t.get("citation"):
            print(f"  Cite:                  {t['citation']}")
        if t.get("obtain_url"):
            print(f"  License / source:      {t['obtain_url']}")
        if t.get("license_quote"):
            tag = "verbatim" if t.get("license_quote_kind") == "verbatim" else "summary"
            print(f"  License notice ({tag}):")
            for ln in t["license_quote"].splitlines():
                print(f"      {ln}" if ln else "")
    print("\n" + "=" * 74)


def _interactive_accept(current_version: str) -> dict:
    terms = catalog.load_terms()
    cat = catalog.load_catalog()
    gated = catalog.gated_tool_ids()
    by_id = {t["id"]: t for t in cat.get("tools", [])}

    print("\n" + "=" * 74)
    print("MARGIE PIPELINE — LICENSING TERMS")
    print("=" * 74 + "\n")
    print(terms["text"])
    _print_license_disclosure(cat)
    print("\n" + "-" * 74)
    print("Your acknowledgments:")
    for item in catalog.ACK_ITEMS:
        print(f"  - {item['label']}")
    print("-" * 74)
    if not _ask_yes_no("\nDo you accept ALL of the acknowledgments above?"):
        raise LicenseError("You must accept the licensing terms to run MARGIE.")

    # Who is accepting (for the record).
    name = _ask_required("\nYour name: ")
    email = _ask_required("Your email: ")

    # Usage type.
    print("\nHow will you use MARGIE?")
    for i, u in enumerate(catalog.USAGE_TYPES, start=1):
        print(f"  {i}. {u['label']}")
    usage_type = None
    while usage_type is None:
        choice = input(f"Choose 1-{len(catalog.USAGE_TYPES)}: ").strip()
        if choice.isdigit() and 1 <= int(choice) <= len(catalog.USAGE_TYPES):
            usage_type = catalog.USAGE_TYPES[int(choice) - 1]["id"]
        else:
            print("  Please enter a valid number.")

    # Which license-required tools has the user obtained themselves?
    to_ask = sorted(gated["blocked"])
    if usage_type == catalog.USAGE_COMMERCIAL:
        to_ask += sorted(gated["commercial_restricted"])
    licensed: list[str] = []
    if to_ask:
        print("\nSome tools need a license you obtain yourself. For each, say whether")
        print("you have already obtained it. Tools you have NOT licensed are disabled.")
        for tid in to_ask:
            tool = by_id.get(tid, {})
            label = tool.get("name", tid)
            url = tool.get("obtain_url")
            suffix = f"  (obtain: {url})" if url else ""
            if _ask_yes_no(f"  Do you have a license for {label}?{suffix}"):
                licensed.append(tid)

    # Build + persist the record.
    accepted_at = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    timestamp_dir = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    import getpass
    import socket
    os_user = getpass.getuser()
    record = {
        "name": name,
        "email": email,
        "os_user": os_user,
        "hostname": socket.gethostname(),
        "operator": os.getenv("MARGIE_OPERATOR", "lindems"),
        "accepted_at_utc": accepted_at,
        "terms_version": terms["version"],
        "terms_sha256": terms["sha256"],
        "catalog_version": cat.get("catalog_version"),
        "accepted_acknowledgments": [item["id"] for item in catalog.ACK_ITEMS],
        "usage_type": usage_type,
        "licensed_tools": sorted(set(licensed)),
        # Exact snapshot of the per-tool license details shown, recorded to depot.
        "license_catalog": cat,
        "source": "cli",
    }

    # Best-effort mirror to the shared depot record dir (legal record-keeping).
    depot_path = catalog.save_depot_record(
        username=re_safe(name or os_user), record=record,
        terms_text=terms["text"], timestamp=timestamp_dir,
    )
    record["depot_record_path"] = depot_path
    # Durable, timestamped copy under the user's own data dir (always kept).
    local_archive = catalog.save_local_record(
        username=re_safe(name or os_user), record=record,
        terms_text=terms["text"], timestamp=timestamp_dir,
    )
    record["local_archive_path"] = local_archive

    # Local record (this is what future runs read to skip re-prompting).
    path = _config_path()
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(record, indent=2, ensure_ascii=False), encoding="utf-8")
    except OSError as exc:
        raise LicenseError(f"Could not save your acceptance to {path}: {exc}") from exc

    print("\nThank you — acceptance recorded:")
    print(f"  local record:   {path}")
    if local_archive:
        print(f"  local archive:  {local_archive}")
    if depot_path:
        print(f"  shared record:  {depot_path}")
    else:
        print("  (shared/depot copy not written here — your local copies are kept.)")
    print()
    return {
        "usage_type": usage_type,
        "licensed_tools": record["licensed_tools"],
        "source": "cli",
    }


def re_safe(value: str) -> str:
    import re
    return re.sub(r"[^A-Za-z0-9._-]", "_", value) or "user"
