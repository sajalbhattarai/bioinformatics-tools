"""License-gate support for the MARGIE analyze page.

Holds the terms and tool catalog (terms.md, licensing_catalog.json) and the
helpers that version them, record acceptance and check it.
"""
from .catalog import (  # noqa: F401
    ACK_ITEMS,
    USAGE_TYPES,
    build_terms_payload,
    disabled_tool_ids,
    gated_tool_ids,
    get_entitlement,
    has_accepted_current_terms,
    load_catalog,
    load_terms,
    record_acceptance,
    revoke_current_acceptance,
    save_depot_record,
    save_local_record,
)
