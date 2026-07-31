"""Licensing / license-gate support for the MARGIE analyze page.

This package holds the source-of-truth data (terms.md, licensing_catalog.json)
and the helpers that load them, compute the current terms version/hash, record
a user's acceptance, and check whether a user has accepted the current terms.
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
)
