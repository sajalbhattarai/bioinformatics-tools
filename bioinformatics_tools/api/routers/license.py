"""License-gate endpoints.

The analyze page is gated: a user must accept the current licensing terms before
running any workflow. Acceptance is recorded once per account per terms version
(re-prompted when the terms change), stored in the app database, and an exact
copy is written to the depot records directory with the user's IP (disclosed in
the terms).

Routes (all under /v1/license, auth required):
  GET  /terms   → terms text + version/hash + tool catalog + acknowledgments
  GET  /status  → whether the current user has accepted the current terms
  POST /accept  → record acceptance (validates version/hash, captures IP), returns status
"""
import logging

from fastapi import APIRouter, Depends, HTTPException, Request, status
from pydantic import BaseModel

from bioinformatics_tools.api.auth import get_current_user
from bioinformatics_tools.api import licensing

LOGGER = logging.getLogger(__name__)

router = APIRouter(prefix="/v1/license", tags=["license"])


class LicenseAccept(BaseModel):
    accepted_items: list[str]          # acknowledgment ids the user checked
    terms_version: str                 # version the user was shown (must match current)
    terms_sha256: str | None = None    # optional integrity check of the shown terms
    usage_type: str                    # "academic" | "commercial"
    licensed_tools: list[str] = []     # tool ids the user holds their own license for


def _client_ip(request: Request) -> str | None:
    """Best-effort client IP. Honors X-Forwarded-For (first hop) behind a proxy."""
    fwd = request.headers.get("x-forwarded-for")
    if fwd:
        return fwd.split(",")[0].strip()
    real = request.headers.get("x-real-ip")
    if real:
        return real.strip()
    return request.client.host if request.client else None


@router.get("/terms")
def get_terms(current_user: dict = Depends(get_current_user)):
    """Return everything the frontend needs to render the license gate."""
    return licensing.build_terms_payload()


@router.get("/status")
def get_status(current_user: dict = Depends(get_current_user)):
    """Whether the current user has accepted the CURRENT terms version, plus the
    entitlement that drives which tools the analyze page must disable."""
    terms = licensing.load_terms()
    accepted = licensing.has_accepted_current_terms(current_user["username"])
    entitlement = licensing.get_entitlement(current_user["username"])
    disabled = sorted(
        licensing.disabled_tool_ids(
            entitlement.get("usage_type"), entitlement.get("licensed_tools")
        )
    )
    return {
        "accepted": accepted,
        "current_terms_version": terms["version"],
        "usage_type": entitlement.get("usage_type"),
        "licensed_tools": entitlement.get("licensed_tools", []),
        "disabled_tools": disabled,
    }


@router.post("/revoke")
def revoke_terms(current_user: dict = Depends(get_current_user)):
    """Revoke the current user's acceptance of the current terms, so the
    license gate re-prompts them next time they open the analyze page.

    Self-service: a user can revoke their own acceptance at any time (e.g.
    if their usage/eligibility status changes). The immutable depot record
    of the original acceptance is preserved for legal history -- this only
    clears the app-side gate state."""
    removed = licensing.revoke_current_acceptance(current_user["username"])
    terms = licensing.load_terms()
    LOGGER.info(
        "License revoked: user=%s version=%s rows_removed=%d",
        current_user["username"], terms["version"], removed,
    )
    return {"revoked": True, "rows_removed": removed, "current_terms_version": terms["version"]}


@router.post("/accept")
def accept_terms(
    body: LicenseAccept,
    request: Request,
    current_user: dict = Depends(get_current_user),
):
    """Record the current user's acceptance of the current terms."""
    terms = licensing.load_terms()

    # The user must be accepting the terms version we are actually serving.
    if body.terms_version != terms["version"]:
        raise HTTPException(
            status_code=status.HTTP_409_CONFLICT,
            detail=(
                "The licensing terms have changed since this page loaded. "
                "Please reload and review the current terms before accepting."
            ),
        )
    if body.terms_sha256 and body.terms_sha256 != terms["sha256"]:
        raise HTTPException(
            status_code=status.HTTP_409_CONFLICT,
            detail="Terms integrity check failed. Please reload and try again.",
        )

    ip_address = _client_ip(request)
    user_agent = request.headers.get("user-agent")

    try:
        result = licensing.record_acceptance(
            username=current_user["username"],
            accepted_items=body.accepted_items,
            ip_address=ip_address,
            user_agent=user_agent,
            usage_type=body.usage_type,
            licensed_tools=body.licensed_tools,
        )
    except ValueError as exc:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(exc))

    LOGGER.info(
        "License accepted: user=%s version=%s ip=%s depot=%s",
        current_user["username"], result["terms_version"], ip_address,
        result.get("depot_record_path"),
    )
    if result.get("depot_record_error"):
        # Acceptance is still valid (DB row written); flag that the depot copy failed.
        LOGGER.error(
            "License depot record write FAILED for user=%s: %s",
            current_user["username"], result["depot_record_error"],
        )

    return {"accepted": True, **result}
