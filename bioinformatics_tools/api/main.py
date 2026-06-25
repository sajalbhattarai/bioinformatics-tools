"""
FastAPI application entry point for Bioinformatics Tools API

Usage:
    Development: uvicorn bioinformatics_tools.api.main:app --reload
    Production: dane-api (after installing with pip install .[api])
"""
import logging
from pathlib import Path

from dotenv import load_dotenv

# Load project-root .env first, then allow shell environment to override values.
_PROJECT_ROOT = Path(__file__).resolve().parents[2]
load_dotenv(dotenv_path=_PROJECT_ROOT / '.env', override=False)
load_dotenv(override=False)  # Must run before local imports that read env vars at import time

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
import uvicorn

from bioinformatics_tools.api.database import init_db
from bioinformatics_tools.api.routers import auth, dane, fasta, ssh, workflows

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s %(name)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
LOGGER = logging.getLogger(__name__)

# Create FastAPI app
app = FastAPI(
    title='Bioinformatics Tools API',
    version="0.0.1",
    description="API for bioinformatics file processing and analysis",
    docs_url="/docs",
    redoc_url="/redoc"
)

# Add CORS middleware (adjust origins for production)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In production, specify your frontend domain
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.on_event('startup')
def startup_event():
    init_db()


# Include routers
app.include_router(auth.router)
app.include_router(fasta.router)
app.include_router(dane.router)
app.include_router(ssh.router)
app.include_router(workflows.router)

# Root endpoint
@app.get("/")
async def root():
    """Root endpoint with API information"""
    return {
        "status": "success",
        "message": "Bioinformatics Tools API",
        "version": "0.0.1",
        "docs": "/tbd",
        "endpoints": {
            "fasta": "/v1/fasta",
            "ssh_upload": "/v1/ssh"
        }
    }

@app.get("/health")
async def health():
    """Health check endpoint"""
    return {"status": "success", "message": "API is healthy"}


def _ensure_remote_deployment_symlink() -> None:
    """The /v1/ssh/run_workflow endpoint SSHes into the user's own cluster
    account and runs `uvx --from ~/bioinformatics-tools/ dane_wf ...` -- a
    deliberately separate code path from whatever process dane-api itself
    is running from, since the API and the SSH target may be different
    machines in a multi-user deployment. When they happen to be the same
    machine (a single-developer setup, or local testing), ~/bioinformatics-
    tools needs to actually point at this checkout, or the SSH-invoked
    workflow silently runs against a stale, disconnected copy instead of
    whatever's actually being worked on.

    Self-heals the common, safe case (missing, or a symlink pointing
    somewhere else) by linking it to this same checkout, every time the API
    starts. Never touches an existing REAL directory there -- only logs a
    warning, so a deliberate, separate deployment is never silently
    destroyed.
    """
    target = Path.home() / "bioinformatics-tools"
    try:
        if target.is_symlink():
            if target.resolve() == _PROJECT_ROOT.resolve():
                return
            LOGGER.info("Relinking stale %s -> %s (was -> %s)", target, _PROJECT_ROOT, target.resolve())
            target.unlink()
            target.symlink_to(_PROJECT_ROOT)
        elif target.exists():
            LOGGER.warning(
                "%s exists as a real directory, not a symlink to %s -- leaving it alone. "
                "If it should track this checkout automatically instead, replace it with: "
                "rm -rf %s && ln -s %s %s",
                target, _PROJECT_ROOT, target, _PROJECT_ROOT, target,
            )
        else:
            LOGGER.info("Creating %s -> %s", target, _PROJECT_ROOT)
            target.symlink_to(_PROJECT_ROOT)
    except OSError as exc:
        LOGGER.warning("Could not verify/create %s -> %s: %s", target, _PROJECT_ROOT, exc)


def serve(host: str = "0.0.0.0", port: int = 8000, reload: bool = False):
    """
    Entry point for running the API server
    """
    _ensure_remote_deployment_symlink()
    LOGGER.info(f"Starting Bioinformatics Tools API server on {host}:{port}")
    uvicorn.run(
        "bioinformatics_tools.api.main:app",
        host=host,
        port=port,
        reload=reload,
        log_level="info"
    )


if __name__ == "__main__":
    # For development: python -m bioinformatics_tools.api.main
    serve(reload=True)
