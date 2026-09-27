#!/usr/bin/env bash
#
# Runs the API on 127.0.0.1:8000 with this computer as its cluster (after ./setup.sh)
# and creates a dev account that reaches it over SSH with the setup key.
#
#   ./start.sh          start (or reuse a running API) and print the sign-in
#   ./start.sh --stop   stop the API it started
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
STATE="$HOME/.margie-dev"
KEY="$HOME/.ssh/margie_dev_ed25519"
API="http://127.0.0.1:8000"
PIDFILE="$STATE/api.pid"
LOG="$STATE/api.log"
ACCOUNT="$STATE/account"          # username and password of the dev account
mkdir -p "$STATE"; chmod 700 "$STATE"

# Finds the API by its port, since $! from a backgrounded subshell is not the server's pid.
api_pid() {
    local p
    p="$(lsof -nP -iTCP:8000 -sTCP:LISTEN -t 2>/dev/null | head -1)"
    [[ -n "$p" ]] && ps -o command= -p "$p" | grep -q 'bioinformatics_tools.api.main' && echo "$p"
}

if [[ "${1:-}" == "--stop" ]]; then
    pid="$(api_pid || true)"
    if [[ -n "$pid" ]]; then
        kill "$pid"; rm -f "$PIDFILE"
        # Waits for the process to exit so an immediate restart does not reuse it.
        for _ in $(seq 1 40); do kill -0 "$pid" 2>/dev/null || break; sleep 0.25; done
        echo "API stopped (pid $pid)."
    else
        echo "No API from this checkout is listening on 8000."
    fi
    exit 0
fi

[[ -f "$KEY" ]] || { echo "Run ./setup.sh first (no $KEY)." >&2; exit 1; }
[[ -f "$REPO/.env" ]] || { echo "No $REPO/.env: copy .env.example and set BSP_SECRET_KEY and BSP_ENCRYPTION_KEY." >&2; exit 1; }

# 1. The API on localhost; BSP_SKIP_DANE_WF_SYNC keeps submissions from resetting
#    ~/bioinformatics-tools (this checkout) to origin.
if curl -sf "$API/health" >/dev/null 2>&1; then
    echo "API already running at $API"
else
    echo "Starting the API at $API (log: $LOG)"
    ( cd "$REPO" && unset BSP_LOCAL_MODE && BSP_SKIP_DANE_WF_SYNC=1 nohup .venv/bin/python -c \
        "from bioinformatics_tools.api.main import serve; serve(host='127.0.0.1', port=8000)" \
        > "$LOG" 2>&1 & )
    for _ in $(seq 1 40); do curl -sf "$API/health" >/dev/null 2>&1 && break; sleep 0.5; done
    curl -sf "$API/health" >/dev/null || { echo "The API did not come up; see $LOG" >&2; tail -20 "$LOG" >&2; exit 1; }
    api_pid > "$PIDFILE" || true
fi

# 2. A dev account whose cluster is this computer; registering SSHes in with the key.
if [[ -f "$ACCOUNT" ]]; then
    read -r USERNAME PASSWORD < "$ACCOUNT"
else
    USERNAME="dev-$USER"
    # Uses python rather than a tr|head pipeline, whose SIGPIPE would trip pipefail.
    PASSWORD="$(python3 -c 'import secrets; print(secrets.token_urlsafe(15))')"
    body="$(python3 -c 'import json,sys; print(json.dumps({"username":sys.argv[1],"password":sys.argv[2],"cluster_host":"localhost","cluster_username":sys.argv[3],"private_key":open(sys.argv[4]).read()}))' \
        "$USERNAME" "$PASSWORD" "$USER" "$KEY")"
    code="$(curl -s -o "$STATE/register.out" -w '%{http_code}' -X POST "$API/v1/auth/register" -H 'Content-Type: application/json' -d "$body")"
    if [[ "$code" != 201 ]]; then
        echo "Registering the dev account failed (HTTP $code):" >&2; cat "$STATE/register.out" >&2; echo >&2; exit 1
    fi
    printf '%s %s\n' "$USERNAME" "$PASSWORD" > "$ACCOUNT"; chmod 600 "$ACCOUNT"
    echo "Registered the dev account $USERNAME (its cluster: $USER@localhost)"
fi

# 3. A default config on the "cluster" if it has none yet.
TOKEN="$(curl -s -X POST "$API/v1/auth/login" -H 'Content-Type: application/json' \
    -d "{\"username\":\"$USERNAME\",\"password\":\"$PASSWORD\"}" | python3 -c 'import json,sys; print(json.load(sys.stdin).get("access_token",""))')"
[[ -n "$TOKEN" ]] || { echo "Signing in as $USERNAME failed." >&2; exit 1; }
# A missing config returns 200 with an empty body.
cfg="$(curl -s "$API/v1/ssh/config" -H "Authorization: Bearer $TOKEN")"
if [[ -z "$(printf '%s' "$cfg" | python3 -c 'import json,sys; d=json.load(sys.stdin); print(d.get("main_database") or (d.get("config") or {}).get("main_database") or "")' 2>/dev/null)" ]]; then
    curl -sf -X POST "$API/v1/ssh/config/create-default" -H "Authorization: Bearer $TOKEN" >/dev/null \
        && echo "Wrote a default config at ~/.config/bioinformatics-tools/config.yaml"
fi

cat <<MSG

Ready. In the front end, sign in (cluster mode) as:
  username  $USERNAME
  password  $PASSWORD
(kept in $ACCOUNT)

Jobs run on this computer: sbatch starts them as local processes, squeue and
sacct report them, and their output lands where the API tells them to write.
Stop the API with: ./start.sh --stop
MSG
