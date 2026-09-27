#!/usr/bin/env bash
#
# Make this computer its own "cluster", so the API can be developed without HPC.
#
# The API does its work by SSHing into a login node, running dane_wf from
# ~/bioinformatics-tools, and submitting to SLURM. On a laptop all three can be
# the laptop itself:
#
#   ssh        Remote Login (sshd) on, and a key kept only for this, allowed in
#              to your own account (~/.ssh/authorized_keys, tagged margie-dev)
#   dane_wf    ~/bioinformatics-tools linked to this checkout -- the same link
#              the API makes for itself on start (main.py) when it is safe to
#   SLURM      sbatch / squeue / sacct / scancel from ./bin, on the PATH that
#              non-interactive SSH commands see (~/.zshenv, tagged margie-dev)
#
# Everything it adds is tagged, and --undo takes it all out again.
#
#   ./setup.sh            set it up (idempotent)
#   ./setup.sh --dry-run  say what it would change, change nothing
#   ./setup.sh --undo     remove everything it added
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
BIN="$HERE/bin"
KEY="$HOME/.ssh/margie_dev_ed25519"
AUTH="$HOME/.ssh/authorized_keys"
ENVFILE="$HOME/.zshenv"
LINK="$HOME/bioinformatics-tools"
TAG="margie-dev"
DRY=no
MODE=setup
for a in "$@"; do
    case "$a" in
        --dry-run) DRY=yes ;;
        --undo)    MODE=undo ;;
        -h|--help) sed -n '2,24p' "$0"; exit 0 ;;
        *) echo "unknown option: $a (see --help)" >&2; exit 2 ;;
    esac
done

say()  { printf '  %s\n' "$*"; }
# In a dry run, prints the action and returns 1 so the following && message is skipped.
act()  { if [[ $DRY == yes ]]; then say "would: $*"; return 1; else "$@"; fi; }

# ---- undo ----
if [[ $MODE == undo ]]; then
    echo "Removing the local dev cluster:"
    if [[ -f "$AUTH" ]] && grep -q " $TAG\$" "$AUTH"; then
        act sed -i '' "/ $TAG\$/d" "$AUTH" && say "key removed from $AUTH"
    fi
    if [[ -f "$ENVFILE" ]] && grep -q "# >>> $TAG" "$ENVFILE"; then
        act sed -i '' "/# >>> $TAG/,/# <<< $TAG/d" "$ENVFILE" && say "PATH line removed from $ENVFILE"
    fi
    if [[ -L "$LINK" && "$(readlink "$LINK")" == "$REPO" ]]; then
        act rm "$LINK" && say "$LINK removed (it pointed here)"
    fi
    for f in "$KEY" "$KEY.pub"; do [[ -f "$f" ]] && { act rm "$f" && say "$f removed"; }; done
    say "done. Remote Login is left as it was -- turn it off in System Settings if you only had it on for this."
    exit 0
fi

# ---- setup ----
echo "Setting up this computer as its own cluster:"

# 1. sshd: only checked, since turning it on needs an administrator.
if ! nc -z -G 2 127.0.0.1 22 2>/dev/null; then
    echo "  Remote Login is off. Turn it on in System Settings > General > Sharing > Remote Login," >&2
    echo "  allowing only your own account, then run this again." >&2
    exit 1
fi
say "Remote Login: on"

# 2. A dedicated key, allowed into the user's own account.
if [[ ! -f "$KEY" ]]; then
    act ssh-keygen -q -t ed25519 -N '' -C "$TAG" -f "$KEY" && say "made $KEY"
else
    say "key: $KEY (already there)"
fi
if [[ $DRY == no || -f "$KEY.pub" ]]; then
    PUB="$(cat "$KEY.pub" 2>/dev/null || true)"
    if [[ -n "$PUB" ]] && ! grep -qF "$PUB" "$AUTH" 2>/dev/null; then
        act mkdir -p "$HOME/.ssh"
        if [[ $DRY == yes ]]; then say "would: add the key to $AUTH"; else
            printf '%s\n' "${PUB% *} $TAG" >> "$AUTH"; chmod 600 "$AUTH"; say "key allowed in $AUTH"; fi
    else
        say "key already allowed in $AUTH"
    fi
fi

# 3. dane_wf: links ~/bioinformatics-tools to this checkout, never replacing a real folder.
if [[ -L "$LINK" && "$(readlink "$LINK")" == "$REPO" ]]; then
    say "$LINK -> this checkout (already)"
elif [[ -e "$LINK" && ! -L "$LINK" ]]; then
    echo "  $LINK is a real folder, not a link to this checkout. Leaving it alone;" >&2
    echo "  move it aside if it should track $REPO instead." >&2
    exit 1
else
    act ln -sfn "$REPO" "$LINK" && say "$LINK -> $REPO"
fi
[[ -x "$REPO/.venv/bin/dane_wf" ]] || { echo "  no .venv/bin/dane_wf in $REPO: run 'uv sync' there first" >&2; exit 1; }

# 4. SLURM: puts the stand-ins on the PATH that non-interactive SSH commands see.
if [[ -f "$ENVFILE" ]] && grep -q "# >>> $TAG" "$ENVFILE"; then
    say "SLURM stand-ins already on PATH ($ENVFILE)"
else
    if [[ $DRY == yes ]]; then say "would: put $BIN on PATH in $ENVFILE"; else
        printf '\n# >>> %s: the local dev cluster'"'"'s SLURM stand-ins (margie-backend/scripts/dev-local)\nexport PATH="%s:$PATH"\n# <<< %s\n' "$TAG" "$BIN" "$TAG" >> "$ENVFILE"
        say "SLURM stand-ins on PATH ($ENVFILE)"; fi
fi

# 5. Checks the setup over SSH with that key, as the API will use it.
if [[ $DRY == no ]]; then
    out="$(ssh -i "$KEY" -o BatchMode=yes -o StrictHostKeyChecking=accept-new -o ConnectTimeout=5 "$USER@localhost" \
        'command -v sbatch && test -x "$HOME/bioinformatics-tools/.venv/bin/dane_wf" && echo ready' 2>&1 || true)"
    if [[ "$out" == *ready ]]; then say "checked over SSH: sbatch and dane_wf both found"; else
        echo "  the SSH check did not pass:" >&2; printf '%s\n' "$out" | sed 's/^/    /' >&2; exit 1; fi
fi
echo "Done. Start the API with ./start.sh"
