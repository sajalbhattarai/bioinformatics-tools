# Running the API with this computer as its cluster

For developing the API and the front end's cluster mode without an HPC
account. Nothing here is used by a real deployment.

The API does its work by SSHing into a login node, running `dane_wf` from
`~/bioinformatics-tools`, and submitting jobs to SLURM. All three can be this
computer:

| The API expects | Here it is |
|---|---|
| a login node it can SSH into | your own account, over Remote Login, with a key kept only for this |
| `~/bioinformatics-tools/.venv/bin/dane_wf` | a link to this checkout (the API makes the same link itself on start) |
| `sbatch`, `squeue`, `sacct`, `scancel` | the stand-ins in `bin/`: jobs are local processes, their records in `~/.margie-dev/slurm/` |
| `setsid` (Linux only) | a stand-in in `bin/`: the API starts login-node runs as `nohup setsid bash -c ...`, which dies at once on macOS without it |

This runs the API's **real** code path — SSH, the job runner, SLURM
submission and polling, the job history, the result cache — which is what
makes it useful for development. (`BSP_LOCAL_MODE`, by contrast, only changes
registration; every job endpoint still needs an SSH key it then does not
have.)

## Use

```sh
./setup.sh --dry-run   # what it would change
./setup.sh             # key, ~/bioinformatics-tools link, SLURM stand-ins on PATH
./start.sh             # API on 127.0.0.1:8000 and a dev account; prints how to sign in
./start.sh --stop
./setup.sh --undo      # take all of it out again
```

Remote Login must be on (System Settings › General › Sharing › Remote Login,
allowed for your account only). `setup.sh` checks, but does not change it.

## What it touches

Everything it adds is tagged `margie-dev` and removed by `--undo`:

- `~/.ssh/margie_dev_ed25519` — the key, used only for this
- one line in `~/.ssh/authorized_keys` — that key, allowed into your own account
- a block in `~/.zshenv` — puts `bin/` on the PATH non-interactive SSH commands see
- `~/bioinformatics-tools` — a link to this checkout (never replaces a real folder)
- `~/.margie-dev/` — the API's log and pid, the dev account, the stand-in job records

## What the stand-ins do

They implement exactly what `bioinformatics_tools/utilities/ssh_slurm.py`
calls — `sbatch --parsable`, the `squeue` and `sacct` formats it reads, and
`scancel` — and refuse any other flag rather than guess. A submitted script
runs detached in its submit folder with `SLURM_JOB_ID` set, writes to its
`--output`, and records `COMPLETED`, `FAILED` (non-zero exit) or `CANCELLED`.
There is no queueing, no partitions and no resource limits: everything runs
at once, on this machine.

## What it does and does not cover

It runs everything **around** the science for real: sign-in, SSH, the job
runner, submission and polling, logs, the job history, config, file browsing.
A `quick_example` self-test goes from submit to *completed* in about ten
seconds.

It does not run the annotation itself. The `margie` workflow's tools are
Apptainer images, and Apptainer does not run on macOS; a real run belongs on
the cluster (or on the local pipeline, which uses Apple's container).

A `margie_sb` submission does get through the API and into SLURM, and fails
on the Mac as it should (`Read-only file system: '/depot'`), which exercises
the job list, logs and failure reasons. For the API to accept one, its
pre-flight checks need, on this computer:

- a job-history database at `main_database` (an empty SQLite file will do:
  the API creates its table),
- the shared reference paths under `margie_sb` (operon reference,
  fingerprint database, genome pool, ...) pointing at files and folders that
  exist -- dummies under `~/.margie-dev/shared/` are enough,
- a SLURM account (any name: the stand-ins ignore it) and accepted licence
  terms.

## Your checkout is the cluster's checkout

`~/bioinformatics-tools` is a link to this repo, and every `margie_sb`
submission runs `sync_remote_dane_wf`: a `git fetch`, then -- when the tree
is clean and not at `origin/for-website-deployment` -- a checkout of that
branch and `git reset --hard` to it. On a real cluster that keeps users
current; here it would throw away the branch you are working on. `start.sh`
starts the API with `BSP_SKIP_DANE_WF_SYNC=1`, which switches the sync off
(the API log says "version-sync check disabled"). Start the API some other way
and set it yourself.

## Found while building it

- Completion of a login-node run was never detected over SFTP: the run's
  files are named `$HOME/...`, which a shell expands and SFTP does not. Fixed in
  `utilities/ssh_slurm.py`. SFTP expands `$HOME` on no server, so the check
  cannot have worked on the cluster either; there something else must have
  been settling finished jobs (not verified).
- `tests/test_api.py::TestJobStatusEndpoint::test_job_status_404` fails
  (500, not 404) with or without that fix.
- Something under `tests/` other than `test_api.py` hangs; `test_api.py` alone
  runs in under a second.
- Self-test runs complete but do not appear in `GET /v1/ssh/jobs`; that list
  reads the history database, which the self-tests may not write to.
- 13 tests in `tests/test_api.py` fail on the committed code as well (mocks
  that no longer match `list_jobs_and_count` and the resume/restart paths,
  and three auth tests); none are from the changes above.
