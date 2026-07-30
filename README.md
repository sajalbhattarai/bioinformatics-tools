<h1 align="center">MARGIE</h1>
<p align="center"><b>M</b>ostly <b>A</b>utomated <b>R</b>apid <b>G</b>enome <b>I</b>nference <b>E</b>nvironment</p>
<p align="center"><i>A pipeline that annotates prokaryotic genomes, finds operons, and scores how sure it is about each gene.</i></p>

<p align="center"><img src="docs/img/ecoli-s10-ribosomal-operon.png" width="900" alt="E. coli S10 operon with per-gene confidence"></p>
<p align="center"><sub><i>Example: an</i> E. coli <i>operon, with each gene's confidence shown underneath.</i></sub></p>

---

## What it does

MARGIE runs in four steps:

**1. Annotate genes → 2. Find operons → 3. Score confidence (C1–C4) → 4. Draw operon figures**

Steps 1–2 use existing tools; MARGIE joins them together and adds the confidence scoring. It works best on many genomes at once, but one genome is fine too. You decide what the results mean — MARGIE does not decide for you.

---

## Quickstart

**1. Install (once)**

```bash
git clone <github-url> ~/bioinformatics-tools
cd ~/bioinformatics-tools
uv sync
```

This creates the tool at `~/bioinformatics-tools/.venv/bin/dane_wf`.

**2. Set up your config (once)**

The first time you run `dane_wf`, it creates a config file at `~/.config/bioinformatics-tools/config.yaml`. That file controls everything — where results are saved, your cluster settings, and all database/container paths. Open it and set at least:

```yaml
main_database: <path-to-your-results>.db   # caches results so re-runs are fast
# db_root  defaults to /depot/lindems/data/margie/db  if not set
# sif_path defaults to /depot/lindems/data/margie/sif if not set
```

Re-running the same genome (same FASTA) restores the saved result instead of redoing the work.

Every database has a `/depot/lindems/...` default — point any of them at your own folder instead, and it builds itself up as you run genomes (no lab access needed).

**3. Run**

```bash
~/bioinformatics-tools/.venv/bin/dane_wf margie sb \
    input: <folder of genome FASTA files, or one FASTA> \
    output_dir: <folder for this run's results> \
    run_full_operon_map: true
```

- Type `margie sb` as two words (config keys still use `margie_sb`).
- Each option is `key: value` with a space after the colon.
- `output_dir:` is the exact folder — MARGIE does not add a timestamp.
- `run_full_operon_map: true` draws the full operon figure atlas (off by default).

**Run it unattended.** `dane_wf` must stay running the whole time, so don't tie it to your SSH session. Submit it as a small SLURM job (1 CPU, long time limit — the heavy work runs in its own jobs), or use `nohup`:

```bash
nohup ~/bioinformatics-tools/.venv/bin/dane_wf margie sb \
    input: <...> output_dir: <...> run_full_operon_map: true > margie.log 2>&1 &
```

Watch it: `tail -f margie.log`.  Stop it: `pkill -f dane_wf`.

**Extras**

| To do this | Add |
|---|---|
| Resume a crashed run | `margie_sb.resume: true` (point `output_dir:` at the old folder) |
| Run only some tools | `margie_sb.selected_tools: prodigal,rast,pfam` |

For the full toolkit (file tools, API server, front-end), see [README-detailed.md](README-detailed.md).

---

## Citation

```yaml
cff-version: 1.2.0
message: "If you use this software, please cite it as below."
authors:
  - family-names: "Bhattarai"
    given-names: "Sajal"
    orcid: "https://orcid.org/0000-0002-3143-5483"
  - family-names: "Deemer"
    given-names: "Dane"
    orcid: "https://orcid.org/0000-0002-4485-0280"
  - family-names: "Lindemann"
    given-names: "Stephen"
    orcid: "https://orcid.org/0000-0002-3788-5389"
title: "bioinformatics-tools"
version: 1.0.0
date-released: 2026-07-27
url: "https://github.com/Diet-Microbiome-Interactions-Lab"
```
