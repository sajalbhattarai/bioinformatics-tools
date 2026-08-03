<div align="center">

# MARGIE

**Mostly Automated Rapid Genome Inference Environment**

Annotates prokaryotic genomes, finds operons, and scores how confident it is about each gene.

<img src="docs/img/ecoli-s10-ribosomal-operon.png" width="820" alt="E. coli operon with per-gene confidence">

[![Hosted app](https://img.shields.io/badge/Hosted_app-Open-2ea44f?style=for-the-badge)](https://bsp.anvilcloud.rcac.purdue.edu/)
[![Front-end (GUI)](https://img.shields.io/badge/Front--end-biolab--fe-1f6feb?style=for-the-badge)](https://github.com/sajalbhattarai/biolab-fe)
[![Built with Snakemake](https://img.shields.io/badge/Built_with-Snakemake-039475?style=for-the-badge)](https://snakemake.readthedocs.io/)

<a href="#install"><b>Install</b></a> &nbsp;|&nbsp;
<a href="#configure"><b>Configure</b></a> &nbsp;|&nbsp;
<a href="#run"><b>Run</b></a> &nbsp;|&nbsp;
<a href="#licensing"><b>Licensing</b></a> &nbsp;|&nbsp;
<a href="#prefer-a-web-app"><b>Web app</b></a> &nbsp;|&nbsp;
<a href="#ai-usage-in-the-project"><b>AI usage in the project</b></a>

</div>

<details>
<summary><b>Tip: Prefer the web app?</b></summary>

Prefer clicking to typing? MARGIE also has a **web interface**. Use the hosted app at **[bsp.anvilcloud.rcac.purdue.edu](https://bsp.anvilcloud.rcac.purdue.edu/)** with nothing to install, or run your own from **[biolab-fe](https://github.com/sajalbhattarai/biolab-fe)**. This page covers the **command-line** backend.

</details>

## Repo scope

This repository is the **backend** repository for MARGIE.
It contains backend code and backend operations only (pipeline, CLI, API, licensing gate, and backend runtime/config behavior).
It does **not** contain frontend GUI code.

For frontend setup and browser usage, use **[biolab-fe](https://github.com/sajalbhattarai/biolab-fe)**.

## Install

```bash
git clone path-to-repository ~/bioinformatics-tools
cd ~/bioinformatics-tools
uv sync
source .venv/bin/activate
```

## Configure

Run it once to create your config file:

```bash
dane_wf                        # creates ~/.config/bioinformatics-tools/config.yaml
```

Open that file and set where results are saved:

```yaml
main_database: path-to-sqlite-database.db
```

`main_database` is the SQLite database file used to store run metadata and results.
It is different from the `output_dir` folder used for generated run files, and different from fingerprint-related databases.

Databases default to `/depot/lindems/...`. To use your own, set those paths in the config — empty folders fill up as you run.

## Run

```bash
dane_wf margie sb \
  input: path-to-input-folder-or-fasta \
  output_dir: path-to-output-folder \
  run_full_operon_map: true
```

Run it in the background so it survives logout:

```bash
nohup dane_wf margie sb input: path-to-input-folder-or-fasta output_dir: path-to-output-folder run_full_operon_map: true > margie.log 2>&1 &
```

Watch it: `tail -f margie.log`  —  Stop it: `pkill -f dane_wf`

## Licensing

On the **first run**, MARGIE shows the licensing terms and **every tool/database's license**, then asks you to accept and to tell it:

- how you'll use MARGIE — **academic/non-profit** or **commercial**, and
- which license-required tools you've obtained yourself (Phobius, SignalP 4/6, MEROPS; plus TMbed/TCDB/KEGG for commercial use).

Tools you're not licensed for are **automatically disabled**.

> [!IMPORTANT]
> Accept it once in a real terminal — a background/`nohup` first run will stop and ask you to accept interactively.

Your acceptance is saved to `~/.config/bioinformatics-tools/` and archived (with the exact licenses shown) under `~/.local/share/bioinformatics-tools/licensing-records/`. The backend's secret keys are also created automatically in a git-ignored `.env` on first run — there is no manual key step.

## Prefer a web app?

This backend can be used from the command line or from the separate frontend repository, **[biolab-fe](https://github.com/sajalbhattarai/biolab-fe)**.

- **Just use it** — the hosted app at **[bsp.anvilcloud.rcac.purdue.edu](https://bsp.anvilcloud.rcac.purdue.edu/)**. Nothing to install; sign in and submit.
- **Run your own** — clone `biolab-fe` and run its one-time `setup.sh`. That frontend launcher starts this backend on your HPC over SSH, opens a tunnel, and launches the app in your browser.

Both routes drive this backend's `dane-api`, and both use the same licensing acceptance described above.

## Notes

- `margie sb` is two words. `output_dir:` is the exact folder (no timestamp added).
- Resume a crashed run: add `margie_sb.resume: true`.
- Run only some tools: add `margie_sb.selected_tools: prodigal,rast,pfam`.
- Full toolkit (file tools, API, front-end): [README-detailed.md](README-detailed.md).

## AI usage in the project

Phase 9-12 scripts were designed and implemented by **Sajal Bhattarai**.
During script development, **Claude Sonnet 4.6** was used in interactive mode to improve robustness and debug issues.
The core ideas, architecture, and intended behavior were defined by Sajal Bhattarai.
These scripts were manually validated for intended behavior.

Visualization and LLM work, including the operon circular diagram page, HTML creation, and interactive chat mode, were refined with interactive-mode assistance from **Claude Opus 4.8**.
These components were also manually checked and validated for intended purpose.

## Disclaimer

This software is provided "as is", without warranty of any kind, express or implied, including but not limited to warranties of merchantability, fitness for a particular purpose, and noninfringement.

## Citation

APA 7th (software):

Bhattarai, S., Deemer, D., & Lindemann, S. (2026). *bioinformatics-tools* (Version 1.0.0) [Computer software]. https://github.com/sajalbhattarai/bioinformatics-tools

Please also cite the individual tools and databases you use in the MARGIE pipeline, in accordance with their licensing and citation requirements. The licensing gates during MARGIE runs provide the relevant licensing details, but you should still cross-check and confirm the requirements before publication.

For machine-readable citation metadata, see [CITATION.cff](CITATION.cff).
