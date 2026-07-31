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
<a href="#prefer-a-web-app"><b>Web app</b></a>

</div>

> [!TIP]
> Prefer clicking to typing? MARGIE also has a **web interface**. Use the hosted app at **[bsp.anvilcloud.rcac.purdue.edu](https://bsp.anvilcloud.rcac.purdue.edu/)** — nothing to install — or run your own from **[biolab-fe](https://github.com/sajalbhattarai/biolab-fe)**. This page covers the **command-line** backend.

## Install

```bash
git clone <github-url> ~/bioinformatics-tools
cd ~/bioinformatics-tools
uv sync
source .venv/bin/activate      # now the `dane_wf` command works
```

## Configure

Run it once to create your config file:

```bash
dane_wf                        # creates ~/.config/bioinformatics-tools/config.yaml
```

Open that file and set where results are saved:

```yaml
main_database: <path-to-your-results>.db
```

Databases default to `/depot/lindems/...`. To use your own, set those paths in the config — empty folders fill up as you run.

## Run

```bash
dane_wf margie sb \
    input: <folder of FASTA files, or one FASTA> \
    output_dir: <folder for results> \
    run_full_operon_map: true
```

Run it in the background so it survives logout:

```bash
nohup dane_wf margie sb input: <...> output_dir: <...> run_full_operon_map: true > margie.log 2>&1 &
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

The command line is one way to run MARGIE; the same pipeline is also driven by a browser front-end, **[biolab-fe](https://github.com/sajalbhattarai/biolab-fe)**.

- **Just use it** — the hosted app at **[bsp.anvilcloud.rcac.purdue.edu](https://bsp.anvilcloud.rcac.purdue.edu/)**. Nothing to install; sign in and submit.
- **Run your own** — clone `biolab-fe` and run its one-time `setup.sh`. It starts this backend on your HPC over SSH, opens a tunnel, and launches the app in your browser.

Both routes drive this backend's `dane-api`, and both use the same licensing acceptance described above.

## Notes

- `margie sb` is two words. `output_dir:` is the exact folder (no timestamp added).
- Resume a crashed run: add `margie_sb.resume: true`.
- Run only some tools: add `margie_sb.selected_tools: prodigal,rast,pfam`.
- Full toolkit (file tools, API, front-end): [README-detailed.md](README-detailed.md).

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
