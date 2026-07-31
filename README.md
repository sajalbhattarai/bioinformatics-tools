<h1 align="center">MARGIE</h1>
<p align="center"><b>M</b>ostly <b>A</b>utomated <b>R</b>apid <b>G</b>enome <b>I</b>nference <b>E</b>nvironment</p>
<p align="center"><i>Annotates prokaryotic genomes, finds operons, and scores how sure it is about each gene.</i></p>

<p align="center"><img src="docs/img/ecoli-s10-ribosomal-operon.png" width="900" alt="E. coli operon with per-gene confidence"></p>

---

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

## Licensing (first run)

The first time you run an analysis, MARGIE shows the licensing terms and **every tool/database's license**, then asks you to accept and to tell it:

- how you'll use MARGIE — **academic/non-profit** or **commercial**, and
- which license-required tools you've obtained yourself (Phobius, SignalP 4/6, MEROPS; plus TMbed/TCDB/KEGG for commercial use).

Tools you're not licensed for are **automatically disabled**. Do this once in a real terminal — a background/`nohup` first run will stop and ask you to accept interactively. Your acceptance is saved to `~/.config/bioinformatics-tools/` and archived (with the exact licenses shown) under `~/.local/share/bioinformatics-tools/licensing-records/`.

On first run the backend also creates its secret keys automatically in a git-ignored `.env` — there's no manual key step.

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
