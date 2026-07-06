<!-- terms_version: 2026-07-01 -->
# MARGIE Pipeline — Licensing Terms & Acknowledgment

**Please accept the licensing terms before use.**

The MARGIE pipeline orchestrates many third-party bioinformatics tools and
databases, each governed by its own license. Some are permissive; several are
**not**. Before you run any analysis, please read and accept the terms below.

---

## 1. Tools you must license yourself

The following components are **academic-use-only and cannot be redistributed by
MARGIE**. You must obtain your own copy or permission **directly from the
provider** before using them:

- **Phobius** — register at <https://phobius.sbc.su.se/>
- **SignalP 4.x** — register at <https://services.healthtech.dtu.dk/services/SignalP-4.1/>
- **SignalP 6.0** — register at <https://services.healthtech.dtu.dk/services/SignalP-6.0/>
- **MEROPS** (pepunit database) — register at <https://www.ebi.ac.uk/merops/>

## 2. Tools restricted for commercial use

The following are free for academic / non-profit use but **restricted or
prohibited for commercial use**. If your use is commercial, you must obtain the
appropriate permission or license first:

- **TMbed / ProtT5-XL-U50 model weights** — CC-BY-NC-SA 4.0, **NonCommercial only**;
  commercial use prohibited (contact Rostlab).
- **TCDB** — commercial use requires permission (<https://www.tcdb.org/>).
- **Full KEGG database** — commercial use requires a KEGG license
  (<https://www.kegg.jp/kegg/legal.html>). The public KOfam subset used by MARGIE
  is free.

## 3. Citation

MARGIE and every underlying tool are the product of others' labour. You are
required to cite the MARGIE developers and authors, and the authors of each tool
whose output you use, in any resulting work or publication.

## 4. Source availability (open-source / AGPL)

MARGIE is a network service built on open-source software. In accordance with the
AGPL-3.0 license of eggNOG-mapper, its corresponding source is offered to you: the
upstream repository (https://github.com/eggnogdb/eggnog-mapper) and the build
recipe maintained in the MARGIE repository under `build-here/`. Corresponding
source and exact versions for all GPL/AGPL components are likewise available under
`build-here/`.

---

## Your acknowledgments

By accepting, you confirm that:

1. **Lawful use.** I will not use the MARGIE pipeline, or its outputs, for any
   purpose prohibited by applicable law or by the licenses of the underlying
   tools and databases.

2. **Attribution.** I will cite the MARGIE developers and authors, and the
   authors of the underlying tools, whenever this work or its outputs are used
   or published.

3. **Third-party licenses obtained.** For every tool that requires a separate
   license or permission (Phobius, SignalP 4.x, SignalP 6.0, and MEROPS; and,
   for commercial use, TMbed/ProtT5 weights, TCDB, and the full KEGG database), I
   have already obtained the necessary license(s) or permission(s) directly from
   the provider.

4. **Record-keeping & privacy notice.** I understand and agree that an exact copy
   of these accepted terms — together with my username, the date and time (UTC)
   of acceptance, and the IP address from which I accepted — will be recorded and
   sent to the MARGIE developers for legal record-keeping. This information is
   collected solely to document license acceptance and is retained for that
   purpose.

Accepting these terms is required to use the analyze page. Your acceptance is
recorded against your account; if these terms are updated, you will be asked to
review and accept the new version.

---

*This acknowledgment exists to keep the MARGIE developers, and the authors of the
tools MARGIE builds upon, protected while making better science openly available.
It is not a substitute for reading each tool's own license.*
