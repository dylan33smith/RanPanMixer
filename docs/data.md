# data.md — the map

Read before touching data, runs or paths.

**Rule.** Every number is traceable to a (`<ARTIFACT>`, `<GRAPH+COHORT>`,
`<TAU>`, `<UTILITY>`, n) tuple. **If it is not in this file, it is not quotable.**

**State values:** `OK` (present, current) · `PLANNED` (not built) · `BUILDING`
(a running job is writing it — the row must name its `logs/` sentinel; a
half-written file is NOT quotable) · `MISSING` (referenced but absent) · `STALE`
(present but superseded) · `DEPRECATED` (renamed, retained, never quote).

Verified against disk 2026-09-18. `tests/test_docs_contract.py` re-verifies.

> **STANDING NOTE (2026-09-18).** This project still has **no mechanism code and
> no released paths**. What now exists is the upstream PanMixer checkout, its
> input data, and a verified map of it. Rows describing OUR outputs remain
> `PLANNED`.

> **STANDING BLOCKER — BUILD MISMATCH (2026-09-18).** PanMixer's released
> preprocessing pipeline computes LD blocks and variant mappings against the
> **GRCh37** 1000G Phase 3 panel while the pangenome is **GRCh38**. Measured on
> chr21: only **1,025** of 340,824 pangenome records match a Phase 3 record
> exactly on (POS, REF, ALT), and 2,176 on (POS, REF). Nothing derived from that
> match set is usable. The repo ships an **unwired** GRCh38 path
> (`external/PanMixer/starting_data/scripts/src/get_blocks_grch38.sh`) which is what we must use.
> The PUBLISHED paper pins the GRCh37 Phase 3 URL in its Data availability
> section while stating a GRCh38 backbone, so this inconsistency is in the paper,
> not only the code. The official GRCh38 liftover of Phase 3 was withdrawn in 2021.
> See `docs/bugs.md` and the 2026-09-18 entries in `docs/memory.md`.

---

## 1. Storage layout

| root | contents | state |
|---|---|---|
| `docs` | the six-file documentation set | `OK` |
| `paper` | the proposal PDF and manuscript drafts | `OK` |
| `archive_docs` | read-only source material: the framework PDF and the PanMixer paper | `OK` |
| `tests` | the docs contract | `OK` |
| `scripts` | our operational scripts: the local `sbatch` shim and pipeline runners. NOT mechanism code. | `OK` |
| `primer` | the onboarding primer: `primer/src/*.md` (five files: `00-whatsnew.md` plus the four sections; the editable source of truth), `primer/build.py` (assembles them into HTML), `primer/primer.html` (the built page, published as an Artifact). Teaching material, NOT one of the six working docs and not a source of project facts. | `OK` |
| `logs` | tmux run logs and status sentinels (gitignored) | `OK` |
| `external/PanMixer` | symlink -> `/data/ds85/RanPanMixer/external/PanMixer`. **PanMixer writes all its data inside its own tree** (`BASE_PATH = dirname(constants.py)`, no env override), and `/home` is 95% full, so the checkout lives on `/data`. Gitignored. | `OK` |
| `external/genetic_maps` | symlink -> `/data/ds85/RanPanMixer/genetic_maps`. GRCh38 genetic maps, acquired 2026-09-27. Working file `external/genetic_maps/chr21.b38.gmap`. Gitignored. | `OK` |
| `src/ranpanmixer` | the mechanism: HMM, utility, sampler, attacks, eval | `PLANNED` |

### The build (consolidated into `main`, 2026-09-29)

The `pipeline` orphan branch and its separate `RanPanMixer-pipeline/` worktree were
folded into `main` on 2026-09-29 — one directory, one branch. The branch is retained
in git as history.

| path | contents | state |
|---|---|---|
| `src/` | the code. `src/cohort/` builds from the graph and cohort alone and **takes no target argument**; `src/target/` is strictly downstream. `src/sites.py` owns the site axis and its digest, `src/store.py` is the only module that writes an artifact | `OK` |
| `configs/` | the input manifest and the four sample lists | `OK` |
| `tests/tiny/` | a 24-record hand-written fixture, the only one committable | `OK` |
| `data/` | symlink -> `/data/ds85/RanPanMixer/data`, gitignored | `OK` |

**Data tree.** `temp/` holds intermediates and is deletable; `cohort/` holds the
durable artifacts and no target ever touches it; `dev/` is a complete parallel root
reached with `--root data/dev`, so the 2 Mb slice runs the identical code path.

```
data/
  raw/            the two staged inputs
  temp/chr21/     deconstruct.vcf.gz, cohort44.vcf.gz   (deletable)
  cohort/chr21/   haplotypes reason missing_runs positions genetic_pos in_chain
                  support allele_lengths sites.tsv haplotype_ids manifest.json
  dev/            the same tree, built from the 14-16 Mb slice
  targets/<id>/chr21/   target.chr21.vcf.gz, path.npy, representability.json
  attack/         1000g_30x_chr21_attackdb.vcf.gz  (evaluation only)
```

⚠ **chr21 is a TESTING scope, not the product.** The tool must work genome-wide.
Preprocessing is already per-chromosome and parallel across chromosomes, so
extending it is a driver loop — but see `docs/plan.md` for the one thing that is
NOT just a loop: `beta_v` must sum to 1 over whatever scope `tau` covers, so a
genome-wide `tau` means a genome-wide normaliser, not 22 per-chromosome ones.

## 3b. OUR pipeline's input set (model C, decided 2026-09-28)

Three inputs, and only three. What model C removed is as important as what it kept.

| # | input | source | size | state |
|---|---|---|---|---|
| 1 | **HPRC v1.0 PGGB graph**, chr21 GFA | `.../freeze1/pggb/chroms/chr21.hprc-v1.0-pggb.gfa.gz` (note: graph files are under `chroms/`, NOT the `vcfs/` prefix) | **434,772,428 B**, sha256 `85222d4f...` | `OK` — staged 2026-09-29 |
| 2 | **GRCh38 genetic map** | SHAPEIT4, corroborated byte-identical by Beagle | `chr21.b38.gmap`, 971,965 B, sha256 `95557a75...`, 44,618 rows | `OK` — staged 2026-09-29 |
| 3 | **The target** | a held-out 1000G 30x sample, until a real external genome exists | one sample column | `OK` (pool selected) |

Staged on the `pipeline` branch under `data/raw/`, with the full manifest at
`pipeline:configs/inputs.chr21.yaml`. PGGB freeze1 publishes **no `.gbz`, `.og` or
`.snarls`** — GFA is the only form.

### Target selection (decided 2026-09-29)

Until a real external genome is available, the target is a **held-out sample from the
1000 Genomes 30x GRCh38 panel**. That panel is phased and GRCh38-called, so a sample
column *is* a path once joined to our site list — which is exactly the format the
pipeline wants, and it avoids `vg giraffe` entirely (HPRC marks short-read mapping
**"untested"** for the PGGB graph). The project already used this for the 2026-09-27
representability measurement.

**Eligibility — all four must hold:**

| rule | count |
|---|---|
| in the 3,202-sample 30x panel | 3,202 |
| NOT one of the 44 HPRC cohort members | −39 (the 39 present in the panel) |
| NOT a first-degree relative of any HPRC member | −78 |
| excluded from the attack database at evaluation time | — |
| **eligible pool** | **3,085** — `pipeline:configs/eligible_targets.chr21.txt` |

⚠ **The relatedness rule is not a formality.** The 39 HPRC donors in the panel are the
**children** of trios, and their **78 parents sit in the unrelated-2504 set** — so a
target drawn at random from the 2504 has a **~3.1% chance of being a cohort member's
parent**, sharing about half its genome with a donor haplotype. That would inflate
apparent fidelity and confound the privacy measurement. Verified against the 1000G
pedigree (`1kGP.3202_samples.pedigree_info.txt`), 2026-09-29.

⚠ **A FIFTH RULE, added 2026-09-29: no first-degree relative anywhere in the 3,202.** The
eligibility table above only excludes relatives of *HPRC* members. But if a target's own
parent or child sits in the **attack database**, the attack can match the relative instead
of the target, which confounds the measurement in the optimistic direction. Measured:
**1,686 of the 3,085 eligible samples (55%) have a first-degree relative in the panel.**
The strict pool is therefore **1,399**, at `pipeline:configs/clean_pool.chr21.txt`.

**THE TEST TARGET SET (selected 2026-09-29): 10 samples, two per super-population.** All
are in the strict 1,399 pool. HG00096 and HG00097 are carried over from the 2026-09-27
representability run for continuity; the other eight were drawn with a fixed seed
(`random.Random(20260929)`) so the selection is reproducible.

| sample | super-pop | pop | | sample | super-pop | pop |
|---|---|---|---|---|---|---|
| HG00096 | EUR | GBR | | HG03757 | SAS | STU |
| HG00097 | EUR | GBR | | HG04211 | SAS | ITU |
| HG01323 | AMR | PUR | | NA19056 | EAS | JPT |
| HG01468 | AMR | CLM | | NA19350 | AFR | LWK |
| HG02389 | EAS | CDX | | NA20294 | AFR | ASW |

Population spread is deliberate: the HPRC cohort is multi-continental, so a target's
ancestry relative to the cohort should strongly affect how well a mosaic of cohort
haplotypes can reproduce it. An all-GBR target set would hide that entirely.

⚠ **A held-out panel sample UNDERSTATES the representability problem.** Measured on these
six: ~80.4% of a target's chr21 non-reference calls have an exact record in the graph;
**9.9% have no record at all and carry 20.5% of the target's `-log f` information**. That
9.9% is a LOWER bound — a 1000G sample's variants are by construction already catalogued,
so a genuinely novel genome fares worse. Always report `target_fidelity` against this
ceiling.

**Dropped from the input set, with the reason:**

| dropped | why |
|---|---|
| PanGenie callset (`PG.vcf.gz`, 5.06 GB) | fed only the anchor set and one allele-frequency branch, both deleted by model C. No evaluator reads it (measured 2026-09-28) |
| the external 1000G panel, for the MECHANISM | fed only PLINK LD block boundaries, which model C does not use. Survives as the attack database in the evaluation only |
| `allele_frequencies.npy` and `get_af.py` | no privacy score and no allele-frequency sampling branch to feed |

⚠ **`vg deconstruct` is pipeline step 1, and our VCF will not match PanMixer's.** Measured
on chrY: `vg` 1.68 recovers 99.55% of the published VCF's sites by POS but only **92.7% by
(POS, REF, ALT)**, concentrated at multi-allelic sites, because the published file was built
with `vg` 1.36. Accepted knowingly. We must also apply chrX and `chm13` removal ourselves —
the published VCF carries 45 samples, we need 44. **If the head-to-head is ever run, both
arms must run on the SAME VCF.**

⚠ **The graph is published only as GFA for PGGB freeze1** — no `.gbz`, no `.snarls`. The two
`.gbz` files on disk under `/data/ds85/RanPanMixer/runs/` are outputs of our own
`vg autoindex` read-mapping baseline, built from the VCF plus `chr21_ucsc.fa`, so they carry
no information the VCF does not.

## 4. Record schema — released path

**PLANNED.** Fields fixed now so releases are self-describing from the first run.
A release that cannot answer "at what tau, under what utility" is not quotable.

| field | type | notes |
|---|---|---|
| `release_id` | str | ⚠ **not a unique key on its own.** The unique key is (`target_id`, `tau`, `utility_id`, `seed`). Collisions must raise, never overwrite. |
| `target_id` | str | pseudonymous target identifier |
| `tau` | float | the privacy parameter. Absent = the record is unusable. |
| `eta_tau` | float | stored, not recomputed at read time, so a calibration change is detectable |
| `utility_id` | str | identifies phi_t, its weights, and the beta_t scheme |
| `graph_id` | str | graph + cohort + chain segmentation. ⚠ NOT a leave-one-out id: our target is external, so nothing is held out |
| `loo_id` | str/null | the held-out stand-in, **comparison arm only**. `null` for a genuine external-target release |
| `seed` | int | the private random draw |
| `state_path` | int[T] | sampled Z_1:T |
| `log_z_p` | float | must satisfy 0 <= log_z_p <= eta_tau |
| `released_at` | date | one release per genome; re-releases compose |

**Integration contract (COMPARISON ARM ONLY).** Our sampler MAY also emit PanMixer's
artifact — `new_haplotypes.npy`, shape `(n_sites, 2)` int16 with `-1` for missing,
row-aligned to the chr VCF — because every PanMixer evaluator consumes exactly that file.
⚠ It is NOT a requirement of the mechanism, and emitting it does NOT make their whole
stack meaningful on our output: `af_loss`, `ld_loss` and the Beagle attack will run and
return numbers that mean nothing for a release which never edits the cohort. Only the
gap-score attack transfers, and only after its "self" becomes the external target's true
path and the target is excluded from the attack database.
⚠ The row-order coupling is UNCHECKED upstream — the shape assertion in their gap-score
evaluators is vacuous (both sides equal `site_mask.sum()` by construction). Ship a
`sites.tsv` (chrom, pos, ref, alt) beside any emitted copy and verify it before scoring.

## 5. Datasets — LIVE

| name | n | construction | leakage control |
|---|---|---|---|
| — | — | none of ours yet | — |

## 6. Datasets — DEPRECATED (DO NOT USE)

| path | why | state |
|---|---|---|
| anything derived from `pangenome_mask.npy` on the Phase 3 panel | The GRCh37/GRCh38 match set (1,025 chr21 sites). It feeds the attack DB, the `af_loss` MAF strata and the entire `ld_loss` site set. | `DEPRECATED` |
| `external/PanMixer/starting_data/chr21/DEPRECATED_1000g_phase3_grch37.vcf.gz` | The GRCh37 panel the shipped pipeline wires in. Retained so the deviation is visible and reversible, never to be quoted. | `DEPRECATED` |

## 7. Artifact registry

**Status key:** LIVE (quotable) · SUPERSEDED · INVALID (never quote)

| artifact | date | contents | status |
|---|---|---|---|
| `paper/Private_Genome_Path_Release.pdf` | 2026-08-29 | The proposal: mechanism, Theorem 1, Corollary 1, Algorithm 1, runtime analysis. | LIVE |
| `archive_docs/The_Six_File_Lab_Record.pdf` | 2026-08-29 | The documentation framework this repository follows. | LIVE |
| `paper/s41467-026-77591-0_reference.pdf` | 2026-09-18 | **The PUBLISHED PanMixer paper** (Nature Communications, DOI 10.1038/s41467-026-77591-0, Article in Press, 12 pages, sha256 `2bbba736...`), supplied by Dylan. **The version of record — all reproduction targets come from here.** | LIVE |
| `archive_docs/Blindenbach2026_PanMixer.pdf` | 2026-09-18 | The bioRxiv PREPRINT, 24 pages. Superseded: it says 47 individuals (published: 44), `eps_private` 0.001 (published: 0.002), AF divergence 0.004/0.004/0.002 (published: 0.006/0.006/0.005), and it has no Data availability section and no membership-inference analysis. Retained because this project's first analysis was made against it. | SUPERSEDED |
| `logs/pm_download.sha256` | 2026-09-18 | sha256 of the three downloaded VCFs. | LIVE |
| `primer/primer.html` + `primer/src/` | 2026-09-20 | Teaching primer for onboarding: a what-changed preface plus the landscape, the problem, PanMixer dissected, and our mechanism. ~57k words, built from `primer/src/` by `primer/build.py`. v2 expanded the hidden-Markov-model subsection (1.2.2) with full derivations of the forward algorithm and backward sampling, worked arithmetic, and a brute-force cross-check. URL: https://claude.ai/artifact/HzUWmUmxAPwXUSq3oFy8sb | LIVE |

## 8. Reproducibility contract

A result is reproducible here only if all five are recorded: the graph+cohort id,
`tau`, the utility id, the number of draws, and the seed. Two results are comparable only
when the first three match. ⚠ A held-out individual is recorded **only for comparison-arm
runs** (`loo_id`); a genuine external-target release holds out nobody, and requiring one
would make a condition no real release can satisfy into a gate on every result.

**The private input is the target path.** The `PLANNED` target directory must
never be committed, copied into an artifact directory, or quoted in any document.
Nothing in it is covered by the release guarantee; only the sampled output path
is. Written as a rule now, before the directory exists, because the first run is
where it would be broken.
