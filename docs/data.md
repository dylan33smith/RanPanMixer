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
| `primer` | the onboarding primer: `primer/src/*.md` (the four sections, the editable source of truth), `primer/build.py` (assembles them into HTML), `primer/primer.html` (the built page, published as an Artifact). Teaching material, NOT one of the six working docs and not a source of project facts. | `OK` |
| `logs` | tmux run logs and status sentinels (gitignored) | `OK` |
| `external/PanMixer` | symlink -> `/data/ds85/RanPanMixer/external/PanMixer`. **PanMixer writes all its data inside its own tree** (`BASE_PATH = dirname(constants.py)`, no env override), and `/home` is 95% full, so the checkout lives on `/data`. Gitignored. | `OK` |
| `src/ranpanmixer` | the mechanism: HMM, utility, sampler, attacks, eval | `PLANNED` |
| `artifacts/releases` | sampled sanitized paths, one directory per (graph, tau) | `PLANNED` |
| `artifacts/eval` | utility and attack measurements | `PLANNED` |

## 2. Upstream dependencies

| source | what | state |
|---|---|---|
| `external/PanMixer` | G2Lab/PanMixer pinned at **`c182c38d5bc8bb6f00f4b0b101207c4a009ca045`** ("Merge pull request #1 from G2Lab/release_v1", 2026-06-06, MIT). 5 commits, ~8,660 Python lines. The only release; no tags, no test suite. | `OK` |
| conda env `panmixer` | Built from PanMixer's own `environment.yaml` at `/home/ds85/miniconda3/envs/panmixer`. Resolved: python 3.11.16, bcftools 1.24, htslib 1.24, PLINK v1.9.0-b.8, numpy 2.4.6, pandas 3.0.6, scipy 1.17.1, ortools 9.15.6755. ⚠ **The pip deps are unpinned upstream**, so these are 2026-09-18 resolutions, not the authors' versions; pandas 3.x and numpy 2.x post-date the paper. | `OK` |
| `paper/s41467-026-77591-0_reference.pdf` | The PUBLISHED PanMixer paper (Nature Communications, DOI 10.1038/s41467-026-77591-0). The version of record. | `OK` |
| `archive_docs/Blindenbach2026_PanMixer.pdf` | The PanMixer PREPRINT, bioRxiv DOI 10.64898/2026.02.16.706152, 24 pages, CC-BY-NC-ND 4.0, sha256 `99ac1bd2...`. **Superseded by the published version**; numbers differ. | `STALE` |
| Eizenga et al. 2020, *Pangenome graphs* | Background on graph construction and path representation. | reference |

## 3. PanMixer input data (inside the checkout, on `/data`)

Paths are written in full so the docs contract can check them. None are committed.

| path | what | state |
|---|---|---|
| `external/PanMixer/starting_data/pangenome.vcf.gz` | HPRC v1.0 PGGB GRCh38 deconstructed VCF, 5.70 GB, sha256 `ead25541...`. **45 samples, one of which is `chm13`** -> 44 individuals = 88 haplotypes after the pipeline drops it. Contigs are PanSN-named `grch38#chr1`..`grch38#chr22`, `grch38#chrX`. chr21: 340,824 unique (POS,REF,ALT), POS max 46,699,788. | `OK` |
| `external/PanMixer/starting_data/PG.vcf.gz` | PanGenie callset, Zenodo 7669083 `grch38_all-samples_bi_all`, 5.06 GB, sha256 `e8c0c48d...`, 368 samples. Source of the "anchor SNP" set. | `OK` |
| `external/PanMixer/starting_data/chr21/1000g_phased.vcf.gz` | ⚠ **Now a SYMLINK to the 30x GRCh38 panel.** PanMixer's code hardcodes this filename, so re-pointing the slot lets the pinned checkout run UNMODIFIED on the build-correct data. The original GRCh37 file is retained beside it as `DEPRECATED_1000g_phase3_grch37.vcf.gz`. | `OK` |
| `external/PanMixer/starting_data/chr21/DEPRECATED_1000g_phase3_grch37.vcf.gz` | 1000G **Phase 3** chr21, 0.21 GB, sha256 `1942e070...`, 2,504 samples. ⚠ **`assembly=b37`, `hs37d5.fa`, contig `21` length 48,129,895 — GRCh37.** This is what the released pipeline wires in. See the standing blocker. | `OK` |
| `external/PanMixer/starting_data/chr21/1000g_30x_phased.vcf.gz` | 1000G 30x **GRCh38** chr21, 0.43 GB, sha256 `a925c112...`, 3,202 samples, contig `chr21` length 46,709,983, 1,002,752 unique (POS,REF,ALT). The build-correct panel, from the URL in the repo's own unwired `get_blocks_grch38.sh`. **258,610 exact matches against the pangenome (75.88%) vs 1,025 (0.30%) for Phase 3.** | `OK` |
| `external/PanMixer/starting_data/pangenome_no_X.vcf.gz` | Autosome-only pangenome (pipeline step 2, `remove_X`). Being written by the local run; log `logs/pm_prep1.log` (its .status sentinel is written only on exit, so its absence means still running). | `BUILDING` |
| `external/PanMixer/starting_data/chr21/pangenome.npy` etc. | The preprocessed numpy model: `(n_subjects, n_sites, 2)` int16 allele codes, `-1` = missing, site axis row-aligned to the chr21 VCF record order. | `PLANNED` |
| `external/PanMixer/downloaded_tools/vg` | vg v1.68.0 binary, 49 MB. Needed for read-mapping utility only. | `OK` |
| `external/PanMixer/downloaded_tools/beagle.27Feb25.75f.jar` | Beagle 27Feb25, 0.3 MB. Needed for the reconstruction attack. ⚠ `environment.yaml` ships no `java`. | `OK` |
| `external/PanMixer/read_fastqs` | chr21 read slices for the paper's five external donors (HG00138 EUR, HG00635 EAS, HG01112 AMR, HG02698 SAS, NA18853 AFR), extracted from the 1000G 30x CRAMs. ~880 MB per donor paired, 4.3 GB total. ⚠ `constants.py` substitutes HG01600 for NA18853; we follow the PAPER. | `OK` |
| `external/PanMixer/starting_data/references` | GRCh38 analysis-set FASTA (for CRAM decoding), plus chr21 extracted from it and from UCSC hg38. ⚠ The two chr21 copies are the same length but the analysis set masks 2.05 Mb more (8,671,409 N vs 6,621,364 N); measured to make no difference to read mapping here. | `OK` |
| the 1000G Phase 3 panel for chr1-chr20, chr22 | Only chr21 was downloaded; the paper's tradeoff curves are all-autosome. | `PLANNED` |

⚠ **Not downloadable from this repo:** the five read FASTQs the mapping utility
needs, `chr21.fa` / `hg38_cleaned.fa`, and the three `.npy` files of the 30x
attack database. No script in the checkout produces any of them.

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
| `graph_id` | str | graph + cohort + block segmentation + **which leave-one-out cohort** |
| `seed` | int | the private random draw |
| `state_path` | int[T] | sampled Z_1:T |
| `log_z_p` | float | must satisfy 0 <= log_z_p <= eta_tau |
| `released_at` | date | one release per genome; re-releases compose |

**Integration contract.** Our sampler must ALSO emit PanMixer's artifact —
`new_haplotypes.npy`, shape `(n_sites, 2)` int16 with `-1` for missing, row-aligned
to the chr VCF — because every PanMixer evaluator consumes exactly that file.
Emitting it makes their whole attack and utility stack run on our output unchanged.

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
| `primer/primer.html` + `primer/src/` | 2026-09-20 | Teaching primer for onboarding: the landscape, the problem, PanMixer dissected, and our mechanism. ~57k words, built from `primer/src/` by `primer/build.py`. v2 expanded the hidden-Markov-model subsection (1.2.2) with full derivations of the forward algorithm and backward sampling, worked arithmetic, and a brute-force cross-check. URL: https://claude.ai/artifact/HzUWmUmxAPwXUSq3oFy8sb | LIVE |

## 8. Reproducibility contract

A result is reproducible here only if all five are recorded: the graph+cohort id
(**including which individual was held out**), `tau`, the utility id, the number
of draws, and the seed. Two results are comparable only when the first three match.

**The private input is the target path.** The `PLANNED` target directory must
never be committed, copied into an artifact directory, or quoted in any document.
Nothing in it is covered by the release guarantee; only the sampled output path
is. Written as a rule now, before the directory exists, because the first run is
where it would be broken.
