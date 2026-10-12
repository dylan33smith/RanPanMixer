# data.md — the map

Read before touching data, runs or paths.

**Rule.** Every number is traceable to a (`<ARTIFACT>`, `<GRAPH+COHORT>`,
`<TAU>`, `<UTILITY>`, n) tuple. **If it is not in this file, it is not quotable.**

**State values:** `OK` (present, current) · `PLANNED` (not built) · `BUILDING`
(a running job is writing it — the row must name its `logs/` sentinel; a
half-written file is NOT quotable) · `MISSING` (referenced but absent) · `STALE`
(present but superseded) · `DEPRECATED` (renamed, retained, never quote).

Verified against disk 2026-09-18. The PanMixer input table (§3), the `data/` tree
and `paper/` were re-verified against disk **2026-10-08** during the scope-change
sweep. `tests/test_docs_contract.py` re-verifies.

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
| `logs` | tmux run logs and status sentinels (gitignored). ⚠ Also where the scoring script writes its emitted tables and release records — see the artifact registry. | `OK` |
| `external/PanMixer` | symlink -> `/data/ds85/RanPanMixer/external/PanMixer`. **PanMixer writes all its data inside its own tree** (`BASE_PATH = dirname(constants.py)`, no env override), and `/home` is 95% full, so the checkout lives on `/data`. Gitignored. | `OK` |
| `external/genetic_maps` | symlink -> `/data/ds85/RanPanMixer/genetic_maps`. GRCh38 genetic maps, acquired 2026-09-27. Working file `external/genetic_maps/chr21.b38.gmap`. Gitignored. | `OK` |
| `src` | our code. Preprocessing BUILT 2026-09-29: `src/cohort/` (no target argument exists), `src/target/`, `src/sites.py`, `src/store.py`, `src/cli.py`. **The MECHANISM is BUILT 2026-09-30:** `src/mechanism/` — `hmm.py` (`baseline_model`; no target parameter, asserted from the AST), `utility.py`, `ffbs.py`, `calibration.py`, plus `toy.py` and `exact.py`, which read NOTHING from disk and are the falsification target. ⚠ **[CORRECTION - 2026-10-10]** the superseded text read "Attacks and evaluation are still `PLANNED`". `src/attacks.py` was BUILT 2026-10-10 (the Neyman-Pearson optimal statistic, a BALANCED `empirical_accuracy`, the arm guard, the pre-registration lock) with `src/runner.py` beside it; EVALUATION is still blocked, but on the SIGNATURE rather than on code. | `OK` |

### The build (consolidated into `main`, 2026-09-29)

The `pipeline` orphan branch and its separate `RanPanMixer-pipeline/` worktree were
folded into `main` on 2026-09-29 — one directory, one branch. The branch is retained
in git as history.

| path | contents | state |
|---|---|---|
| `src/` | the code. `src/cohort/` builds from the graph and cohort alone and **takes no target argument**; `src/target/` is strictly downstream. `src/sites.py` owns the site axis and its digest, `src/store.py` is the only module that writes an artifact | `OK` |
| `configs/` | the input manifest, three sample lists and a README | `OK` |
| `tests/tiny/` | a 24-record hand-written fixture, the only one committable | `OK` |
| `data/` | symlink -> `/data/ds85/RanPanMixer/data`, gitignored | `OK` |

**Artifacts not previously registered** (added 2026-09-29): `an_info.npy` (INFO/AN as the
caller wrote it, cross-checked against `support`), `allele_lengths_offsets.npy` (the ragged
index into `allele_lengths`), `target_reason.npy` and `readable.npy` (see `docs/plan.md`).
⚠ `data/temp/chr21/deconstruct.vcf` is UNCOMPRESSED, not `.vcf.gz`.

**Data tree.** `temp/` holds intermediates and is deletable; `cohort/` holds the
durable artifacts and no target ever touches it; `dev/` is a complete parallel root
reached with `--root data/dev`, so the 2 Mb slice runs the identical code path.

⚠ **INDEPENDENCE (2026-09-29).** Nothing in `src/`, `configs/` or `tests/` reads a path
inside PanMixer's checkout — verified by grep and by tracing every provenance sidecar. The
`vg` binary was relocated to `data/tools/vg` (it is third-party, v1.68.0, pinned by sha256 in
the manifest) and the 1000G panel to `data/raw/` (third-party, evaluation-only). PanMixer's
own chr21 VCF sha256 `5b880a3a…` appears in **no** artifact's inputs.

```
data/                                                      measured 2026-10-08
  tools/          vg v1.68.0 (49,244,936 B) and plink v1.9.0-b.8 (2,241,176 B),
                  both pinned by sha256                            51,486,153 B
  raw/            the staged inputs, plus chrY.hprc-v1.0-pggb.gfa.gz, which is a
                  DIAGNOSTIC input only -- decision D excludes chrY  882,917,881 B
  temp/chr21/     deconstruct.vcf (UNCOMPRESSED), cohort44.vcf.gz,
                  blocks_build/ (panel.bed alone is 699,392,739 B)
                  (deletable)                                   1,053,621,059 B
  cohort/chr21/   haplotypes reason missing_runs positions genetic_pos in_chain
                  support allele_lengths sites.tsv haplotype_ids manifest.json
                  an_info allele_lengths_offsets blocks_idx blocks.blocks.det
                                                                  246,146,173 B
  dev/            the same tree, built from the 14-16 Mb slice      51,531,838 B
  targets/<id>/chr21/   target.chr21.vcf.gz, path.npy, target_reason.npy,
                        readable.npy, representability.json      1,008,111,552 B
  attack/         1000g_30x_chr21_attackdb.vcf.gz  (evaluation only),
                  excluded_targets.txt -- a CONFIG that lives in the data tree
                  and is consumed by the attack-DB build (see configs/README.md)
                                                                  431,990,032 B
```

⚠ Paths inside the block above carry no backticks, so `test_every_referenced_path_exists`
never sees them. Four of them went unregistered for weeks for exactly that reason
(`data/tools/plink`, `data/temp/chr21/blocks_build/`, `data/attack/excluded_targets.txt`,
`data/raw/chrY.hprc-v1.0-pggb.gfa.gz`); all four are now in the block, and the chrY GFA has
its own row in §3. Add new paths HERE as well as in a row.

⚠ **chr21 is a TESTING scope, not the product.** Decision D makes the product the
**22 autosomes** — no chrX, no chrY — so the `chr21/` level of `temp/`, `cohort/` and
`targets/<id>/` becomes 22 chromosome directories. The sizes above are the chr21
figures; **the 22-autosome total is UNMEASURED** and must be measured, not scaled by
a factor carried over from a different scope, before the rerun is sized.
Preprocessing is already per-chromosome and parallel across chromosomes, so
extending it is a driver loop — but see `docs/plan.md` for the one thing that is
NOT just a loop: under decision B `beta_v` sums to 1 **genome-wide over the 22
autosomes**, each chromosome taking a SHARE, so a genome-wide `tau` means a
genome-wide normaliser, not 22 per-chromosome ones.

⚠ **Decision E adds no second tree per target.** Arm 1's leave-one-out is applied by
**masking the target's two haplotype rows in memory at release time**, not by
materialising 44 cohort trees — so `data/cohort/<chrom>/` stays single, the
`axis_digest` stays single (there is no second digest per target), and the thing that
records which rows were dropped is the **release record**, not the data tree. See §4
(`loo_masked_haplotype_ids`) and §8.

## 2. Upstream dependencies

| source | what | state |
|---|---|---|
| `external/PanMixer` | G2Lab/PanMixer pinned at **`c182c38d5bc8bb6f00f4b0b101207c4a009ca045`** ("Merge pull request #1 from G2Lab/release_v1", 2026-06-06, MIT). 5 commits, 8,661 Python lines (exact: `git ls-files` over `.py`; see `docs/memory.md` 2026-09-26). The only release; no tags, no test suite. | `OK` |
| `external/genetic_maps` | GRCh38 genetic map, the last unacquired input for `A-IMP-cohort-hmm`. TWO independent distributions of the SAME map: SHAPEIT4 `genetic_maps.b38.tar.gz` (github.com/odelaneau/shapeit4, pinned at commit `43a9a49703f15e37b4d20e703e0d2b04e0fa5cec`, 23,440,558 B, sha256 `04f97acc6524d75e1dbc397e72cf6776b2f1e33f72a06ee477ef69359a97c69e`) and Beagle `plink.GRCh38.map.zip` (Browning lab, 47,469,980 B, sha256 `521549889b9ce0236142a4fb7db45d3f00035ec645e01465490d55f5b8ef26d6`). **chr21 verified byte-identical between them** — 44,618 positions, identical cM, max abs diff 0 — so the Beagle copy is corroboration, not a second input. ⚠ PanMixer ships `external/PanMixer/starting_data/scripts/get_genetic_maps.sh`, which fetches a **GRCh37** map and is unwired. | `OK` |
| conda env `ranpanmixer` | **OURS** (`environment.yml`, created 2026-09-29). python 3.11.16, numpy 2.4.6, scipy 1.17.1, pyyaml 6.0.3, pytest 9.1.1, bcftools 1.24. Nothing taken from PanMixer's `environment.yaml`. `src/` imports only stdlib, numpy and pyyaml. `requirements.txt` is the venv equivalent. **Use this env, not `panmixer`.** | `OK` |
| conda env `panmixer` | Built from PanMixer's own `environment.yaml` at `/home/ds85/miniconda3/envs/panmixer`. Resolved: python 3.11.16, bcftools 1.24, htslib 1.24, PLINK v1.9.0-b.8, numpy 2.4.6, pandas 3.0.6, scipy 1.17.1, ortools 9.15.6755. ⚠ **`pyyaml` 6.0.3 added by us 2026-09-29** — absent upstream, and `src/cli.py` needs it to read the input manifest. ⚠ **`pytest` is NOT in this env**, so both test suites run under plain `python3` and `tests/test_preprocess.py` is self-running for that reason. ⚠ **The pip deps are unpinned upstream**, so these are 2026-09-18 resolutions, not the authors' versions; pandas 3.x and numpy 2.x post-date the paper. | `OK` |
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
| `external/PanMixer/starting_data/pangenome_no_X.vcf.gz` | Autosome-only pangenome (pipeline step 2, `remove_X`). **The job finished 2026-09-18** — `logs/pm_prep1.status` is `0` and `logs/pm_prep1.log` ends `### done Fri Sep 18 01:09:30 PM EDT 2026`; the file is 4,802,080,188 B (+ `.csi`). It sat marked `BUILDING` for three weeks after it stopped being written, which is the state column failing in the one direction the key cannot catch. ⚠ Marked `STALE` rather than `OK`: it derives from the PUBLISHED vg-1.36 `pangenome.vcf.gz`, and **decision F requires our own 22-chromosome vg-1.68 deconstruction**. Retained, like the `DEPRECATED` Phase 3 file, so the build deviation stays visible. | `STALE` |
| `external/PanMixer/starting_data/pangenome_no_X_no_chm13.vcf.gz` | ADDED 2026-10-08 — it was **unregistered and therefore unattributable**, 4,800,670,620 B (+ `.csi` 1,583,940 B), written 2026-09-18 by the same run. The 44-sample autosome-only VCF, i.e. the file an Arm 1 run would reach for first, which is exactly why leaving it off the registry was dangerous. ⚠ `STALE` for the same reason as the row above: vg-1.36-derived, and decision F requires our own build. | `STALE` |
| `external/PanMixer/starting_data/chr21/pangenome.npy` etc. | The preprocessed numpy model: `(n_subjects, n_sites, 2)` int16 allele codes, `-1` = missing, site axis row-aligned to the chr21 VCF record order. ⚠ **State corrected 2026-10-08:** this was marked `PLANNED` but the files exist (`pangenome.npy` 59,985,152 B, plus `pangenome_positions.npy` and `pangenome_subjects.npy`, all written 2026-09-18). It is `STALE`, not `OK`: row-aligned to PanMixer's vg-1.36 record order, which decision F replaces. The row escaped `test_data_md_layout_rows_marked_ok_exist_on_disk` because the ` etc.` after the closing backtick stops the row regex matching. | `STALE` |
| `external/PanMixer/downloaded_tools/vg` | vg v1.68.0 binary, 49 MB. Needed for read-mapping utility only. | `OK` |
| `external/PanMixer/downloaded_tools/beagle.27Feb25.75f.jar` | Beagle 27Feb25, 0.3 MB. Needed for the reconstruction attack. ⚠ `environment.yaml` ships no `java`. | `OK` |
| `external/PanMixer/read_fastqs` | chr21 read slices for the paper's five external donors (HG00138 EUR, HG00635 EAS, HG01112 AMR, HG02698 SAS, NA18853 AFR), extracted from the 1000G 30x CRAMs. ~880 MB per donor paired, 4.3 GB total. ⚠ `constants.py` substitutes HG01600 for NA18853; we follow the PAPER. | `OK` |
| `external/PanMixer/starting_data/references` | GRCh38 analysis-set FASTA (for CRAM decoding), plus chr21 extracted from it and from UCSC hg38. ⚠ The two chr21 copies are the same length but the analysis set masks 2.05 Mb more (8,671,409 N vs 6,621,364 N); measured to make no difference to read mapping here. | `OK` |
| the 1000G **Phase 3** panel for chr1-chr20, chr22 | ⚠ **WAS `PLANNED`; now `DEPRECATED` (2026-10-08).** This row planned to download 21 more chromosomes of the **GRCh37** Phase 3 panel — the panel whose every derivative is already `DEPRECATED` in §6 and whose chr21 copy sits on disk renamed `DEPRECATED_1000g_phase3_grch37.vcf.gz` (verified 2026-10-08). **Decision F — one consistent GRCh38 build — makes it exactly the wrong thing to download.** Never acquire it. The panel actually needed for 21 more chromosomes is the 30x GRCh38 one, in the row below. | `DEPRECATED` |
| the 1000G **30x GRCh38** panel for chr1-chr20, chr22 | ADDED 2026-10-08 under decisions D and F: the all-autosome build needs the build-correct panel for the other 21 autosomes, for the EVALUATION roles only (target pool, attack database) — never as a mechanism input. Same URL pattern as chr21, with the chromosome substituted: `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.<chrom>.filtered.SNV_INDEL_SV_phased_panel.vcf.gz` (the chr21 URL is pinned in `configs/inputs.chr21.yaml`). Size: chr21 is **427,136,458 B**; the 21-chromosome total is **UNMEASURED** — read the EBI directory listing before committing disk, do not scale the chr21 figure. Stage under `data/raw/`, never read from PanMixer's tree. ⚠ **[CORRECTION - 2026-10-10] THIS ACQUISITION PLAN IS FULFILLED** — all 21 are downloaded and verified; the artifact row above carries the measured figures (29.3 GB, 69,689,262 records, 3,202 samples) and supersedes this row's forecast. The "UNMEASURED" warning stands only as the record of what was unknown beforehand. | `DONE` |
| the other 21 **HPRC v1.0 PGGB** chromosome GFAs (chr1-chr20, chr22) | ✅ **DOWNLOADED 2026-10-10** by `scripts/get_gfas.sh`. All 22 autosomes now on disk: **22 files, 15.4 GB**, every one with a **sha256 recorded in a `.sha256` sidecar** — a graph with no recorded checksum is not quotable and re-deriving one means re-reading 15 GB. Zero failures; the script asserts its own completeness over all 22 and exits non-zero otherwise. chrY is NOT among them (decision D). | `OK` |
| the per-chromosome **genetic maps** for chr1-chr20, chr22 | ✅ **EXTRACTED 2026-10-10** (chunk 3). All 21 remaining autosomes pulled from `external/genetic_maps/genetic_maps.b38.tar.gz`, which was already on disk — this was never a download. ⚠ **chrX was deliberately NOT extracted** (decision D, autosomes only), and a test asserts it stays that way, because an extracted chrX is what a later `range(1, 24)` picks up without anyone deciding to. Each map is asserted strictly position-increasing and non-decreasing in cM at read time, and the per-autosome coverage is measured in `logs/B_map_coverage.md`. | `OK` |
| `data/raw/chrY.hprc-v1.0-pggb.gfa.gz` | ADDED 2026-10-08 — **it was on disk since 2026-09-29 in no row and no configs manifest, so by this file's own rule it was not quotable, while two numbers quoted all over the repo were measured on it.** 20,000,273 B, downloaded 2026-09-29 (`logs/chry_dl.log`, status `0`) and deconstructed (`logs/chry_deconstruct.status` = `0`). Its two roles, both DIAGNOSTIC: (1) the vg-1.68-vs-1.36 site-agreement measurement — 99.55% of the published VCF's sites by POS, **92.7% by (POS, REF, ALT)** — quoted in §3b and in `configs/inputs.chr21.yaml`; (2) the 34,014-record smoke test that vg 1.68 groups PanSN samples correctly, cited in `src/cohort/vcf.py`. ⚠ **Decision D excludes chrY from the product** — autosomes only — so this will never be promoted to a mechanism input and no release ever covers it. Kept `OK` rather than `DEPRECATED` deliberately: `DEPRECATED` means *never quote*, and the two figures above are live and must stay traceable to this file. | `OK` |
| `external/genetic_maps/chr21.b38.gmap` | chr21 extracted, 971,965 B, columns pos/chr/cM with a header. Span **10,326,676 - 46,680,243 bp**, **0.5841 - 62.7865 cM**, monotonic, no duplicate positions. ⚠ **26,018 chr21 variants (7.63%) fall BELOW the map start** and 592 (0.17%) above. The denominator for those two counts is PanMixer's vg-1.36 **340,824**; OUR build has **340,849** records and neither count has been re-measured on it — the 2026-09-29 recount covered the `chain_span` cut (`dropped_below` 34,371, `dropped_above` 592), not the map-span edges. ⚠ **POLICY CHANGED — decision C (2026-10-08).** These map-silent regions, where the map gives **no cM at all**, are no longer excluded: they are sampled under a **declared constant recombination rate**. The old description — interpolation clamps, cM is flat, the sampler cannot recombine — describes the SUPERSEDED policy. **The rate's numeric value is NOT yet chosen** and must not be invented (still-open question 1 in `docs/plan.md`; the per-autosome map-coverage measurement it rests on does not exist yet). Separately, 618 interior intervals inside the mapped span (1.4%, 2.71 Mb, chr21 only) have zero cM change; **decision C does NOT override those** — `P(switch) = 0` between two variants with no genetic separation remains CORRECT, and the 22-autosome figure for them is UNMEASURED. ⚠ **One residual needs Dylan and is flagged, not resolved**, in the `chain_span` block of `configs/inputs.chr21.yaml`: the superseded cut at 12,968,320 reached ABOVE the map's first marker at 10,326,676, so part of what it dropped is an **in-span plateau** that the declared rate does NOT reach and the plateau rule says to leave — the split of those positions into map-silent and in-span is UNMEASURED. | `OK` |

⚠ **Not downloadable from this repo:** the five read FASTQs the mapping utility
needs, `chr21.fa` / `hg38_cleaned.fa`, and the three `.npy` files of the 30x
attack database. No script in the checkout produces any of them.

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

⚠ **All three rows are chr21-only, and decision D makes the product the 22 autosomes.**
Inputs 1 and 2 therefore need 21 more instances each; both are registered as `PLANNED`
rows in §3 (the GFAs are downloads, the maps are an extraction from a tarball already on
disk). Input 3 is arm-dependent and the row above describes **Arm 2's** stand-in only: in
**Arm 1 the target is one of the 44 cohort members** and no external sample is drawn at
all, so the eligibility rules and pools below apply to Arm 2 alone.

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
the published VCF carries 45 samples, we need 44. ⚠ **REVISED 2026-10-08: this is no longer
conditional.** The head-to-head IS Arm 1, and decision F requires PanMixer's starting data
to be **rebuilt on our GRCh38 build** so that **both arms run on ONE consistent VCF** — the
thing upstream's own release got wrong. Every `STALE` row in §3 is stale for this reason.

⚠ **The graph is published only as GFA for PGGB freeze1** — no `.gbz`, no `.snarls`. The two
`.gbz` files on disk under `/data/ds85/RanPanMixer/runs/` are outputs of our own
`vg autoindex` read-mapping baseline, built from the VCF plus `chr21_ucsc.fa`, so they carry
no information the VCF does not.

## 4. Record schema — released path

**LIVE.** ⚠ Corrected 2026-10-10: this read "**PLANNED.** Fields fixed now so releases
are self-describing from the first run." Releases exist and carry these fields —
`logs/A_HG00096_chr21_v2.releases.json` (2026-09-30, one record per release) is
registered LIVE in section 7. What is still PLANNED is routing them through
`store.save()`: today they are written by a bare `write_text` from
`scripts/run_mechanism_chr21.py`, with no `.prov.json`, no `git_rev` and no
`axis_digest`, so by this file's own rule the released paths are not yet quotable
artifacts. ⚠ And the genome-wide fields are not fixed: a release spanning 22
autosomes needs a summed `log_z_p`, a weighted `u_achieved`, an arm, a `loo_id` and
per-chromosome digests, none of which this schema carries (`release_id` currently
collides across all 22 chromosomes of one release).
A release that cannot answer "at what tau, under what utility" is not quotable.

| field | type | notes |
|---|---|---|
| `release_id` | str | ⚠ **not a unique key on its own.** The unique key is (`arm`, `target_id`, `loo_id`, `tau`, `utility_id`, `seed`). Collisions must raise, never overwrite. ⚠ **EXTENDED 2026-10-08:** without `arm` and `loo_id` an Arm 1 release of a cohort member could not be distinguished from an Arm 2 release at the same `tau` and utility, and the schema says collisions raise — so the two arms would have collided. |
| `arm` | str | ADDED 2026-10-08. `"arm1"` or `"arm2"`, and **never absent** — the arm decides what is true of the target, so a record without it is uninterpretable, not merely imprecise. Arm 1: the target IS one of the 44 cohort members, released under leave-one-out at K = 86. Arm 2: the target is outside `G`, `D`, every panel and the attack database, full cohort as panel, K = 88, nothing held out. |
| `target_id` | str | pseudonymous target identifier |
| `tau` | float | the privacy parameter. Absent = the record is unusable. |
| `eta_tau` | float | stored, not recomputed at read time, so a calibration change is detectable |
| `utility_id` | str | identifies phi_t, its weights, and the beta_t scheme |
| `graph_id` | str | graph + **cohort composition** + chain extent. ⚠ **In Arm 1 this IS a leave-one-out id** — the cohort differs per target, so `graph_id` carries the LOO identity; in Arm 2 nothing is held out. ⚠ The chain extent must name the **chromosome SET**, not one chromosome: the emitted form on record (`chr21/hprc-v1.0-pggb/K88/T305886`) silently claims a per-chromosome release, which decision B forbids — a release is genome-wide over the 22 autosomes. |
| `loo_id` | str/null | the held-out individual. ⚠ **MANDATORY and non-null when `arm == "arm1"`**, `null` when `arm == "arm2"`. Decision E makes leave-one-out a defining property of Arm 1, not an option: it is half the paper, not a "comparison arm". |
| `loo_masked_haplotype_ids` | str[2]/null | ADDED 2026-10-08, and load-bearing. The construction itself is `loo_cohort` in `docs/terms.md`; this field is its provenance in the record. Decision E applies leave-one-out by **masking two haplotype rows in memory at release time**, not by materialising 44 cohort trees (identical artifacts, ~160-380 GB cheaper). Consequence: **there is no second `axis_digest` per target** and nothing in `data/cohort/` records which rows were dropped, so the only place the provenance can live is here. Both of the target's haplotype ids, exactly as they appear in `haplotype_ids.npy`. `null` in Arm 2. A record without this is not reproducible: the K = 86 panel cannot be reconstructed from the artifacts alone. |
| `seed` | int | the private random draw |
| `state_path` | int[T] or int[T,2] | sampled Z_1:T; (T,2) when `ploidy` = 2 |
| `released_path` | int16[T] or int16[T,2] | the alleles actually released. The quotable object |
| `log_z_p` | float | must satisfy 0 <= log_z_p <= eta_tau. ⚠ For a diploid release it is the SUM of the two haplotype values, each in [0, eta_tau/2] — so the bracket holds for the release as a whole, at the real `tau` |
| `u_achieved` | float | `utility_retained` for this draw. Meaningless without the next field |
| `u_ceiling` | float | ADDED 2026-09-30. The largest `u` this target could reach at any `tau` under this `beta` scope. **0.692509** on HG00096 chr21 at scope `chain`, not 1 — ⚠ **[CORRECTION - 2026-10-12]** this read **0.815455**, which was measured on the superseded T = 305,886 coordinate-cut chain; the live decision-K chain is T = 332,521 and 30.749104% of its weight is on positions the target does not call |
| `ploidy` | int | ADDED 2026-09-30. 1 or 2. A diploid target's two haplotypes are ONE private input and get ONE release — one `tau`, one `log_z_p`, one bracket |
| `released_at` | date | one release per genome; re-releases compose |

**Integration contract (ARM 1).** Our sampler MAY also emit PanMixer's
artifact — `new_haplotypes.npy`, shape `(n_sites, 2)` int16 with `-1` for missing,
row-aligned to the chr VCF — because every PanMixer evaluator consumes exactly that file.
⚠ It is NOT a requirement of the mechanism. **What transfers depends on the arm, and the
two cases are opposite** — never state one as if it covered both.

**Arm 2 (external target).** `af_loss`, `ld_loss` and the Beagle attack will run and
return numbers that mean nothing for a release which never edits the cohort. Only the
gap-score attack transfers, and only after its "self" becomes the external target's true
path and the target is excluded from the attack database.

**Arm 1 (the target IS a cohort member).** Both clauses above invert. The release edits a
cohort column, so `af_loss` and `ld_loss` ARE defined here. ⚠ **We still report
`utility_loss_panmixer` only (decision H) — but the reason is a SCOPE CHOICE, not a
definitional impossibility.** Any text saying they are omitted "because we never modify
the cohort" is false in Arm 1 and must not be reinstated. On the attack side:
`diploid_gap_score`, `MIA_privacy`, `sampled_gap_score` and the Beagle reconstruction
attack all become meaningful once a cohort column is edited, so "only one attack
transfers" is false here. The target is a cohort member, so **excluding it from the attack
database would remove the thing being attacked**; the database retains it under the Arm 1
rule, and `self` is scored from **the member's own assembly-derived path** (defined for all
44), not its 1000G panel row (defined for only 39 of 44) — decision I; PanMixer's
panel-row semantics run as a labelled parity column. The exact wording of the retention
rule is fixed in the pre-registration — ⚠ **[CORRECTION - 2026-10-10]**
**`paper/prereg_2026-10-10.md` EXISTS** and its checklist item 7 settles this rule as
TWO variants (V1 conservative carries the claim, V2 permissive reported beside it). The
superseded text (SUPERSEDED) read: "`paper/prereg_<date>.md`, which is `PLANNED` and
does not exist yet ... cite it as PLANNED, never as settled text". ✅ **SIGNED 2026-10-10** (Dylan), `git_rev` b5a220f122ca, so it IS in force — cite it as
pre-registered. The fields are frozen under constraint 6 from that date.
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

**Status key:** LIVE (quotable) · SUPERSEDED · INVALID (never quote) · PLANNED (not
built; not quotable, and every reference to it must say PLANNED)

| artifact | date | contents | status |
|---|---|---|---|
| `paper/addendum_2026-10-12_chain_rebuild.md` | 2026-10-12 | **ADDENDUM 2 to the signed pre-registration, SIGNED 2026-10-12.** Corrects two DESCRIPTIVE figures the signed file quotes at `:94-95` and that were measured on the superseded T = 305,886 coordinate-cut chain: the share of `beta_v` the target cannot call is **30.749104%**, not ~18.5%, and `u_ceiling` is **0.692509**, not ~0.815 (HG00096, chr21, `beta_scope` = `chain`, ARM 2, live T = 332,521). ⚠⚠ **NO PRE-REGISTERED FIELD CHANGES** — `tau`, the utility function, the endpoints, the attacker set, the database and the decoy design are untouched, and R1's disposition (declined on deployment-realism grounds) stands; only the size of the accepted cost moved. ⚠ The cause is demonstrated rather than asserted: restricting the live arrays to the old span reproduces T = 305,886 and 0.815455 EXACTLY, so the axis moved and the arithmetic did not. ⚠⚠ States explicitly that no privacy claim moves, because `delta_u` = 1 by construction and `u_ceiling` is a utility DENOMINATOR that must never enter the calibration — doing so is local-sensitivity calibration and the safe form gains exactly zero, since sup `u_ceiling` = 1. ⚠ Deliberately NOT named `prereg_*.md`; verified that `attacks.require_signed_prereg` still resolves to the pre-registration. | LIVE |
| `paper/addendum_2026-10-11_decision_K.md` | 2026-10-11 | **ADDENDUM 1 to the signed pre-registration, SIGNED 2026-10-11.** Decision K post-dated the signature and changed the released object: 172,142 positions left the chain, so genome-wide `T` went 22,306,704 -> **22,134,562**, `mean beta_v` 4.4830e-08 -> **4.5178e-08**, and every `beta_share` moved (chr9 **-6.43%**). It corrects TWO things in the signed file without touching it: section 7 item 6's claim that the interior plateaus are "emitted as verbatim donor runs" (20 blocks are now EXCLUDED), and the fact that `beta_scope` = "chain" now covers a different set of positions under the same name. ⚠ The frozen fields are untouched — `tau` grid, attacker set, endpoints, attack database, decoy design, `delta` — and no attack outcome had been inspected. ⚠⚠ **DELIBERATELY NOT NAMED `prereg_*.md`**: `require_signed_prereg` globs that pattern and returns the first signed match, so a second matching file could shadow the real one if it were ever renamed or unsigned. Three tests pin that, including a structural one that fails the moment anything else matches the glob. | LIVE |
| `paper/prereg_2026-10-10.md` | 2026-10-10 | **The pre-registration: ONE dated file covering BOTH arms.** Fixes `tau`, the utility function, the endpoints, the attacker set and the attack database before any attack outcome is inspected (constraint 6), plus the Arm 1 attack-database retention rule and the `self` semantics of decision I. ⚠ **[CORRECTION - 2026-10-10] IT EXISTS: `paper/prereg_2026-10-10.md`** (23,921 B), seven of eight fields filled, field 6 `[OPEN]` again since 2026-10-10, and **UNSIGNED**. The superseded text read "IT DOES NOT EXIST YET — verified against `paper/` 2026-10-08, which holds only the two PDFs". ⚠ Nothing may still be quoted as "pre-registered" — not because the file is absent, but because it is not in force until its status line reads `SIGNED` with a `git_rev`. ⚠ **[CORRECTION - 2026-10-10]** references to it must now be marked **UNSIGNED**, not PLANNED — the file is dated and committed, so "PLANNED" would understate what exists and misdirect anyone looking for it. The superseded instruction read "Any reference to it, anywhere in the repo, must be marked PLANNED; it is the work item `B-LCK-preregister` ... and it is blocked". The work item is no longer blocked on WRITING; it is blocked on SIGNING, and on field 6, which was reopened 2026-10-10 because the pinned attack design was one-sided. | PLANNED |
| `paper/Private_Genome_Path_Release.pdf` | 2026-08-29 | The proposal: mechanism, Theorem 1, Corollary 1, Algorithm 1, runtime analysis. | LIVE |
| `archive_docs/The_Six_File_Lab_Record.pdf` | 2026-08-29 | The documentation framework this repository follows. | LIVE |
| `paper/s41467-026-77591-0_reference.pdf` | 2026-09-18 | **The PUBLISHED PanMixer paper** (Nature Communications, DOI 10.1038/s41467-026-77591-0, Article in Press, 12 pages, sha256 `2bbba736...`), supplied by Dylan. **The version of record — all reproduction targets come from here.** | LIVE |
| `archive_docs/Blindenbach2026_PanMixer.pdf` | 2026-09-18 | The bioRxiv PREPRINT, 24 pages. Superseded: it says 47 individuals (published: 44), `eps_private` 0.001 (published: 0.002), AF divergence 0.004/0.004/0.002 (published: 0.006/0.006/0.005), and it has no Data availability section and no membership-inference analysis. Retained because this project's first analysis was made against it. | SUPERSEDED |
| `logs/pm_download.sha256` | 2026-09-18 | sha256 of the three downloaded VCFs. | LIVE |
| `data/cohort/DEPRECATED_chr21_coordinate_cut/blocks_idx.npy` + `blocks.blocks.det` | 2026-10-07 | **OUR PLINK LD blocks, the model-B chain unit.** Built by `scripts/build_blocks_chr21.sh` with PLINK v1.9.0-b.8 (staged at `data/tools/plink`, sha256 `4674107b...`) on our 30x panel with **the 10 test targets REMOVED** (3,192 samples, SNPs only) — 8,773 called haploblocks, longest span 172.7 kb. Assigned to the site axis by `scripts/assign_blocks_chr21.py`; orphans become their own singleton units, never merged. ⚠ PanMixer's own block set used all 3,202 samples, which CONTAINS our targets, so it is weakly target-dependent and is NOT used. ⚠ Residual, same shape as R4: first-degree relatives of a target are not excluded. OFF THE MAIN PATH — model C does not read this. ⚠⚠ **MOVED AND SUPERSEDED 2026-10-10.** It was assigned to the PRE-decision-C site axis, so it carries the OLD `axis_digest` (`251b554f...`) and would have flipped chr21's `axis_consistent` to false had it been left beside the rebuilt arrays. A model-B unit set on the current axis must be rebuilt from scratch; this one cannot be reused. | SUPERSEDED |
| `logs/A_mismatch_map_chr21.md` + `.json` | 2026-10-08 | **Whole-chromosome mismatch density per `tau`, binned.** Emitted by `scripts/mismatch_map_chr21.py`. All 304,664 chain positions in 540 bins of 564 variants, density PER CELL so it averages to the share quoted elsewhere. Headline: regional density spans **0.0% to 40.7%** while the largest shift `tau` causes anywhere is **0.0077**, and flips touch **4 of 540 bins**. ⚠ Replaces a 200-variant window figure that covered 0.07% of the chromosome and was chosen for maximum churn, so its local density (~20%) was unrepresentative of the 14.5% chromosome-wide. | LIVE |
| `paper/prereg_2026-10-10.md` (duplicate row) | 2026-10-10 | ⚠ **DUPLICATE — see the fuller row above for this same file.** Two rows for one artifact is two places for its state to drift, and it did: this one still read "DRAFT, UNSIGNED ... 8 fields OPEN" after all eight were filled and the file was signed. ✅ **SIGNED 2026-10-10** (Dylan), `git_rev` b5a220f122ca. | SUPERSEDED BY THE ROW ABOVE |
| `logs/B_toy_genome_demo.md` + `.json` | 2026-10-10 | **The whole mechanism on a two-chromosome toy genome, enumerated exactly.** Emitted by `scripts/toy_genome_demo.py`. chrA = the main toy (96 output paths), chrB = the small segment (4), shares 7/9 and 2/9 by **decision J**. Mean `beta_v` = 0.111 against chr21's 3.27e-06, about 34,000x larger, which is why the tilt is VISIBLE here: `E_Q[u]` moves 0.4502 -> 0.6477 across the grid against chr21's 2.94e-06. `ratio_spread` equals `exp(epsilon_tau)` exactly at every `tau`. The sampler matches the enumeration within 4 sigma at 40,000 draws, and `log_z_p` matches to 9 decimals. ⚠ A DEMONSTRATION, not a gate — `tests/test_diploid_gate.py` is the gate. | LIVE |
| `src/attacks.py` + `tests/test_attacks.py` | 2026-10-10 | **The arm guard and the optimal attacker.** `assert_well_posed` REFUSES `db_search` and `membership` in arm 2, where the first has no correct answer (the target is outside the database by construction, so 0% is the experimental design) and the second inverts (a donor mosaic matches cohort members, identifying DONORS). `AttackSpec` makes `positive_control` a mandatory field, enforcing standing constraint 7 at construction. `lrt_statistic` is the Neyman-Pearson optimum from the proposal's own algebra; `exact_advantage` is the ceiling no attacker can exceed. | LIVE |
| the 21 generated per-chromosome manifests under `configs/` | 2026-10-10 | **GENERATED** by `scripts/make_chrom_configs.py` from one template, with checksums read from the `.sha256` sidecars rather than hand-copied. ⚠ Decision C is what makes a generator viable: the per-chromosome `chain_span` was the only genuinely per-chromosome parameter, so 22 hand-written files would now be 22 places for a shared value to drift silently. The generator REFUSES to write a partial set. ⚠ `configs/inputs.chr21.yaml` is NOT generated — it carries the superseded `chain_span` block and the reasoning behind decision C. | LIVE |
| the 9 orphaned deconstruct VCFs, one per chromosome under `data/temp/chr10/deconstruct.vcf` and its siblings (chr10-15, chr17, chr21, chr22) | 2026-10-10 | **COMPLETE vg OUTPUT WITH NO PROVENANCE SIDECAR**, left by the `8d56346` ordering bug — vg wrote each VCF and the run then died attaching its sidecar. The BYTES are vg's own and unmodified; what is missing is the record. ⚠ NOT QUOTABLE until `vcf.attach_deconstruct_prov` has run on them, which recomputes the GFA digest from the `.gz` without re-running vg. Any sidecar it writes carries `prov_recovered_after_crash: true`, so a reconstructed lineage is never mistaken for one recorded at the time. | NEEDS PROV |
| the 21 **1000G 30x GRCh38 phased panel** VCFs, `data/raw/1000g_30x_chr1.vcf.gz` and its siblings (chr1-chr20, chr22) | 2026-10-10, chr10 refetched 2026-10-11 | ✅ **DOWNLOADED**, sentinel 0. ⚠⚠ **chr10 ARRIVED CORRUPT** — a CRC32 mismatch in the bgzf block at offset 393,732,438, which passed every check in `scripts/verify_panel.py` including the sha256, because the sidecar is computed FROM the download and the truncation test only SEEKS to the far end. Found by `bcftools` three steps downstream. Refetched and verified with `bgzip -t`; all 22 panels now swept, 21 originally sound. The corrupt file is kept as `data/raw/CORRUPT_chr10.vcf.gz.bad` as evidence. 29.3 GB, **69,689,262 records**, **3,202 samples on every chromosome**, each with a `.tbi` and a `.sha256` sidecar. Fetched in parallel by `scripts/get_panel.sh` (4 connections; a single one measured ~700 KB/s, which would have been 11 h) with `curl -C -` so a dropped connection resumes rather than restarts. ⚠ chr21 is NOT among them — it was already on disk as `data/raw/1000g_30x_phased.vcf.gz` (0.43 GB, 1,002,753 records, 3,202 samples), so the 22-chromosome panel is complete. ⚠ It serves TWO purposes that must not be conflated: Arm 2's external TARGET source and the ATTACK DATABASE; the eligibility rules that separate them are decision-gated and pre-registered. ⚠ Decision F: this is the GRCh38-called 30x panel, NOT Phase 3 — PanMixer's released pipeline wires in GRCh37 Phase 3 against a GRCh38 graph, and 0.30% of chr21 pangenome records match Phase 3 on (POS, REF, ALT) against 75.88% here, a 252x difference attributable to the build alone. | OK |
| `scripts/verify_panel.py` + `tests/test_verify_panel.py` + `logs/B_verify_panel.md` | 2026-10-10 | **THE PANEL GATE, and it checks the REFERENCE BUILD against our own pangenome rather than against a typed-in constant.** Each panel VCF's `##contig` length must equal the length OUR `vg deconstruct` output declares for the same chromosome — two independent files on disk, both describing the assembly they were built on. ⚠⚠ This is the decision-F gate, and it matters because a wrong-build panel does NOT fail loudly: it produces a low match rate that reads as biology. A GRCh37 panel would disagree on every chromosome and no constant in the script would need to be right for that to be caught. Also checks truncation by SEEKING to the far end (a cut-short download still opens and still indexes against its own partial content), the three required files per chromosome, a leftover `.part`, and one shared sample set across all 21. **21 of 21 pass; the build was verified on all 21.** 9 tests; deleting the build comparison fails exactly its two tests. ⚠ `--deep` re-reads 29 GB to verify every sha256 and was NOT run — the sidecars are recorded, not re-verified. | LIVE |
| `src/releases.py` + `tests/test_releases.py` + `data/releases/` | 2026-10-10 | **THE RELEASE WRITER — the chokepoint that made a release quotable at all.** `src/mechanism/ffbs.py` never imported `src/store.py`, so a release lived only in memory and the one script that persisted releases used a bare `json.dumps` — the first release on record has no sidecar, no `git_rev` and no axis digest. Layout `data/releases/<arm>/<target>/tau_0p10/<release_id>/` holding `state_path.npy`, `released_path.npy` and `release.json`, each with provenance. ⚠ Seven gates, all mutation-tested: a COLLISION on the key (target, `tau`, utility, seed, chrom) RAISES rather than overwriting; the ARM is explicit and must agree with `loo_id` (arm 1 requires leave-one-out, arm 2 forbids it); the AXIS DIGEST is required; a `tau` that would not survive its two-decimal directory name is refused (0.995 and 0.99 would collide); and **genome scope is asserted against the 22-autosome set** — ⚠⚠ a `GenomeRelease` of ONE chromosome with share 1.0 passes every in-memory check it makes, so only the write can tell. Standing constraint 5 is made VISIBLE: every sidecar records the sibling seed count, and `src/releases.py:assert_one_release_per_genome` is what anything quoting a release as DEPLOYED must call. 20 tests; removing any of the seven gates fails its own test. ⚠ Verified end-to-end on the pre-decision-C chr21 tree, because the REBUILT tree currently dead-ends the forward pass — see the blocker in `docs/plan.md`. | LIVE |
| `scripts/rebuild_chain.py` + `logs/B_rebuild_chain.log` | 2026-10-10 | **THE DECISION-K CHAIN REBUILD, 22 autosomes in parallel.** Re-runs ONLY `genmap` (forced) and `manifest`: decision K changes `in_chain.npy` and nothing upstream, so `deconstruct`, `cohort-vcf`, `arrays` and `support` are untouched and hours of work is not redone. ⚠ It ASSERTS THE ARITHMETIC CLOSES — `before - dropped == after` on every chromosome — because a silent chain rebuild is the most dangerous edit available here: every denominator in the project counts chain positions, and `beta_share` is allocated by variant count under decision J, so a changed `T` moves the budget split on every OTHER chromosome too. Result: **20 blocks, 172,142 positions, 0.772%**, exit 0. | LIVE |
| `src/cohort/genmap.py:untraversable_blocks` + `tests/test_untraversable.py` | 2026-10-10 | **DECISION K's implementation and its 15 tests.** Fixtures are hand-checkable on purpose: the real chromosomes have 52,870 zero-distance runs of which only 20 are untraversable, so a rule that dropped too much or too little would be invisible against real data. Three mutations caught: removing the `+1` that makes a zero STEP span TWO positions, checking adjacent PAIRS instead of the run's intersection, and dropping every zero-distance run rather than only the untraversable ones. ⚠ One test asserts from the function signatures that the rule takes no target-shaped argument — standing constraint 1 in the type, not in a comment. | LIVE |
| the 220 target artifacts under `data/targets/HG00096/chr1/path.npy` and its siblings (10 targets x 22 autosomes) | 2026-10-11 | **EVERY PREPARED TARGET ON EVERY AUTOSOME.** Four artifacts per pair — `path.npy`, `readable.npy`, `target_reason.npy`, `representability.json` — each with provenance. **187,268,720 called cells** over a genome-wide `in_chain` of **22,134,562**, mean called fraction **0.846047**, by chromosome 0.7141 (chr22) to 0.8778 (chr3). ⚠ `path.npy` is FULL-AXIS, not chain-masked, which is why decision K's chain change did not invalidate the existing chr21 paths and why a future chain revision will not either. ⚠ Targets differ from each other at 15-21% of chain positions, with same-population pairs closer than cross-population ones (HG00096/HG00097 84.50%, HG00096/NA20294 79.14%) — the check that rules out ten copies of one person, which a presence check cannot. ⚠ chr21's ten pairs carry NO `extracted_by` provenance: they were hand-made before the producer existed and resume correctly skips them. | LIVE |
| `scripts/prep_all_targets.py` + `logs/B_prep_all_targets.log` | 2026-10-11 | **THE PRODUCER THE TARGET VCFs NEVER HAD.** chr21's were made by hand with no sidecar, so the only target set that existed was irreproducible. This records the `bcftools view -s` command and the panel digest on an artifact that SURVIVES, then deletes the single-sample VCF (~6.8 GB per target genome-wide, 68 GB for ten) as a pure extract of a checksummed panel. ⚠ The digest is taken BEFORE the unlink — the ordering whose reverse cost nine chromosomes on 2026-10-10. ⚠ Resume treats a target VCF present WITHOUT its `path.npy` as evidence of an interrupted run and DISCARDS it: a half-written bgzf is non-empty, so an `st_size == 0` test would have fed a truncated VCF to the next stage. Exercised for real on the 10 chr10 pairs. | LIVE |
| `logs/B_beta_concentration.md` + `logs/B_beta_concentration.json` | 2026-10-12 | **WHAT CONCENTRATING THE BUDGET BUYS, EXACTLY — the answer to Dylan's accumulated-weight question.** Emitted by `scripts/beta_concentration.py`, 14 allocations x 3 `tau`, exact marginals, chr21, `beta_share` = 1 (CHR21-SCOPE, NOT genome-wide — 66.6x more generous per position than the real release), HG00096, ARM 2, n = 1. Across `tau` = 0 to 0.99 the production allocation moves `e_q_u` by +3.6368e-06; **10 CONTIGUOUS positions move it by +1.6546e-02, about 4,500x**. Contiguity is worth 6.7x at 10 positions and 60x at 1,000 against evenly-spaced weight at matched count and matched `sum beta_v`. ⚠⚠ **PRIVACY IS UNCHANGED ON EVERY ROW** — each allocation sums to 1 so `delta_u` = 1, and `tau`, `eta_tau`, `epsilon_tau` and `p_succ_bound` are identical throughout. Reallocation is privacy-neutral; it redefines what is called useful. ⚠⚠ **READ THE HEADROOM COLUMN, NOT THE MOVEMENT**: `window_1` captures 73.3% of its headroom but that headroom is 0.004164, because the baseline already matches that position 99.6% of the time — the positions where the tilt has power are mostly the ones that needed it least. ⚠ Every allocation is selected by POSITION INDEX and so cannot see the target; `spread_1` returns identically zero because that index is a no-call for HG00096, which is the honest cost of target-independent selection. ⚠ Constraint 8: each row is a DIFFERENT utility function, so levels are not comparable across rows. ⚠ Carries a dated CORRECTION block — the tilt column first reported full `eta_tau` where the diploid release applies `eta_tau`/2, a 2x overstatement of the realised exponent; regenerated from the same JSON, measured values bit-identical. | LIVE |
| `logs/B_exact_utility.md` + `logs/B_exact_utility.json` | 2026-10-12 | **THE GENOME-WIDE UTILITY FIGURE, EXACT — and it is the result, not a defect.** Emitted by `scripts/exact_utility.py` via the state marginals, 44 exact passes (22 autosomes x 2 `tau`), ~30 min, zero Monte-Carlo error and no finite difference. HG00096, ARM 2, scope "chain", shares by decision J summing to **1.000000** over T = 22,134,562, n = 1 target. `e_q_u` = **0.708967646** at `tau` = 0 and **0.708967729** at `tau` = 0.99; `u_ceiling` = 0.825708741; `e_q_u_over_ceiling` = 0.858617102 and 0.858617203. **THE MOVEMENT IS +8.345e-08** (+0.000012%) and that is the only figure that carries a claim. ⚠⚠ **DO NOT QUOTE THE 0.8586 AS "WE RETAIN 86% OF UTILITY"** — it is reached at `tau` = 0, where the release is a pure cohort mosaic carrying ZERO target information, so the level is not attributable to the mechanism; that is the free-lunch error `target_fidelity` exists to prevent. ⚠⚠ **TWO FIRST-ORDER PREDICTIONS OF THE MOVEMENT BOTH FAILED, IN OPPOSITE DIRECTIONS:** measured shrink against chr21's 2.94e-06 is **35.2x**, while `mean beta_v` scaling predicts 72x and sum(`beta_v`^2) scaling predicts 128.6x (9.290392e-06 chr21 against 7.223240e-08 genome-wide). So the movement is NOT a function of budget dilution alone. Positive correlation of matches along the HMM copying path is the UNTESTED explanation and needs `Var_Q(u)` from PAIR marginals, which this script does not compute — do not quote a mechanism for the 35.2 until it is measured. ⚠ The table was first emitted with the ratio row mislabelled `utility_retained` (which is the SAMPLED `u_path`, not a ratio); the file carries a dated RELABEL block and a regenerated table from `--from-json`, values bit-identical. | LIVE |
| `logs/B_exact_kl.md` + `logs/B_exact_kl.json` | 2026-10-12 | **THE ARM 2 PRIVACY CEILING, EXACT.** Emitted by `scripts/exact_kl.py`. `acc_ceiling` = **0.501555** at `tau` = 0.99 against a `p_succ_bound` of 0.9950, and 0.500029 to 0.500150 across the operating range; chr21, T = 332,521, K = 88, D1_extremal, 10 targets, zero Monte-Carlo error, 929 s. ⚠ An UPPER bound on EVERY attacker, including unwritten ones — strictly stronger than `attacker_accuracy`, which constraint 7 makes a LOWER bound on distinguishability. ⚠⚠ Validated by a check I did not plan for: `kl_exact / eta_tau^2` is constant to **0.31%** across a 53-fold range of `eta_tau`, which is the second-order scaling theory requires and would break if the normalisers or marginals were wrong. Gate `acc_ceiling <= p_succ_bound` passes 70/70 cells. ⚠ Pinsker is an inequality: the true ceiling is closer to 0.5, so these flatter us and the caveat must be quoted with them. | LIVE |
| `logs/B_verify_cohorts.md` | 2026-10-10 | **THE 22-CHROMOSOME COHORT TABLE, and the gate on quoting any genome-wide number.** Emitted by `scripts/verify_cohorts.py`. 22 of 22 pass every check; `K_states` = 88 throughout. **22,307,009 records, 22,134,562 chain positions.** ⚠ TWO exclusions, both target-independent: 305 with `support` == 0 (0.0014%), and 172,142 in 20 untraversable blocks (0.772%) under decision K. The pre-decision-K chain was 22,306,704. ⚠⚠ **`mean beta_v` per position is 4.5178e-08 against chr21-alone's 3.27e-06 — 72x smaller** (4.4830e-08 before decision K), so every chr21-scope per-position utility figure on record shrinks at genome scope. ⚠ **[CORRECTION - 2026-10-12]** This read "shrinks by that factor", i.e. 72x. MEASURED, the `e_q_u` movement shrinks **35.2x**, not 72x — and a sum(`beta_v`^2) argument predicts 128.6x, so the error is not a refinement in one direction but a failure of both first-order scalings. `mean beta_v` is still 72x smaller; it just does not set the movement. See `logs/B_exact_utility.md`. Declared rate spans **434x** (chr21 7.81% to chr6 0.02%) and the top four are four of the five acrocentrics. | LIVE |
| `data/cohort/DEPRECATED_chr21_coordinate_cut/` | 2026-09-29, retired 2026-10-10 | **chr21's ENTIRE pre-decision-C tree, 31 files**, moved aside rather than deleted because many recorded numbers were measured on it. Built under the coordinate cut (`12,968,320 < POS <= 46,680,243`, 305,886 of 340,849 positions) and carries axis digest `251b554f...`. ⚠ NOT comparable with the current build, which keeps 340,848 positions — 34,963 more — and whose three cM-source classes reconcile with this tree's dropped count EXACTLY (26,044 below the map + 592 above it on the declared rate, plus 8,327 interpolated plateau positions). | SUPERSEDED |
| `logs/B_build_all_cohorts_final.log` + `logs/B_build_all_cohorts.run1_failed.log` | 2026-10-10 | **The driver's own completeness verdicts, both kept.** The final run exits **0** with "all 22 chromosome(s) built and complete". ⚠ The `run1_failed` log is RETAINED deliberately: it is the record of nine chromosomes lost to a commit that landed mid-run, and it shows the driver refusing to report a partial genome (exit 1, "9 of 9 batch-3 chromosomes FAILED"). Do not delete it — it is the evidence that the fan-out asserts its own completeness. | LIVE |
| `scripts/verify_cohorts.py` + `tests/test_verify_cohorts.py` | 2026-10-10 | **CORRECTNESS of the built cohort tree, as distinct from its COMPLETENESS.** `build_all_cohorts.py` asserts that every chromosome ran and every expected file exists and is non-empty; a row-misaligned `support.npy` passes all of that and then silently scores the wrong variants, because the arrays are positional and a shift produces numbers rather than an error. Eight checks, each corresponding to a way this pipeline can actually go wrong: per-site array lengths agree; `support` equals BOTH the called-genotype count and `INFO/AN` (two independent routes to one quantity, so they CAN disagree); positions strictly increasing; no NaN and no decreasing cM inside the chain; `in_chain` exactly `support > 0` (decision C retired the coordinate cut, so any other exclusion is an accidentally reintroduced span); declared-rate positions only at the ENDS; and the manifest's `axis_consistent` and `chain_rule`. ⚠ All ten injected corruptions are CAUGHT and the unmutated control is clean; deleting any one check fails exactly its own test. 13 tests, which SKIP with a stated reason when no chromosome is built rather than reporting green on checks that never ran. | LIVE |
| `scripts/build_all_cohorts.py` + `src/runner.py` | 2026-10-10 | **The 1..22 driver and the parallel runner.** Fan-out is bounded by MEASURED memory, not cores: 22.3 MB of RSS per MB of gzipped GFA, from chr21's observed 9.48 GB peak, so chromosomes are bin-packed largest-first into batches that fit the budget (5 wide for the big ones, 9-10 for the small). `runner.run_jobs` RAISES on any failed job rather than returning a short list — a parallel map that drops a failure yields a result set that looks complete. 8 tests. | LIVE |
| `data/raw/ucsc_gap/gap.txt.gz` | 2026-10-10 | **The hg38 assembly gap table**, 12,368 B, sha256 `c466135b...`, from UCSC goldenPath. 829 rows, types scaffold/contig/telomere/heterochromatin/short_arm. ⚠ Acquired to answer the N-gap question for `declared_rate`: **79.8% of the map-silent sequence is gap** (68.8 of 86.2 Mb), so accumulating the declared rate across it would assert 109.5 cM of which **87.4 cM sits on placeholder bp**. The five `short_arm` records ARE the acrocentric p arms. ⚠ **No variant of ours falls inside any gap** — chr21: 0 of 340,849 — so gaps are pure placeholders whose only effect is inflating bp distance between variants that straddle them. | LIVE |
| `logs/B_exact_endpoint.md` + `logs/B_exact_endpoint.json` | 2026-10-10 | **The primary endpoint, exactly, at three candidate genome-wide shares.** Emitted by `scripts/exact_endpoint.py`. Two independent routes: marginals (linear functional, no differencing) and finite-difference derivatives of `log Z_p`. ⚠⚠ **The derivative route returns a NEGATIVE `Var_Q(u)` at chr21's genetic-length share**, so every recorded `Var_Q(u)` is a share-1 figure and must not be carried to genome scope. The marginals route is monotone with 3.45e+09x float64 headroom at the smallest share, and the endpoint's movement scales exactly linearly in the share. | LIVE |
| `logs/B_map_coverage.md` + `logs/B_map_coverage.json` | 2026-10-10 | **Per-autosome genetic-map coverage, the input to decision C's declared rate.** Emitted by `scripts/map_coverage.py`. Three classes reported SEPARATELY and never added: mapped, map-silent (no cM at all — what the declared rate covers, 3.00% of autosomal bp, 96% of it the five acrocentric p arms), and interior zero-cM (inside the span, 3.27%, which decision C does NOT override). Autosomal mean mapped rate **1.2707 cM/Mb**; per-chromosome spread 1.84x. ⚠ Position counts are chr21-only — the other 21 need their own VCFs. ⚠ The N-gap half is BLOCKED: no reference FASTA, AGP or assembly report is on disk. | LIVE |
| `scripts/map_coverage.py` | 2026-10-10 | The measurement above. Asserts each map's monotonicity at read time and RAISES if chr21's map-silent + in-span split stops equalling the recorded 34,371 — the gate that resolved two irreconcilable totals (ours 26,044 + 8,327, PanMixer's 26,018 + 8,327). | LIVE |
| `logs/A_churn_across_tau.md` | 2026-10-08 | **Churn of the mismatch set across `tau`**, both models, fixed seeds, on the COMMON 304,664-position axis both can emit. Emitted by `scripts/churn_across_tau.py`. Decomposes the net change into positions fixed and broken, so the question "does the tilt move a few or reshuffle many" is answered rather than inferred. ⚠ Supersedes the model-B shares in `A_models_B_vs_C_chr21.md`, which used the full chain's denominator against a 304,664-position numerator — see `docs/bugs.md`. | LIVE |
| `logs/A_models_B_vs_C_chr21.md` | 2026-10-07 | **Model C against model B, exact.** Emitted by `scripts/compare_models_chr21.py`. Model B: 59,472 units, `expected_switches` 3,394.1 against the 3,532 predicted independently on 2026-09-28, floor 14.5699% against C's 14.8077%, slope −1.7738 against −1.2381 cells per unit `eta_tau`. ⚠ Supersedes a first run whose transitions were start-to-start rather than inter-unit gaps — see `docs/bugs.md`. | LIVE |
| `logs/A_fidelity_frontier_chr21.md` + `logs/pmsweep.log` | 2026-10-07 | **The `target_fidelity` frontier: ours across `tau`, PanMixer's across `capacity`.** Emitted by `scripts/fidelity_frontier.py`, which fails loudly if the sweep log is incomplete. PanMixer swept at 8 capacities by `scripts/run_panmixer_capacity_sweep.sh` (sentinel 0): **1.20% to 16.41%, saturating from capacity 0.3**. Ours spans 14.8077% to 14.8070% across `tau` = 0 to 0.99 — a point, not a curve. Matched fidelity at **capacity ~= 0.17**. ⚠ Different targets and unequal denominators: read the SHARES, as a frontier. | LIVE |
| `logs/A_exact_fidelity_HG00096.md` | 2026-10-03 | **Exact expected released-vs-target mismatch count per `tau`, zero Monte-Carlo error.** Emitted by `scripts/exact_fidelity.py` from the posterior state marginals (`ffbs.state_marginals`), validated against the repo's own O(K^T) enumerator on the toy to **4.8e-15**. Carries BOTH the correct half-`eta_tau`-per-chain column and the full-`eta_tau` column a buggy sampling script used, so the 2x error stays visible. Headline: baseline 75,177.40 of 507,692 known cells (14.8077%), and the tilt reduces it by **1.2381 cells per unit `eta_tau`**, exactly linearly. **Supersedes every sampling estimate of this quantity.** | LIVE |
| `logs/A_results_toy_vs_chr21.md` | 2026-10-01 | **The full metric set, toy and chr21 as the two COLUMNS.** Emitted by `scripts/results_table.py`, never hand-assembled. Every row in `docs/terms.md`'s vocabulary, gate rows marked `*`, and `n/a` carrying its reason rather than being omitted. ⚠ Two deliberate departures, both stated in the file: the `tau` axis is folded into ROWS (floor and ceiling become rows) because the column axis is the configuration; and the two columns are NOT at matched ploidy — the toy gate is haploid, every chr21 release is diploid, and rows are marked `[hap]` or `[dip]` accordingly. That asymmetry IS the open gap. | LIVE |
| `logs/A_HG00096_chr21_v2.md` + `.releases.json` | 2026-09-30 | **The first RanPanMixer mechanism result, re-run after the RNG fix.** Same script, same seeds, same everything except `SeedSequence.spawn` streams per haplotype. `utility_retained` 0.6931, `donor_switches` 13,908, `sd(u)` at `tau` = 0.5 = 1.173e-03 against the 1.422e-03 the exact value predicts for a diploid draw. **Quote this one.** | LIVE |
| `logs/A_HG00096_chr21_v1.md` + `.releases.json` | 2026-09-30 | **The first RanPanMixer mechanism result.** Metrics-as-rows table over `tau` = 0 / 0.05 / 0.10 / 0.25 / 0.50 / 0.90, HG00096 diploid, chr21, `indicator/w=inv_support/beta=chain`, 3 draws per `tau`, seeds 1000-1002. Emitted by `scripts/run_mechanism_chr21.py`, never hand-assembled. The `.releases.json` carries one record per draw under the section 4 schema. ⚠ **A v1 DIAGNOSTIC RUN, NOT A PRE-REGISTERED RESULT** — `A-LCK-preregister` has not run and deferred decisions R1-R5 are unsettled. ⚠ **SUPERSEDED THE SAME DAY — DO NOT QUOTE.** Its diploid seeds were `seed` and `seed + 1`, so each release's haplotype 1 was bit-identical to the next release's haplotype 0: the three draws per `tau` were four correlated haplotype draws, and the resulting `sd(u)` was inflated about 2.2x. Retained so the deviation is visible and the correction in `docs/memory.md` 2026-09-30 is checkable. | SUPERSEDED |
| `src/mechanism/toy.py` (the toy cohort) | 2026-09-30 | The falsification fixture: K = 4, T = 7, hand-written, reads nothing from disk. 16,384 state paths -> 2,304 with mass -> 96 output paths in `Y_D`. Registered as an artifact because `A-THY-toy-enumeration` is a result that depends on its exact structure, and `toy_cohort()` asserts each feature at construction so a degenerate edit cannot pass silently. | LIVE |
| `primer/primer.html` + `primer/src/` | 2026-09-20 | Teaching primer for onboarding: a what-changed preface plus the landscape, the problem, PanMixer dissected, and our mechanism. ~57k words, built from `primer/src/` by `primer/build.py`. v2 expanded the hidden-Markov-model subsection (1.2.2) with full derivations of the forward algorithm and backward sampling, worked arithmetic, and a brute-force cross-check. URL: https://claude.ai/artifact/HzUWmUmxAPwXUSq3oFy8sb | LIVE |

## 8. Reproducibility contract

A result is reproducible here only if all six are recorded: the **arm**, the graph+cohort
id, `tau`, the utility id, the number of draws, and the seed. Two results are comparable
only when the first four match — and **never across arms**, which differ in K (86 vs 88),
in what is held out, and in whether the target is in the cohort at all.

⚠ **REVISED 2026-10-08 (decision E).** The held-out individual (`loo_id`) is **mandatory
in Arm 1 and forbidden in Arm 2**, not "comparison-arm only". Arm 2 holds out nobody, so
requiring one there would make a condition no external-target release can satisfy into a
gate on every result; Arm 1 holds out both of the target's haplotypes by definition, so a
missing `loo_id` there means the panel is unknown and the release is not reproducible.
⚠ Because Arm 1's leave-one-out is applied by **masking two haplotype rows in memory**
rather than by building a per-target cohort tree, the artifacts on disk are **identical
across all 44 targets and carry a single `axis_digest`** — the digest therefore cannot
witness which rows were dropped. `loo_masked_haplotype_ids` in the §4 record is the only
provenance for that, so a release record missing it is unreproducible even though every
artifact digest matches. ⚠ `tests/test_mechanism.py` currently hard-asserts
`rec["loo_id"] is None` with the reason "a genuine external-target release holds out
nobody"; that assertion is correct for Arm 2 and will fail on **every** Arm 1 release.
Fixing it is code work, outside this file.

**The private input is the target path.** The `PLANNED` target directory must
never be committed, copied into an artifact directory, or quoted in any document.
Nothing in it is covered by the release guarantee; only the sampled output path
is. Written as a rule now, before the directory exists, because the first run is
where it would be broken.
