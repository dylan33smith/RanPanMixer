# plan.md — the board

**Last updated:** 2026-09-28

Read at session start. This file exists so a new session never has to grep
`docs/memory.md` to know where things stand.

---

## Current state

- **Branches:** `main` (docs, the PanMixer investigation). `pipeline` — an ORPHAN
  branch started 2026-09-28 holding the build, with `scripts/ configs/ logs/` and a
  gitignored `data/` symlink to `/data/ds85/RanPanMixer/data`. Remote `origin` ->
  `https://github.com/dylan33smith/RanPanMixer`.
- **Phase:** `A` — foundation. Goal: an end-to-end sampler provably correct on a
  graph small enough to enumerate. A head-to-head against PanMixer is now an
  OPTION for the writeup, not a requirement of the design — see below.
- **Still no mechanism code.** As of 2026-09-28 the architecture is settled
  (model C), the preprocessing input and output sets are specified, and the
  documentation has been corrected where the PanMixer comparison had been forced
  onto a different problem. No RanPanMixer number exists yet.
- **⚠ THE SCOPE CORRECTION OF 2026-09-28.** The v1 scope decision — *"maximise reuse
  of PanMixer, change exactly ONE thing"* — imported PanMixer's threat model, in
  which the target is a cohort member, into a design whose target is EXTERNAL to
  `G`, `D`, every panel and the attack database. An audit of the six docs plus the
  primer filed 91 findings; adversarial verification upheld 76 and found 51 more.
  Dylan's words: *"Boundaries aren't a weak function of the target's genotypes
  because the target is completely outside of any other dataset."* Corrected in
  commit `1f8ea3f` (`CLAUDE.md`, `terms.md`, `plan.md`, `data.md`). **When in doubt,
  the question to ask is: is this sentence about PanMixer, or about us?**
- **How heavily to lean on the comparison is now an open question for the PI.** The
  claim this project makes — a bound on TV between the releases of ANY two input
  paths — is a PROOF, and PanMixer's design structurally cannot make it: it releases
  unselected blocks verbatim, so two inputs differing in an unselected block have
  disjoint supports and TV = 1. A head-to-head on their metrics may not be needed at
  all. What IS needed is the toy enumeration and a utility curve. See
  `two_target_indistinguishability` in `docs/terms.md` — nothing in PanMixer's suite
  tests it.
- **What the investigation established** (all in `docs/memory.md`, all measured on
  chr21): PanMixer's privacy score reads only anchor variants, so in 32% of blocks
  the rarity it ignores exceeds the score it reports; 30.9% of its selected moves
  change nothing while booking 7.0% of the reported privacy credit; it rewrites
  3.13 Mb of sequence that its per-record utility accounting cannot see; `CONFLICT`
  records become missing genotypes and are therefore released verbatim; and `-1` in
  the matrix conflates four different conditions. None of this says the published
  numbers are wrong — it says the metrics measure something narrower than their
  names suggest, which matters because we were going to inherit them.
- **PanMixer is pinned** at `c182c38` with its 16 GB of inputs downloaded and
  checksummed, and a working conda environment. The cohort is **44 individuals /
  88 haplotypes** — so `K_states` ~ 88, not the ~1000 the proposal's runtime
  example assumes. Sampling will be cheap; the cost is all in preprocessing.
- **The headline from the ingest:** *"reuse PanMixer's cohort HMM as `baseline_model`"
  is not available.* Its Li-Stephens sampler, as shipped, never recombines — it
  reduces to a uniform draw over donor haplotypes copied verbatim per block. What
  we inherit is its DATA (the GRCh38-rebuilt allele matrix and data model) and its
  forward recursion. Its BLOCK structure is **not** inherited — under model C the
  chain steps variant to variant and LD blocks play no role in the mechanism — and
  its cohort-level utility metrics (`af_loss`, `ld_loss`) do not apply to a release
  that never edits the cohort. See the 2026-09-18 entry in `docs/memory.md` and the
  upstream section of `docs/bugs.md`; both findings were reproduced by hand.
- **The binding constraint right now** is the GRCh37/GRCh38 build mismatch in
  PanMixer's released pipeline. Nothing derived from their variant-match set is
  usable by the **PanMixer arm or the head-to-head** — not the LD blocks, not the
  attack database, not `af_loss`, not `ld_loss`. ⚠ **Our own sampler depends on
  none of them.** Under model C it consumes only `G`, the cohort allele matrix, the
  GRCh38 genetic map, and the target path; there is no external panel in the
  mechanism's input set at all.
- **The one thing most likely to go wrong** is unchanged and now has a concrete
  shape: an implementation that quietly lets `output_support` or `baseline_model`
  depend on the target. On PanMixer's data that is the DEFAULT, not an accident —
  its allele frequencies, support counts and scoring panel all include the target,
  because in its threat model the target is a cohort member.

## Headline result

No RanPanMixer mechanism result exists. The table below is the ingest measurement
that gates use of PanMixer's data. Metrics are rows; panels are columns.

| metric | 1000G Phase 3 (wired by the pipeline) | 1000G 30x (repo's unwired path) |
|---|---|---|
| reference build | GRCh37 (`assembly=b37`, hs37d5) | GRCh38 (contig chr21, length 46,709,983) |
| samples | 2,504 | 3,202 |
| chr21 records | 1,105,538 | 1,002,752 |
| exact (POS,REF,ALT) matches vs pangenome chr21 * | 1,025 | 258,610 |
| relaxed (POS,REF) matches | 2,176 | 268,851 |
| match rate (exact / 340,824 pangenome records) * | 0.30% | 75.88% |

**Provenance.** HPRC v1.0 PGGB GRCh38 VCF (sha256 `ead25541...`), chr21 extracted
by index, 340,824 unique (POS,REF,ALT) records · 1000G Phase 3 chr21 (sha256
`1942e070...`) · 1000G 30x chr21 (sha256 `a925c112...`) · `bcftools 1.24` ·
measured 2026-09-18. The 30x column is the CONTROL: same pangenome, same
comparison, a same-build panel.

**COLUMNS**
- *Phase 3* — what `get_1000g_phased.sh` downloads and what `get_blocks.sbatch`
  and `get_mappings.py` consume. GRCh37.
- *30x* — what the repo's unwired `get_blocks_grch38.sh` downloads. GRCh38, and
  what `diploid_gap_score.py` prefers when present.

**ROWS**
- *reference build* — the whole finding in one row. The pangenome is GRCh38.
- *samples* — also the attack-database size; a denser panel is a strictly
  stronger attacker, so this row moves the privacy axis.
- *exact matches* — the set that becomes `pangenome_mask.npy`, which drives the
  attack DB, the `af_loss` MAF strata and the entire `ld_loss` site set.
- *match rate* — 0.30% is consistent with coincidental collisions across builds,
  not with two panels describing the same coordinates.

**SYNTHESIS.** The released pipeline bins GRCh38 pangenome positions into LD-block
intervals computed on GRCh37 coordinates and matches variants across the two
builds. The control settles it: the same pangenome against a same-build panel
matches **252x more records** (75.88% vs 0.30%), so the deficit is the build, not
the data. At 0.30% the derived artifacts carry almost no real signal,
and the affected metrics fail quietly — `ld_loss` returns `(0.0, 0)` and every
`af_loss` stratum goes 0/0 rather than raising. The repo ships an unwired GRCh38
path and code branches that prefer 3,202-sample artifacts no released script
builds, which suggests the authors' own runs used a path the release does not
wire in. **This does not show the paper's published numbers are wrong** — we have
not reproduced their figures and cannot say which code produced them. It shows the
released pipeline as documented does not reproduce the paper, and that we must
rebuild preprocessing on the GRCh38 panel before either arm runs.

## Validation against the published paper (2026-09-18)

Metrics are rows; the two arms are columns. Published = Nature Communications
DOI 10.1038/s41467-026-77591-0.

| metric | published reports | we measured | verdict |
|---|---|---|---|
| cohort size | 44 individuals | 44 / 88 haplotypes | match |
| capacity -> normalized utility loss | budget = fraction of max loss | 928.92/9288.85 = 10.000% at capacity 0.1 | exact |
| per-subject cost, longest chromosome | 5 min / 14 GB | 68 s / 1.7 GB (chr21); x5.47 -> 6.2 min / 9.3 GB | consistent |
| determinism, allele-bound invariants | not reported | byte-identical per seed; 0 violations | pass |
| chr21 variant / block / anchor counts | **not reported anywhere** | 340,824 / 100,757 / 273,475 | not comparable |
| read mapping, BASELINE (chr21, 5 external donors) | 77.83 perfect / 95.62 gapless / 77.01 MAPQ60 | 79.05 / 95.96 / 79.23 | **mismatch** (+1.2 / +0.3 / +2.2) |
| `eps_private` and the remaining downstream metrics | 0.002; AF 0.006/0.006/0.005; LD 0.003; obfuscated mapping 77.82/95.61/77.01 | — | not yet run |

**Provenance.** `paper/s41467-026-77591-0_reference.pdf` · `external/PanMixer` @
`c182c38` · chr21 rebuilt on the 30x GRCh38 panel · HG00438 · capacities 0.1 / 0.5
· seeds 123 / 456.

**COLUMNS** · *published* — the version of record; the preprint's numbers differ and
must not be cited. · *we measured* — chr21 only, one subject, on a substituted panel.

**ROWS** · *cohort size* is the only headline entity count the paper states that we
can check, and it matches. · *capacity* confirms the knapsack normalization is
exactly right. · *per-subject cost* is a scope-corrected extrapolation, not a
measurement at their scope. · *counts* have no published counterpart at all — the
paper reports no variant, block or anchor numbers. · *downstream metrics* are the
real reproduction targets and none has been run.

**SYNTHESIS.** Every *structural* check matches, and those matches were not
guaranteed to — particularly the exact capacity normalization and the saturation
behaviour. But "correct installation" here means *we run their code as shipped*,
not *their code implements the published method*: four confirmed divergences
between code and Methods stand. And the one published RESULT we have attempted —
the read-mapping baseline — came out in the right regime but **not** within
tolerance, with the residual offset unexplained after three hypotheses were tested.
No reported result has been reproduced.

## Reporting contract

Two results are comparable only if their config stamps match on:
**graph+cohort id · `tau` · utility id** — plus `loo_id` for comparison-arm runs only
(our target is external; a genuine release holds out nobody).

Metrics are rows, configurations are columns. Every table carries a provenance
line and is followed by one bullet per column, one bullet per row, and a prose
synthesis. Every number in prose traces to a table cell. Emit tables from the
scoring script; never hand-assemble one.

The floor is `tau` = 0 (an untilted draw from `baseline_model`). The ceiling is the
target's own path, `u_path` = 1. Both are columns, not prose asides.

⚠ **`eps_pmi` and `tau` are incomparable, and so are `capacity` and `tau`** (see
`docs/terms.md`). The head-to-head aligns the arms on MEASURED empirical attack
success only, and reports both frontiers.

## Ledger — Phase A

| ID | Intervention | Endpoint | n | Result | Verdict | memory |
|---|---|---|---|---|---|---|
| A-FIX-panmixer-framing | Audit the docs for PanMixer framing forced onto our different problem | every instance found and scoped | 6 files, 12 agents | 91 filed, 76 upheld, 51 more found by verifiers; CLAUDE.md constraint 1 was re-seeding it every session | done | 2026-09-28 |
| A-DAT-genetic-map-parm | Can any map give chr21 genetic distance below 10.3 Mb? | a defensible rule for every variant | chr21, 10 agents | no usable map; `np.interp` clamping + the map's own centromeric plateau give 34,345 zero-distance positions = 7.25 Mb emitted verbatim; cut at 12,968,320 | done | 2026-09-28 |
| A-DAT-input-set | What does the mechanism actually need as input? | a minimal, verified input set | full checkout | model C needs only graph + genetic map + target; PanGenie and the external panel drop out; `INFO/AN` == `support(v)` bit-identical | done | 2026-09-28 |
| A-DAT-evaluation-fit | Which PanMixer evaluators transfer to our release object? | a per-evaluator verdict | full suite | only the gap score, and only after 3 changes; `af_loss`/`ld_loss` undefined for us; nothing in their suite tests our actual claim | done | 2026-09-28 |
| A-FIX-docs-system | Adopt the six-file documentation system and write the verifier | verifier passes on a greenfield tree | n/a | 6 files + 20 checks | done | 2026-08-29 |
| A-DAT-panmixer-ingest | Pin PanMixer, acquire its data, map it against the proposal | a verified reuse map + a runnable upstream | 22 agents, 0 errors | reuse map done; 2 blockers found; 4 proposal claims corrected | done | 2026-09-18 |
| A-DAT-grch38-repin | Rebuild chr21 preprocessing on the GRCh38 30x panel; run PanMixer | get_mappings finds ~258k matches, not ~1k | 1 chr, 4 runs | 258,610 strict / 268,851 relaxed; PanMixer runs in 68 s, deterministic, 0 invariant violations | done | 2026-09-18 |
| A-DAT-missingness | Establish what `-1` actually means before choosing a policy | a taxonomy with shares | chr21 | 4 distinct causes; 66.1% of nested missingness is "parent allele does not contain this bubble"; 349 runs >=100 are assembly gaps | done | 2026-09-25 |
| A-DAT-anchors | Why anchors exist and whether we can drop them | the rationale, tested | chr21 | paper's reason is top-level SNPs (supported: nested are 22.5% uncalled vs 2.25%); the code's PanGenie proxy is neither necessary nor sufficient | done | 2026-09-25 |
| A-DAT-noop-moves | Do selected moves actually change anything? | share of no-ops | 47,070 moves | 30.9% change nothing, carrying 7.0% of reported privacy credit at zero utility cost | done | 2026-09-24 |
| A-DAT-af-table | Where the allele frequencies come from | the branch structure | chr21 | 3-way priority; 2 corrections issued; our repin puts the target in its own frequencies | done | 2026-09-24 |
| A-DAT-conflicts | Characterise CONFLICT records and the block-spanning output semantics | what they do downstream | chr21 | conflicts become missing data -> zero privacy value AND 1.6x cost -> released verbatim; output is per-record allele indices with no coherence check | done | 2026-09-23 |
| A-DAT-hypervariable-sites | Characterise the worst-case multi-allelic structural sites | a measured profile of the class | chr21 | 7 sites where all 88 haplotypes differ; 2,485 matrix rows for one 257 kb bubble; these score eps_pmi = 0 | done | 2026-09-22 |
| A-EVL-readmapping | Reproduce the published chr21 read-mapping BASELINE | within ~0.1 pt of 77.83 / 95.62 / 77.01 | 5 donors, ~11M reads each | 79.05 / 95.96 / 79.23 — right regime, NOT an exact match; offset unexplained | partial | 2026-09-18 |
| A-EVL-published-validation | Check our run against the PUBLISHED paper (Nat Commun) | our numbers vs theirs | 9 agents | installation CORRECT; 3 published numbers changed vs preprint; build inconsistency is in the paper | done | 2026-09-18 |
| A-IMP-local-slurm-shim | A blocking local sbatch so the pinned checkout runs unmodified | array expansion + failure propagation | 4-task self-test | shim drives the whole pipeline; failures exit non-zero | done | 2026-09-18 |

**Provenance.** No mechanism runs. The ingest row's evidence is workflow
`wf_a30688a6-3d7` plus hand verification recorded in `docs/memory.md`.

## In Progress

### [A-IMP-preprocess] THE PREPROCESSING PIPELINE — specified 2026-09-28

**What this is.** Everything up to the point the sampler can run. Written fresh, not
ported: we copy functions from PanMixer where they are right, but the file structure and
the pipeline are ours, so that nothing is inherited without a decision. Under model C this
is a much smaller job than it was under the v1 spec.

**TWO STAGES, TWO DIRECTORY TREES, ENFORCED BY THE SIGNATURE.** Standing constraint 1 says
the baseline is built from `G` and `D` only. Rather than enforce that by discipline — which
`docs/plan.md` open question 4 has flagged as inadequate since the start — the **cohort
stage takes no target argument at all.** Target leakage stops being something to remember
not to do and becomes something the interface cannot express.

#### Inputs — the whole set

| # | input | what it is | state |
|---|---|---|---|
| 1 | **The pangenome graph** | HPRC v1.0 PGGB, published as GFA (no `.gbz`/`.snarls` for freeze1). chr21 = 434,772,428 B, sha256 `85222d4f...`. `vg deconstruct` is **pipeline step 1**, not an assumed input (Dylan, 2026-09-28) | `OK` — staged 2026-09-29 |
| 2 | **A genetic map** | GRCh38 cM table. A DOWNLOADED artifact with cM already computed — not something estimated from our data. Columns `pos / chr / cM`, 44,618 rows for chr21 | `OK` — `external/genetic_maps/chr21.b38.gmap` |
| 3 | **The target** | a phased GRCh38 VCF, handled by the SEPARATE target stage below. Test targets: held-out 1000G 30x samples, pool of 3,085 | `OK` (pool selected) |

**What model C removed from this list.** No LD blocks means no PLINK, so **no external
panel in the mechanism's input set**. No anchors means **no PanGenie callset** (measured
2026-09-28: exactly two live consumers, both deleted by our design; no evaluator reads it;
saves a 5.06 GB download). No allele-frequency branch means **no `allele_frequencies.npy`
and no `get_af.py`** with its three-dataset entanglement. The external panel survives only
as the attack database in the evaluation, where cohort overlap must be removed before it is
loaded.

⚠ **We do our own deconstruct, and it will NOT match PanMixer's VCF.** Measured on chrY,
`vg` 1.68 recovers 99.55% of the published VCF's sites by POS but only **92.7% by
(POS, REF, ALT)**, concentrated at multi-allelic sites, because the published file was built
with `vg` 1.36. Accepted knowingly (Dylan, 2026-09-28). Two consequences: we must apply
chrX removal and `chm13` removal ourselves (the published VCF carries 45 samples, we need
44); and **if the head-to-head is ever run, both arms must run on the SAME VCF.**

#### Outputs — cohort stage (`preprocess-cohort --chrom 21`, no target argument)

| file | shape / dtype | what it is for |
|---|---|---|
| `haplotypes.npy` | `(88, n_sites)` int16 | the donor panel. Flattened from `(44, n, 2)` because the HMM's states are HAPLOTYPES, not people |
| `reason.npy` | `(88, n_sites)` int8 | why each `-1` is missing. Codes 0-4; see `reason_array` |
| `missing_runs.tsv` | `(sample, strand, start_row, end_row)` | every run of consecutive `-1`. Assembly gaps are recorded as INTERVALS, not per-cell codes, so no length threshold is baked into the artifact — see `reason_array` |
| `positions.npy` | `(n_sites,)` int64 | bp position per row |
| `genetic_pos.npy` | `(n_sites,)` float64 | cM per row, interpolated with `left=nan, right=nan` |
| `support.npy` | `(n_sites,)` int32 | non-missing haplotype count. ⚠ Read `INFO/AN` — measured bit-identical, no computation needed |
| `allele_lengths` | ragged int32 | length of every declared allele. 2.9 MB against 118.3 MB for the sequences; needed for length-weighted utility and rewrite accounting |
| `sites.tsv` | `n_sites` rows | row -> `(chrom, pos, ref, alt, n_alt, LV, PS)`. The join back to the VCF, the only place nesting survives, and the row-order check for any artifact handed to another tool |
| `haplotype_ids.npy` | `(88,)` | which donor is which `(sample, strand)` |
| `manifest.json` | — | input URLs + checksums, tool versions, parameters, `chain_span`, and a digest of `sites.tsv` |

#### Outputs — target stage (`preprocess-target --target <id> --chrom 21`)

| file | shape / dtype | what it is for |
|---|---|---|
| `path.npy` | `(n_sites, 2)` int16 | the target's allele at each row. **The private input.** |
| `representability.json` | — | count and total `-log f` of target variants with NO record in the graph. ⚠ AUDIT ONLY — adding a record for a target-specific variant would make `output_support` target-dependent and void Theorem 1 |

#### Layout

```
data/
  raw/            chr21.hprc-v1.0-pggb.gfa.gz , chr21.b38.gmap
  cohort/chr21/   haplotypes.npy reason.npy positions.npy genetic_pos.npy
                  support.npy allele_lengths.npy sites.tsv haplotype_ids.npy manifest.json
  targets/<id>/chr21/   path.npy representability.json manifest.json
```

**Per chromosome, not all in one go.** The chain cannot cross chromosomes — they are
independent — preprocessing is embarrassingly parallel across them, and a driver can loop
1..22. chr21 alone is one invocation today.

⚠ **chr21 IS A TESTING SCOPE. THE PRODUCT IS GENOME-WIDE** (Dylan, 2026-09-29). Working on
one small chromosome is a way to iterate in minutes instead of days; it is not the
deliverable. Most of the extension really is just the driver loop — but three things are
NOT, and they should be visible now rather than discovered later:

1. **`beta_v` must sum to 1 over whatever scope `tau` covers.** Today it sums to 1 over
   chr21, which silently means "`tau` per chromosome". A genome-wide `tau` needs a
   genome-wide normaliser, so each chromosome gets a SHARE of the budget rather than all
   of it. Running 22 chromosomes each normalised to 1 would be 22 independent releases of
   the same person and would compose — exactly what standing constraint 4 forbids. This is
   the same scope algebra already recorded for per-block vs genome-wide: the product of
   the per-scope tilted samplers equals the global one provided the weights sum to 1
   ACROSS the whole scope. **Decide before any multi-chromosome run.**
2. **`chain_span` is per chromosome and must be derived, not copied.** Each chromosome has
   its own map span and its own centromere, and chr13, 14, 15, 21 and 22 are all
   acrocentric, so four more will have the same unassembled p-arm problem. The rule
   ("start where the map carries real markers") generalises; the numbers do not.
3. **The graph is 15.64 GB compressed / 86.5 GB raw genome-wide.** Per-chromosome GFAs
   total 15.89 GB and are the better form, because deconstruct is per-chromosome anyway
   and a single-shot whole-genome run would need several hundred GB of RAM against our
   251 GB.

Compute is otherwise not a constraint: alphas are per chromosome, so memory does not grow
with the genome, and chr21's ~30M operations scale linearly.

**Target as a separate command** (Dylan, 2026-09-28) so a new target does not re-run the
cohort stage. ✅ **ANSWERED 2026-09-29 — the target arrives as a phased GRCh38 VCF**, and
the mapping is a join on `(POS, REF, ALT)` against `sites.tsv`. No alignment step, and no
`vg giraffe` — which matters, because HPRC marks short-read mapping **"untested"** for the
PGGB graph. ⚠ If a real target later arrives as reads or an assembly, that is a materially
harder step and this decision reopens.

**The test target, until a real external genome exists** (Dylan, 2026-09-29): a held-out
sample from the 1000G 30x panel, which is phased and GRCh38-called, so a sample column IS
a path once joined to our site list. Eligibility is **not** simply "panel minus cohort":
the 39 HPRC donors in the panel are trio CHILDREN and their **78 parents sit in the
unrelated-2504 set**, so a random draw has a **~3.1% chance of being a cohort member's
parent**, sharing about half its genome with a donor haplotype. Pool = 3,202 − 39 − 78 =
**3,085**, listed at `pipeline:configs/eligible_targets.chr21.txt`; the rule and the
pedigree check are in `docs/data.md`. ⚠ A panel sample UNDERSTATES representability —
9.9% of its calls have no record in the graph and that is a LOWER bound, since its
variants are already catalogued by construction.

#### What we deliberately do NOT produce

`allele_frequencies.npy`, the anchor mapping, `biallelic_snp_mask.npy`, the PanGenie matrix,
and — under model C — `blocks_dict.json`, `simple_blocks_idx.npy` and the PLINK outputs.
Keep the block build as a SEPARATE optional script so the model-B fallback stays one command
away, but it is off the main path.

#### Exit gates

1. Runs end to end from the graph + the genetic map with no PanMixer artifact on the input
   side, and no target argument reachable by the cohort stage.
2. Every chain position inside `chain_span` carries a genuine interpolated cM; asserts no
   NaN survives.
3. `support.npy` equals `INFO/AN` (a free correctness check that the matrix and the VCF are
   row-aligned), and `num_alleles` bounds every observed allele index.
4. `reason.npy` carries codes 0-4 under the pinned precedence (conflict > not applicable
   > inherited > uncategorised), with a test that fails if the precedence changes;
   `missing_runs.tsv` reproduces the measured chr21 totals (1,056,726 cells in runs
   >= 100, and the >= 10 / >= 1000 endpoints at 92.06% / 68.67%), proving no threshold
   was baked in.
5. `sites.tsv` digest recorded in `manifest.json`, and re-verified by anything that consumes
   a downstream artifact.

---

### [A-IMP-v1-sampler] v1 SAMPLER SPECIFICATION — pinned 2026-09-25, NARROWED 2026-09-28

**Scope decision (Dylan, 2026-09-25; NARROWED 2026-09-28).** The original wording was
*"maximise reuse of PanMixer and change exactly ONE thing, the sampler — keep its blocks,
its data model, its donor panel and its whole evaluation suite."* That went too far: it
imported PanMixer's threat model, in which the target is a cohort member, into a design
whose target is external to `G`, `D`, every panel and the attack database. **Narrowed:**
we reuse PanMixer's DATA and its forward recursion. Blocks as a chain unit, the
donor-panel-minus-target, and the cohort-distortion metrics are COMPARISON-ARM choices,
not properties of our design. **The chromosome-wide chain is no longer deferred — model C
is the architecture (decided 2026-09-28), with model B retained as a fallback if the PI
wants closer parity with PanMixer.**

**INHERIT UNCHANGED**
- The data model: the int16 (subjects, sites, 2) matrix with `-1` for missing,
  row-aligned to the chromosome's VCF record order.
- The forward recursion's log-space stay/switch form, with constants re-derived from
  the published Methods (see REPLACE item 4).

**COMPARISON-ARM ONLY — not properties of the design**
- plink LD block boundaries and the site-to-block assignment. Under **model C** the
  chain unit is the individual variant and blocks play no role in the mechanism. They
  are retained only for the model-B fallback and for parity runs against PanMixer.
  ⚠ Do NOT compute them from the 44-sample cohort: measured 2026-09-28, that gives
  7,858 blocks covering 50.4% of variants against 14,137 / 74.2% from the panel, with
  the largest block 2,562 variants against 618. Only 44% of SNPs clear MAF 0.05 in 88
  haplotypes. If blocks are wanted, they need a large external panel.
- The donor panel is **K = 88, the full cohort.** Our target is external to `D`, so no
  haplotype is removed. The head-to-head arm instead uses 2(N-1) = 86, because there a
  cohort member stands in as target and externality has to be simulated. Removing both
  of that stand-in's haplotypes (never just one) is what PanMixer's departure (c) gets
  wrong — PanMixer-side toy: leaving the target in scores 1.0000, removing one 0.8571,
  removing both the honest 0.5714.
- PanMixer's ATTACK evaluators and its read-mapping utility, by **optionally** emitting
  `new_haplotypes`. ⚠ `af_loss` and `ld_loss` are NOT inherited: they are
  cohort-distortion metrics, undefined for an external-target release that never edits
  the cohort. Our utility axes are `target_fidelity` (always against its ceiling) and
  `utility_retained`. The claim our guarantee actually makes — two-target
  indistinguishability — is tested by NOTHING in their suite and is ours to write.

**REPLACE**
1. **Drop the anchor concept entirely.** Every variant in a block is a chain
   position. No dependency on the PanGenie callset. Rationale and evidence: the
   2026-09-25 entry in `docs/memory.md`.
2. **One code path for every block.** A single-variant block is a chain of length
   T = 1 — a tilted categorical draw over donors — not a special case. **There is no
   allele-frequency branch in v1.** This keeps all 340,824 chr21 variants
   cohort-supported instead of the 71.4% sampled by donor copying today (73.9% sit in multi-variant blocks, but those with <= 1 anchor are frequency-drawn too; `docs/memory.md` 2026-09-26), and removes the paper's "only one top-level
   SNP, not enough information" problem, which is a limitation of INFERRING a chain
   rather than copying donors.
3. **Exact tilted `ffbs`** in place of forward-only ancestral simulation. The forward
   pass must STORE every alpha column; PanMixer keeps only the running one.
4. **Re-derive the transition constants from the published Methods** (Ne = 10,000,
   r = 1.26, distance in centiMorgans, d = distance * Ne * r) — NOT from PanMixer's
   code, whose constant is the reciprocal of Ne and whose distance is a VCF row
   index. ✅ **GRCh38 genetic map ACQUIRED 2026-09-27** — `external/genetic_maps/chr21.b38.gmap`,
   verified byte-identical between the SHAPEIT4 and Beagle distributions. The interim
   base-pair fallback is no longer needed. ⚠ 7.63% of chr21 variants fall below the map
   span and get flat cM, so 296 multi-variant blocks cannot recombine at all.
5. **No knapsack.** Every block is resampled from the tilted distribution; nothing is
   released verbatim by default.

**WHY PER-BLOCK SCOPE IS SAFE.** The product of per-block tilted samplers equals the
global tilted distribution: taking the product over blocks of
`baseline_model` times exp(`eta_tau` * beta_b * phi_b) / Z_b yields
`release_distribution` exactly, provided `baseline_model` factorises over blocks,
the `beta_t` sum to 1, and each `phi_t` lies in [0,1]. Theorem 1 needs nothing more.
**So per-block now and genome-wide later is a scope choice, not a weaker guarantee.**

**⚠ MODEL C DOES NOT REMOVE THE DOUBLE-COUNT (raised 2026-09-28, open).** A parent record
and its nested children are separate ROWS of the allele matrix, so under C they are separate
CHAIN POSITIONS — the hypervariable chr21 region at 14,569,980 contributes **2,485 positions**
(one top-level record plus 2,464 nested children) for what is one stretch of DNA, and each
gets its own `beta_v`. So that locus draws 2,485/340,824 of the budget. C dissolves the
plink-versus-snarl segmentation question and the block-boundary incoherence, but NOT this:
the double-count is a property of the VCF's representation, not of the segmentation.
Options not yet assessed: collapse a parent and its descendants into one position, weight by
`LV`, or length-weight `beta_v`. **Not decided; measure the inflation first.**

**CONSTRAINT THIS IMPOSES.** `phi_t` must decompose along the chain positions inside
a block. A whole-block similarity that does not decompose cannot be sampled exactly
by `ffbs`.

**DECIDED 2026-09-28 — the missing-data policy is RENORMALISE.** A donor carrying `-1`
at a position is dropped from the state distribution there; see `missing_policy` in
`docs/terms.md` for the rejected options and their biases. Dylan accepted the known cost
deliberately: RENORMALISE is CORRECT for "not applicable" (the largest class, 66.1% of
nested missingness — the DNA genuinely is not on that chromosome) and does penalise
poorly-assembled donors for a defect that is not theirs. Released pangenomes should be
well assembled; if they are not, that is something to note in the writeup rather than
model around at this stage.
⚠ It acts on the STATE SPACE, not the transition — our baseline knows nothing about
alleles, so an excluded donor simply cannot be occupied at that position.
⚠ TO MEASURE once the sampler runs: leaving a donor and returning costs TWO switch
events, so `ffbs` will prefer to switch away and STAY away — missingness acts as a switch
TRIGGER, not a blip, and it is spatially clustered (28,052 runs, mean 20, max 23,478).
⚠ GUARD REQUIRED: 1 chr21 site of 340,824 has all 88 donors missing and 89 leave <= 1, so
the renormalised denominator can be zero. Declare the fallback, count it, report it. This matters because the
positions v1 adds are far worse behaved than anchors. Like for like, v1 steps over every variant (mean missingness **4.26%**, 30,800 over 25%) where PanMixer's HMM steps over anchors in >=2-anchor blocks (**0.26%**, 538). The non-anchor records inside multi-variant blocks are worst: **11.31%**, with 7,555
over 25% missing (against 0.35% over all anchors; scopes in `docs/memory.md` 2026-09-26).
**Still a prerequisite, for a cause-aware refinement later:** `-1` conflates **four**
distinct conditions on chr21 (see `missing_policy`; a fifth, haploid genotypes, was
expected and does NOT occur — all 237,594 non-piped GT fields have exactly one distinct
value, `.`). They ARE distinguishable at conversion time — LV/PS, run-length, the
CONFLICT tag — but
`VCFtoNP` keeps only position and genotype, so all five are identical by the time
the mechanism sees them. **v1 preprocessing should emit an auxiliary REASON array
beside the allele matrix**, which keeps a cause-aware policy available without
committing to one now. Cheap now, unrecoverable later.

**SECOND PREREQUISITE — the representability audit.** `Map(g, G)` must emit, beside the
path, a count of target variants it could NOT represent. Measured 2026-09-27: about **9.9%**
of a real external target's chr21 non-reference calls have no record in `G` — and they carry
**20.5%** of its `-log f` information — so they never
become chain positions, never enter `phi_t`, and are invisible to AF loss, LD loss and read
mapping alike. Silently dropped and silently unmeasured. **This count is an audit artifact
only.** Adding a record for a target-specific variant would make `output_support` depend on
the target and void Theorem 1 by the §4.4 disjoint-support argument — the fix is to measure
the loss, never to represent it. See `target_fidelity` for the reporting rule.

**DECIDED 2026-09-28 — THE CHAIN ARCHITECTURE IS MODEL C.** Raised 2026-09-27 from first
principles rather than from PanMixer's implementation; settled by Dylan 2026-09-28.
**Model B is retained as a fallback** if the PI wants closer parity with PanMixer; model A
is ruled out. The text that preceded this decision described **model A**, which was
inherited as a DESCRIPTION OF PANMIXER'S CODE and never argued for on its own terms.
Measurements are chr21, the acquired GRCh38 map, published constants (Ne = 10,000,
r = 1.26, n = 86). Analytic segment length for reference: n/(Ne*r) = 0.0068 cM, ~**4.5 kb**.

| model | chain | E[donor switches] | mean segment | alphas | status |
|---|---|---|---|---|---|
| **A** per-block chains, no inter-block transitions | 100,757 separate | 104,379 | 0.39 kb | 69 MB | RULED OUT |
| **B** one chain, units = **blocks**, blocks atomic | 100,757 | 3,532 | 11.60 kb | 69 MB | fallback |
| **C** one chain, units = **variants** | 340,824 | 7,155 | 5.73 kb | 234 MB | **DECIDED** |

**What choosing C removes from the project**, and the reason it was chosen: with the
individual variant as the chain unit there are no LD blocks in the mechanism, hence no
PLINK run and **no external panel in the mechanism's input set at all**; no anchor
concept, hence **no PanGenie callset** (measured 2026-09-28: it has exactly two live
consumers, both of which our design deletes, and no evaluator reads it); and no
allele-frequency table. C also deletes the singleton special case — there are no T = 1
chains, no per-block `Z_b`, one forward pass and one backward sample. **It is simpler than
what was pinned, not harder.** The external panel survives only as the attack database in
the evaluation, where cohort overlap must be removed before it is loaded.

- **A is the worst of the three and in the wrong direction** — 11x too much recombination.
  Its fresh-uniform start per block is algebraically the `d` -> infinity limit of the
  Li-Stephens switch formula, i.e. it asserts consecutive blocks are infinitely far apart
  genetically. Measured truth: median inter-block switch mass **0.0022** against the 0.9884
  that independence implies, with **98.5%** of block pairs below half of it. Consecutive
  blocks are as tightly linked as adjacent variants within one (median gap 1.52e-05 cM
  versus 1.81e-05 cM).
- **B discards 51% of the recombination** — the within-block half — but block sizes are
  median 1 and p90 3, so the loss concentrates in the tail. The largest block is 618
  variants over 39.5 kb, which the correct model recombines through 73% of the time and
  which B emits as one donor verbatim. That is a mild form of departure (b): the release
  carries another individual's real haplotype intact over tens of kb.
- **C is the correct Li-Stephens model** and lands nearest the analytic segment length.
  Costs 165 MB more and about 3x the operations, on a budget of 29M — not a constraint.
- **Blocks are a poor chain unit on this data regardless**: 89,087 of the 100,757 entries
  (88.4%) are orphan variants that fell in NO called plink block, so a "block boundary" is
  mostly a bookkeeping outcome rather than a recombination site.
- **"Stay close to PanMixer" does not favour any of them.** All three reuse the same data
  model and donor panel; PanMixer contains no backward pass at all, so the chain is our
  code under every option. ⚠ They do NOT all reuse the plink blocks — that was the
  original wording and it is wrong for C, which has no blocks — nor the evaluation suite,
  whose utility half does not apply to any of them.

**What does NOT change between them:** the privacy budget. `eta_tau` is identical and
`beta` still sums to 1 over the whole path, so the §4.10.4 thin-budget problem is untouched
by this choice. Demonstrated on the toy: the per-position tilt and the utility gain both
fall as 1/T regardless of what a position is.

**If C is chosen**, the spec's per-block factorisation argument is no longer needed — it
justified decomposing into independent per-block samplers, and C has one chain, so Theorem 1
applies to the whole path directly. C also DELETES the singleton special case: there are no
T = 1 chains, no per-block `Z_b`, one forward pass and one backward sample. It is simpler
than what is pinned, not harder.

**Numerical requirements.** Normalise every column (`alpha_t` = `r_t`/`c_t`,
`log_z_p` = sum of log `c_t`); an unnormalised forward pass underflows float64 after about
**183 positions** and we need up to 340,824. Use **float64** for the stored alphas — under
model C this is load-bearing, not a preference. float32 resolution is 1.19e-07 while the
per-position tilt at tau = 0.5 is **1.61e-06 under C**, only **13.5x** the noise floor,
against 5.45e-06 and **45.7x** under a block partition. Rounding would consume a meaningful
fraction of the target's entire influence. Alpha storage is **240 MB** under C and **71 MB**
under the model-B fallback, both at K = 88; a K = 86 head-to-head figure is ~2.3% smaller
and must say so. ⚠ Break the chain at genuine discontinuities — outside `chain_span`, the
centromere, long assembly gaps — never at block boundaries.

**EXIT GATES**
1. `A-THY-toy-enumeration` passes: sampler matches exact enumeration; exact TV <= `tau`
   across a grid of `tau` and many input pairs; `log_z_p` within [0, `eta_tau`] on
   every draw; `tau` = 0 reproduces `baseline_model` exactly.
2. One code path, and under model C there is no special case to test: every chain
   position is a variant, so no singleton/T = 1 branch exists at all. What the gate still
   tests is that **no allele-frequency branch exists in the source**. (The original
   wording was written against the block-structured spec, where singletons were the
   special case; retained here only for the model-B fallback.)
3. ⚠ **MIS-SPECIFIED — must be rewritten before it is used as a gate** (see
   `docs/memory.md` 2026-09-27, A-DAT-genetic-map). As written it demands switch mass of
   order 0.1 to 1 between adjacent chain positions, but those sit about 12 bp apart, where
   the biologically correct value is the measured 0.0026 — the gate would reject a correct
   model. Its second clause, "a sampled block is demonstrably NOT a single donor copied
   verbatim", fails for about 91% of blocks even with the correct map, because an LD block
   is by definition a stretch with little recombination and v1 chains within one block.
   Replacement candidates, all measured 2026-09-27: expected donor switches per chromosome
   (3,622.6, against PanMixer's ~0); median P(mosaic) among blocks with >= 20 variants
   (0.2747); or measured segment length against the analytic n/(Ne*r) = 0.0068 cM.
4. The missing-data policy is chosen, implemented, and covered by a test that fails
   if the behaviour changes.
5. Every chain position carries a genuine interpolated cM — no clamped values survive
   into the chain. Interpolate with `left=nan, right=nan` and assert no NaN remains
   after the coordinate cut, so the failure is loud rather than silent.

*(The former gate 5 — "emits `new_haplotypes` so the inherited evaluation suite runs
unmodified" — is not a gate of the mechanism. It is one option for one comparison and
now lives under `[A-EVL-headtohead]`.)*

**MEASURED CONTEXT (chr21).** Under model C the chain is one pass over all 340,824
variants, less the coordinate cut at both ends (305,887 positions). K = 88 states — our target
is external to the cohort, so no haplotype is removed; K = 86 applies only in the
head-to-head arm, where a cohort member stands in as target. `ffbs` cost is about 30M
operations at K = 88. Compute is not a constraint.

**COORDINATE CUT (decided 2026-09-28; closed at both ends 2026-09-29).** The chain runs
over **12,968,320 < POS <= 46,680,243** — exactly where the genetic map carries real
markers — keeping **305,887** of 340,824 variants (89.75%). ⚠ **PROVISIONAL: Dylan to
settle with his PI.** It is the cut that lets v1 run, not a defended answer.
The telomeric end was closed for the same reason as the p arm: 592 variants sat above the
map end forming a 593-position zero-distance run, and extrapolating there would invent
genetic distance in a subtelomeric region. Cost 0.19%. Closing both ends makes
`genetic_pos_valid` and explicit chain segments unnecessary — inside the span every
position has a real cM, and a long variant-free stretch produces a large `Delta_x` and
therefore near-uniform switching, which is already correct. Below that, `np.interp`'s clamping plus the
map's own zero-cM centromeric plateau give **34,345 consecutive positions with identical
cM**, so `P(switch) = 0` and the sampler emits one donor's real haplotype verbatim across
7.25 Mb — departure (b) at ~200x the scale of the worst block. Five independent reasons
support the cut: GRCh38 chr21 5.01-10.81 Mb is placeholder model sequence with 25
fabricated N-gaps; 1000G's strict accessibility mask rules 99.98% of those variants
inaccessible; our cohort matrix is 36.76% missing there against 6.57% elsewhere; deCODE's
pedigree map assigns 0.0 cM/Mb to every 1 Mb bin from 0-13 Mb; and a four-generation
pedigree observed zero p-arm allelic crossovers in 107 transmissions.
⚠ **REVISIT.** This is a pragmatic cut to get v1 running, not a settled answer. Maps that
nominally cover the region exist (pyrho hg38 from 5,088,754 bp; the Eagle redistribution
has 66 markers there) but spend ~10 cM across intervals defined by two markers, i.e.
interpolation across a void. Two smaller zero-distance defects also remain unaddressed:
**592 variants above the map end** (a 593-position zero-distance run at the telomere) and
**618 interior zero-cM map intervals** holding ~1,650 kept variants. The interior ones are
real map plateaus and are correct to leave; the telomeric run is not and needs a decision.

## Backlog

Ordered. Each item names what must pass before it starts and what closes it.

### [A-IMP-local-slurm-shim] A blocking `sbatch` shim — DONE 2026-09-18
- **Why:** PanMixer funnels every job through `sbatch`, directly in the data
  pipeline and via `external/PanMixer/tools/common/slurm_helper.py` in the CLI. One shim unblocks
  the whole toolkit and keeps the pinned checkout byte-identical.
- **Prerequisite:** none.
- **Exit gate:** parses `#SBATCH --array=A-B[%M]`, runs `bash <script>` once per
  index with `SLURM_ARRAY_TASK_ID` exported, writes the `out_`/`error_` files
  `check.py` greps for, BLOCKS until done (which makes `--dependency=afterok` a
  no-op), exits non-zero on any task failure, prints an id for `--parsable`, and
  caps concurrency by the declared `--mem` against 251 GB.

### [A-THY-toy-enumeration] Brute-force correctness gate on a toy graph
- **Why:** the ONLY test that can falsify our implementation of Theorem 1. At
  scale `tv_empirical` is biased upward and cannot falsify anything.
- **Prerequisite:** `A-IMP-cohort-hmm`, `A-IMP-utility`, `A-IMP-ffbs-sampler`.
  Deliberately independent of all PanMixer data.
- **Exit gate:** (a) sampler matches exact enumeration within Monte-Carlo error;
  (b) exact TV <= `tau` across a grid and many input-path pairs; (c) `log_z_p` in
  [0, `eta_tau`] on every draw; (d) `tau` = 0 reproduces `baseline_model` exactly.
- **Kill criterion:** any exact-TV value above `tau` on the toy support.

### [A-IMP-cohort-hmm] The target-independent prior — NEXT
- **Why:** `baseline_model` is what the guarantee rests on.
- **Prerequisite:** none for the toy version.
- **Exit gate:** rho and A built from the cohort alone, with a test that fails if
  any target-derived quantity reaches the constructor. **Transition constants
  re-derived from the paper (Ne = 10,000, r = 1.26, `Delta_x` in cM,
  d = `Delta_x`*Ne*r) — NOT from PanMixer's code**, which uses 1/Ne and a variant
  ROW INDEX as distance. Gate: measured switch mass is O(0.1-1) at realistic
  separations, and a sampled block is NOT a single donor copied verbatim.

### [A-IMP-utility] Bounded utility from PanMixer's support weights
- **Why:** an unbounded `u_path` removes the guarantee rather than weakening it.
- **Prerequisite:** none for the toy version.
- **Exit gate:** `phi_t(a,b) = 1 - (sum_v w_v 1[a_v != b_v]) / W_t` with
  `w_v = 1/support_D(v)`, `phi_t` = 1 when `W_t` = 0; asserted in [0,1] at runtime,
  not documented; `beta_t` summing to 1 asserted; support computed on the FULL cohort
  (our target is external, so there is nothing to hold out — leave-one-out support
  applies only to the head-to-head arm); `-1` masked before comparing;
  `diploid_utility` using MEAN with a test that catches a SUM.
  ⚠ Under model C this collapses: `phi_v` is an indicator (1 if the drawn donor carries
  the target's allele at `v`, else 0) and `beta_v = w_v / sum_u w_u`. Bounded in [0,1]
  and additive by construction, with no per-block normaliser.
  ⚠ `support(v)` needs no computation: measured 2026-09-28, the pangenome VCF's own
  `INFO/AN` field is bit-identical to `np.sum(pangenome != -1, axis=(0,2))` across all
  340,824 chr21 records.

### [A-IMP-ffbs-sampler] Algorithm 1
- **Why:** the mechanism. No backward pass exists anywhere upstream.
- **Prerequisite:** `A-IMP-cohort-hmm`, `A-IMP-utility`.
- **Exit gate:** forward alphas STORED for all t (PanMixer keeps only the running
  vector); tilt potentials multiplied in; per-position normalization; `log_z_p`
  range-asserted; O(TK) `stay_switch` path agreeing with an O(TK^2) dense
  reference; no argmax anywhere. May ALSO emit `new_haplotypes` for the comparison
  arm (see `docs/terms.md`), but that is optional and not part of this gate.
  ⚠ **Never materialise the transition matrix.** Li-Stephens `A` is
  `(1-s_t)*I + (s_t/K)*J` — rank-one plus diagonal — so the forward step is
  `alpha_t(j) = [(1-s_t)*alpha_{t-1}(j) + (s_t/K)*sum_i alpha_{t-1}(i)] * tilt_t(j)`,
  O(K) per position from one scalar sum. The backward draw is O(K) the same way.
  Store only the per-position switch scalar `s_t` (T floats) and the alphas
  (T x K). A materialised `A` per position would be 340,824 x 88 x 88 x 8 B = ~21 TB.

### [A-LCK-preregister] Freeze the evaluation before looking at it
- **Why:** the proposal requires the utility function be chosen before any attack
  outcome is inspected; otherwise every later privacy claim is unfalsifiable.
- **Prerequisite:** `A-THY-toy-enumeration` passing.
- **Exit gate:** a committed, dated file pinning the `tau` grid, the utility id,
  the primary endpoint, the attacker set, and **which attack database** is used.

### [A-EVL-headtohead] The comparison on PanMixer's data
- **Why:** the direct comparison this project is for.
- **Prerequisite:** `A-DAT-grch38-repin`, `A-LCK-preregister`.
- **Exit gate:** for each subject, a `loo_cohort` with donor panel, allele
  frequencies, support counts, site list and block assignment all recomputed
  without the target — applied to BOTH arms, with PanMixer's as-published numbers
  also reported so the cost of the correction is visible. Both arms emit
  `new_haplotypes`; the shared evaluators run unchanged. R >= 10 seeds per
  (subject, level) in both arms — PanMixer's published frontier rests on ONE
  presampled draw reused across capacities. Attacks run PER operating point, never
  pooled, because a frontier of k points is k correlated releases and our guarantee
  is one-release-per-genome.
- **Must include:** a `target_fidelity` axis, or our `tau` = 0 point looks like a
  free lunch — a pure-prior draw is a real cohort haplotype and scores near-perfectly
  on every cohort-level utility metric while carrying zero target information.

### [A-ATK-matched-pair] The discriminating experiment
- **Why:** the only measurement where the two guarantees differ BY CONSTRUCTION.
  PanMixer releases unselected blocks verbatim, so two inputs differing in an
  unselected block have disjoint output supports — TV = 1. Ours bounds every pair.
- **Prerequisite:** `A-EVL-headtohead` plumbing.
- **Exit gate:** a matched 1-vs-1 attacker over (target, decoy) pairs reporting
  WORST-CASE advantage, which our theorem bounds by `p_succ_bound` for every pair.
  Include the shipped `sampled_gap_score.py` (full-resample, no selection) as the
  PanMixer-side counterpart of our `tau` = 0 point.

## Blocked

| Item | Blocked on |
|---|---|
| all-autosome use of PanMixer's LD blocks and attack DB | only chr21 has been rebuilt on GRCh38; chr1-20,22 still need the 30x panel |
| reproducing the paper's privacy numbers | the 30x attack database: three `.npy` files no released script produces |
| read-mapping utility | five read FASTQs, `chr21.fa` and `hg38_cleaned.fa` that no script in the repo produces; also needs a Java runtime, absent from PanMixer's `environment.yaml` |
| all-autosome results | only chr21 of the 1000G panel is downloaded |
| our own data | not yet acquired; deferred deliberately |

## Dropped — with reasons

| Item | Why |
|---|---|
| Viterbi / MAP decoding of the tilted HMM | Not a draw from `release_distribution`; no guarantee. Retired in `docs/terms.md` before use. |
| Per-target pruning of the HMM state space | Makes `output_support` target-dependent and voids Theorem 1. |
| Reusing PanMixer's sampler as `baseline_model` | As shipped it never recombines; see `docs/bugs.md`. We take the transition FORM and the forward recursion, and re-derive the constants. |
| Reusing `eps_pmi` inside the mechanism | It is the self-information of the private input. Fatal to the guarantee; kept only as a plotting axis. |
| `hmm_1000g.py` as a starting point | Dead and broken upstream: nothing imports it, it references an undefined attribute, and its forward pass is an O(A*K^2) double loop. |

## Open questions

1. **ANSWERED 2026-09-18: yes.** The GRCh38 panel matches 75.88% of chr21
   pangenome records against 0.30% for the wired GRCh37 panel. The remaining ~24%
   is expected — the pangenome carries graph-only and structural variants absent
   from a SNV/INDEL/SV panel — but should be characterised before it is assumed benign.
2. **Which 1000G panel is the faithful one?** The published Methods and Data
   availability pin GRCh37 Phase 3, but the graph is GRCh38 and the official
   GRCh38 liftover of Phase 3 was withdrawn in 2021. We substituted the 30x GRCh38
   panel (what the repo's unwired script fetches and what `constants.py`'s 3202
   implies). Every number we produce is conditioned on that choice, so it must be
   stated in every result. Resolving it properly may need to ask the authors.
3. **What `tau` is defensible.** `tau` = 0.10 bounds equal-prior identification at
   0.55. A policy question we can inform but not settle; pre-register either way.
4. **ANSWERED 2026-09-28 — enforced by construction, not discipline.** The
   preprocessing splits into a COHORT stage and a TARGET stage writing to separate
   directory trees, and **the cohort stage takes no target argument at all**, so
   target leakage is not expressible through the interface. `chain_span` is likewise
   a fixed coordinate rule derived from the assembly and the map, frozen before the
   target file is opened. ⚠ Residual to hold: never let an exclusion become
   target-dependent (e.g. "drop positions where the target is missing") — that
   would leak, and it is the one shape of this error the directory split does not
   catch. Note also that on PanMixer's data target inclusion IS the default; that is
   correct for their threat model and is not a constraint we inherit.
5. **Composition across releases.** One release per genome is an operational rule,
   not a theorem. The head-to-head's frontier is k correlated releases and is an
   evaluation device only.
6. **Can `Y_D` be made graph-valid?** PanMixer never checks that a stitched block
   path is a valid traversal, and the paper's hierarchical nested sampling is not
   implemented. If our support must contain only graph-valid paths, we build that.
7. **TRACKED (raised by Dylan 2026-09-22): do hypervariable structural sites make
   or break us?** On chr21 there are 7 sites where all 88 haplotypes carry distinct
   alleles, and 618 with >=20 alleles. They are simultaneously the most identifying
   sites on the chromosome and the ones PanMixer scores at exactly zero — not via
   `-log(0)` zeroing, which never fires on our data (Correction 7), but because the
   privacy score reads only anchor alleles and these sites are not anchors
   (Correction 8). Our `tau` bound does not degrade with
   rarity, so this is where a path-level guarantee should show an advantage — but
   three things must be settled first: whether `output_support` can contain a valid
   traversal there at all, whether `phi_t` should be length-weighted when one block
   spans 257 kb, and how to avoid double-counting a top-level record against its
   2,464 nested children. **MEASURED 2026-09-22:** 9.2% of records sit inside
   another record's span, 15.4% of blocks contain such a record, and 23.2% of
   variants live in an affected block — so block disjointness fails broadly, not
   locally. This is a design decision, not an edge case to document. Leading option:
   derive `T` from the snarl tree rather than from plink intervals, which makes
   blocks disjoint by construction. See the 2026-09-22 entry in `docs/memory.md`.
8. **ANSWERED 2026-09-28 (Dylan): ONE VARIANT PER POSITION — model C.** Model B
   (one position per LD block) is the retained fallback if the PI wants closer
   parity with PanMixer; model A is ruled out. The reasoning that settled it is in
   the chain-architecture table above: A's fresh-uniform restart per block is
   algebraically the infinite-genetic-distance limit and produces 104,379 donor
   switches against ~3,532 for B and ~7,155 for C, where the analytic segment
   length says ~4.5 kb. C lands nearest, is the correct Li-Stephens model, and is
   SIMPLER than what was pinned — no singleton special case, no per-block `Z_b`,
   one forward pass, one backward sample. It also removes the external panel, the
   PanGenie callset and the allele-frequency table from the input set entirely.
   ⚠ Two sub-questions bundled here were also answered: **anchors are dropped**
   (no external callset in the mechanism), and the double-counting worry is moot
   under C, where every row is its own position. ⚠ STILL OPEN, carried forward: a
   path containing a tandem duplication visits a position TWICE, so "the allele of
   path p at position v" is not well-defined for such a path — which `u_path` and
   `phi_v` both assume it is. That is now a `nested_variant` question, not a
   segmentation question.
   *The original framing and the three candidates are preserved below.*
   **[SUPERSEDED — retained for the reasoning]**
   `T` is the number of positions in the hidden Markov chain: the number of steps
   the sampler takes from one end of the chromosome to the other. At each position
   it decides which cohort haplotype to copy from and emits whatever that donor
   carries there. So `T` is not a free parameter — it is fixed by what we declare a
   "position" to BE, and that choice determines the chain length, what "switch
   donors here" means, and whether positions are disjoint.
   The candidates, with what each costs:
   - **plink LD blocks (what PanMixer uses).** T ~ 100,757 on chr21. Boundaries come
     from an external panel in bp coordinates, know nothing about graph structure,
     and are NOT disjoint: 9.2% of records sit inside another record's span and
     15.4% of blocks contain such a record. Breaks the additive decomposition of
     `u_path` and can emit incoherent paths.
   - **Snarl-tree levels.** Blocks disjoint by construction; every boundary is a real
     graph cut point, so switching donors at a boundary is graph-valid. The 257 kb
     region becomes 1 position instead of 370 overlapping ones. Objection — snarls
     are not LD blocks — is weak for us, because in a copying model the linkage is
     carried by the donor panel, not by where boundaries fall. Cost: we stop
     inheriting PanMixer's segmentation.
   - **One variant per position.** Maximal resolution, largest T, worst overlap problem.
   ⚠ Also unresolved within this question: a path containing a tandem duplication
   visits a position TWICE, so "the allele of path p at position t" is not
   well-defined for such a path — which `u_path` and `phi_t` both assume it is.
   **Settle before writing the sampler; it defines what a position is.**
   ⚠ Bundled into this question (Dylan, 2026-09-25): whether we keep any notion of
   "anchor". PanMixer's anchors are a PanGenie overlap used as a proxy for
   "well-behaved chain position"; the paper's justification is about TOP-LEVEL SNPs
   and is supported by the data (nested records are 22.5% uncalled vs 2.25% for
   top-level SNPs), but the code's proxy admits 13,436 nested variants and drops
   18,468 qualifying top-level SNPs. Compute is not a reason (1.25x). Leading
   answer: derive positions from graph structure, not from an external callset.
9. **ANSWERED 2026-09-18 (chr21, GRCh38 repin).** Measured on the rebuilt chr21:
   14,137 plink LD blocks, but `blocks_dict.json` holds **100,757** blocks because
   **88.4% are singletons** — one variant, in no LD block. By BLOCK count only
   **10.4%** have the >= 2 anchors that put them on PanMixer's HMM path; by VARIANT
   mass those blocks hold **71.4%** of the genome. So the HMM governs most variants
   but a small minority of blocks. Consequences: our `T_blocks` is ~100,757 for
   chr21 alone with `K_states` = 88, so FFBS is ~9M operations — trivially cheap;
   and any claim about "PanMixer's cohort HMM" applies to 10.4% of blocks, with the
   other 89.6% drawn from independent per-variant allele frequencies.
