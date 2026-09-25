# plan.md — the board

**Last updated:** 2026-09-25

Read at session start. This file exists so a new session never has to grep
`docs/memory.md` to know where things stand.

---

## Current state

- **Branch:** `main`. Remote `origin` -> `https://github.com/dylan33smith/RanPanMixer`.
- **Phase:** `A` — foundation. Goal: an end-to-end sampler provably correct on a
  graph small enough to enumerate, plus a fair head-to-head against PanMixer.
- **Still no mechanism code**, but as of 2026-09-25 there IS a pinned v1
  specification — see `A-IMP-v1-sampler` under In Progress. The upstream dependency
  is in, pinned, running, and has been taken apart in enough detail to know exactly
  which parts we inherit and which we replace. No RanPanMixer number exists yet.
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
  we do inherit is its DATA, its BLOCK structure, its forward recursion, and its
  EVALUATION suite. See the 2026-09-18 entry in `docs/memory.md` and the upstream
  section of `docs/bugs.md`; both findings were reproduced by hand.
- **The binding constraint right now** is the GRCh37/GRCh38 build mismatch in
  PanMixer's released pipeline. Until preprocessing is rebuilt on the GRCh38 30x
  panel, nothing derived from their variant-match set is usable by either arm —
  not the LD blocks, not the attack database, not `af_loss`, not `ld_loss`.
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
**graph+cohort id (including which individual was held out) · `tau` · utility id**.

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

### [A-IMP-v1-sampler] v1 SPECIFICATION — pinned 2026-09-25

**Scope decision (Dylan).** Maximise reuse of PanMixer and change exactly ONE thing:
the sampler. Keep its blocks, its data model, its donor panel and its whole
evaluation suite. Defer the chromosome-wide chain.

**INHERIT UNCHANGED**
- plink LD block boundaries computed on the 1000 Genomes panel, and the
  site-to-block assignment by position. Blocks stay exactly as PanMixer builds them.
- The data model: the int16 (subjects, sites, 2) matrix with `-1` for missing, and
  the block dictionary.
- The donor panel: 2(N-1) = 86 haplotypes, the target's own two removed.
- The evaluation suite, unchanged, by emitting `new_haplotypes`.

**REPLACE**
1. **Drop the anchor concept entirely.** Every variant in a block is a chain
   position. No dependency on the PanGenie callset. Rationale and evidence: the
   2026-09-25 entry in `docs/memory.md`.
2. **One code path for every block.** A single-variant block is a chain of length
   T = 1 — a tilted categorical draw over donors — not a special case. **There is no
   allele-frequency branch in v1.** This keeps all 340,824 chr21 variants
   cohort-supported instead of 73.9%, and removes the paper's "only one top-level
   SNP, not enough information" problem, which is a limitation of INFERRING a chain
   rather than copying donors.
3. **Exact tilted `ffbs`** in place of forward-only ancestral simulation. The forward
   pass must STORE every alpha column; PanMixer keeps only the running one.
4. **Re-derive the transition constants from the published Methods** (Ne = 10,000,
   r = 1.26, distance in centiMorgans, d = distance * Ne * r) — NOT from PanMixer's
   code, whose constant is the reciprocal of Ne and whose distance is a VCF row
   index. ⚠ Needs a GRCh38 genetic map, **not yet acquired**. Interim fallback: base
   pairs with a constant cM/Mb rate, stated explicitly wherever a number is reported.
5. **No knapsack.** Every block is resampled from the tilted distribution; nothing is
   released verbatim by default.

**WHY PER-BLOCK SCOPE IS SAFE.** The product of per-block tilted samplers equals the
global tilted distribution: taking the product over blocks of
`baseline_model` times exp(`eta_tau` * beta_b * phi_b) / Z_b yields
`release_distribution` exactly, provided `baseline_model` factorises over blocks,
the `beta_t` sum to 1, and each `phi_t` lies in [0,1]. Theorem 1 needs nothing more.
**So per-block now and genome-wide later is a scope choice, not a weaker guarantee.**

**CONSTRAINT THIS IMPOSES.** `phi_t` must decompose along the chain positions inside
a block. A whole-block similarity that does not decompose cannot be sampled exactly
by `ffbs`.

**OPEN DECISION — the missing-data policy.** What the emission does when a donor
carries `-1` at a position. Dylan is leaning **renormalise** (drop that donor from
the state distribution at that position) but has NOT decided; see `missing_policy`
in `docs/terms.md` for the three options and their biases. This matters because the
positions v1 adds are far worse behaved than anchors: mean missingness **11.31%**
versus **0.35%**, with 7,555 positions over 25% missing.
**Prerequisite for deciding it well:** `-1` conflates five distinct conditions and
they do not want the same treatment (see `missing_policy`). They ARE distinguishable
at conversion time — LV/PS, run-length, the CONFLICT tag, GT arity — but
`VCFtoNP` keeps only position and genotype, so all five are identical by the time
the mechanism sees them. **v1 preprocessing should emit an auxiliary REASON array
beside the allele matrix**, which keeps a cause-aware policy available without
committing to one now. Cheap now, unrecoverable later.

**EXIT GATES**
1. `A-THY-toy-enumeration` passes: sampler matches exact enumeration; exact TV <= `tau`
   across a grid of `tau` and many input pairs; `log_z_p` within [0, `eta_tau`] on
   every draw; `tau` = 0 reproduces `baseline_model` exactly.
2. One code path: singleton blocks provably traverse the same sampler with T = 1, and
   no allele-frequency branch exists in the source.
3. Transition constants give measured switch mass of order 0.1 to 1 at realistic
   separations, and a sampled block is demonstrably NOT a single donor copied verbatim
   (PanMixer fails this 300/300).
4. The missing-data policy is chosen, implemented, and covered by a test that fails
   if the behaviour changes.
5. Emits `new_haplotypes` so the inherited evaluation suite runs unmodified.

**MEASURED CONTEXT (chr21).** 11,670 multi-variant blocks plus 89,087 singletons =
100,757 chains covering all 340,824 variants; K = 88 states, 86 after removing the
target; `ffbs` cost about 30M operations. Compute is not a constraint.

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
  not documented; `beta_t` summing to 1 asserted; support computed on the
  leave-one-out cohort; `-1` masked before comparing; `diploid_utility` using MEAN
  with a test that catches a SUM.

### [A-IMP-ffbs-sampler] Algorithm 1
- **Why:** the mechanism. No backward pass exists anywhere upstream.
- **Prerequisite:** `A-IMP-cohort-hmm`, `A-IMP-utility`.
- **Exit gate:** forward alphas STORED for all t (PanMixer keeps only the running
  vector); tilt potentials multiplied in; per-position normalization; `log_z_p`
  range-asserted; O(TK) `stay_switch` path agreeing with an O(TK^2) dense
  reference; no argmax anywhere. Also emits `new_haplotypes` (see `docs/terms.md`).

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
4. **How to certify `output_support` independence.** Currently enforced by
   discipline; it should be enforced by construction. On PanMixer's data, target
   inclusion is the default, which makes this urgent rather than theoretical.
5. **Composition across releases.** One release per genome is an operational rule,
   not a theorem. The head-to-head's frontier is k correlated releases and is an
   evaluation device only.
6. **Can `Y_D` be made graph-valid?** PanMixer never checks that a stitched block
   path is a valid traversal, and the paper's hierarchical nested sampling is not
   implemented. If our support must contain only graph-valid paths, we build that.
7. **TRACKED (raised by Dylan 2026-09-22): do hypervariable structural sites make
   or break us?** On chr21 there are 7 sites where all 88 haplotypes carry distinct
   alleles, and 618 with >=20 alleles. They are simultaneously the most identifying
   sites on the chromosome and the ones PanMixer scores at exactly zero (population
   frequency 0 -> `-log(0)` -> zeroed). Our `tau` bound does not degrade with
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
8. **OPEN — what is a chain position? (`T`) — raised by Dylan 2026-09-22, to
   discuss with his PI; deliberately NOT decided here.**
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
