# plan.md — the board

**Last updated:** 2026-10-12 (later — genome-wide utility measured). Chunks 0-7 DONE. ⚠ Next: chunk 8 (Arm 1 leave-one-out releases). ⚠⚠ The open scientific question is UTILITY, not privacy: genome-wide `e_q_u` MOVES +8.345e-08 across the whole `tau` range, and the 0.8586 level is reached at `tau` = 0 so it is not ours to quote — see the Current state bullet.

Read at session start. This file exists so a new session never has to grep
`docs/memory.md` to know where things stand.

---

## Current state

- **Branches:** `main` (docs, the PanMixer investigation). `pipeline` — an ORPHAN
  branch started 2026-09-28 holding the build, with `scripts/ configs/ logs/` and a
  gitignored `data/` symlink to `/data/ds85/RanPanMixer/data`. Remote `origin` ->
  `https://github.com/dylan33smith/RanPanMixer`.
- **Phase:** `B` — **THE TWO-ARM PAPER** (opened 2026-10-08). Phase A (foundation:
  an end-to-end sampler provably correct on a graph small enough to enumerate)
  **CLOSED 2026-09-30, goal met.**
  - **Arm 1 — the PanMixer head-to-head.** The 44 HPRC cohort members as targets,
    diploid, 22 autosomes, per-target leave-one-out at K = 86, running PanMixer's
    own evaluations so the comparison is direct. **A REQUIRED paper section** — the
    earlier wording, "an OPTION for the writeup, not a requirement of the design",
    is superseded.
  - **Arm 2 — the external target.** The work of Phase A, carried to genome scale,
    with its own attacks and utility axes.
  ⚠ **Every statement about the target must name its arm.** See `CLAUDE.md`
  Project Core; conflating them is what the 2026-09-28 scope correction was about.
- **PREPROCESSING IS BUILT AND HAS RUN (2026-09-29).** `src/` holds our own code —
  the first in this project that is not PanMixer's. The cohort stage and all ten
  target stages completed on chr21; 17 unit tests and 28 consistency checks on the
  real output pass.
- ✅ **THE MECHANISM IS BUILT AND GATED (2026-09-30).** `src/mechanism/` holds
  `baseline_model` (`hmm.py`, no target parameter, asserted from the AST), the bounded
  utility, and an exact `ffbs` sampler. `A-THY-toy-enumeration` passes every clause on a
  hand-written K = 4, T = 7 cohort that reads nothing from disk: `ffbs` matches exact
  enumeration on both the output distribution and the state-path JOINT, exact TV never
  exceeds `tau` over all 3,571,128 input pairs per `tau`, `log_z_p` stays in
  `[0, eta_tau]`, and `tau` = 0 reproduces `baseline_model` to 1e-15. 46 new tests, **87
  in the suite, all passing.** The O(TK) rank-one forward pass agrees with an O(TK^2)
  dense reference to **1.11e-16**.
- ⚠⚠⚠ **THE OPEN SCIENTIFIC QUESTION IS UTILITY, AND IT IS NOT A DEFECT — IT IS THE
  RESULT (raised by Dylan 2026-10-12).** Privacy is settled and strong: `acc_ceiling` =
  **0.501555** at `tau` = 0.99 against a `p_succ_bound` of 0.9950, computed exactly. But
  the SAME thinness that makes the mechanism undetectable makes it useless:
  - `utility_retained` on chr21 at SHARE 1 is **0.6931 at `tau` = 0** and **0.6932 at
    `tau` = 0.9**. The mechanism moves the FOURTH DECIMAL. Exact expected-utility
    movement across the whole grid: **2.94e-06**.
  - ⚠ **THE 0.6931 IS NOT OURS.** It is obtained at `tau` = 0, where the release is a
    pure cohort mosaic carrying ZERO target information — any two humans share most of
    their alleles. Quoting the LEVEL as the mechanism's utility is exactly the free-lunch
    error `target_fidelity` exists to prevent. **The movement is the figure.**
  - ✅ **MEASURED GENOME-WIDE 2026-10-12, EXACTLY** (`scripts/exact_utility.py`, 44
    marginals passes, no sampling, no finite difference; HG00096, ARM 2, 22 autosomes,
    shares by decision J summing to 1.000000 over T = 22,134,562, scope "chain", n = 1):
    `e_q_u` = **0.708967646** at `tau` = 0 and **0.708967729** at `tau` = 0.99,
    `u_ceiling` = 0.825708741, `e_q_u_over_ceiling` = 0.858617102 and 0.858617203.
    **THE MOVEMENT IS +8.345e-08** (+0.000012%). Provenance `logs/B_exact_utility.md`.
  - ⚠ **AND THE 0.8586 IS NOT OURS EITHER** — same free-lunch error as the 0.6931, same
    reason: it is reached at `tau` = 0. Never write "we retain 86% of utility".
  - ⚠⚠ **TWO PREDICTIONS OF THAT MOVEMENT BOTH FAILED, IN OPPOSITE DIRECTIONS.** This
    bullet used to read: `mean beta_v` falls 3.27e-06 -> **4.52e-08**, a factor of 72, so
    the movement "should shrink by about that again". Measured shrink is **35.2x**. The
    better first-order argument — movement ~ `eta_tau` * Var_Q(u), and Var_Q(u) ~
    sum(`beta_v`^2) if matches were independent, measured 9.290392e-06 on chr21 against
    7.223240e-08 genome-wide — predicts **128.6x**. The two bracket the measurement and
    neither predicts it, so **the movement is not a function of budget dilution alone**.
    Untested explanation: matches are POSITIVELY CORRELATED along the HMM copying path, so
    Var_Q(u) exceeds sum(`beta_v`^2). Settling it needs `Var_Q(u)` from PAIR marginals,
    which `scripts/exact_utility.py` does NOT compute. **Quote no mechanism for the 35.2
    until that is measured.** `mean beta_v` is still 72x smaller; it just does not set the
    movement.
  - ✅ **THE MECHANISM IS NOT BROKEN AND THE TOY PROVES IT**: at `mean beta_v` = 0.111
    instead of 4.52e-08, `E_Q[u]` moves **0.4502 -> 0.6477** across the grid. It is BUDGET
    DILUTION, which is a far better finding to write up, and it is what R5 was reopened to
    quantify.
  - ⚠ **NOT YET MEASURED AT ALL:** Arm 1 utility at any scope (needs chunk 8), and
    PanMixer's own four metrics on our build (chunk 9).
  ⚠ A reviewer will put the privacy ceiling and the utility movement side by side
  immediately. They are the same fact from two sides, and the paper has to say so rather
  than be caught by it.

- ✅✅ **CHUNK 6 IS COMPLETE: ALL 22 AUTOSOMES AND THE FULL PANEL, BOTH GATED
  (2026-10-10).** The cohort is **22,307,009 records / 22,134,562 chain positions**;
  driver exit 0, verifier 22 of 22. **TWO exclusions, both target-independent:** 305
  positions with `support` == 0 (0.0014%), and **172,142 positions in 20 untraversable
  blocks** (0.772%) under decision K. ⚠ **[CORRECTION - 2026-10-10]** this read
  "22,306,704 chain positions, with only 305 positions excluded genome-wide" before
  decision K, which was the figure the first genome-wide build produced and which the
  forward pass could not actually release from. The panel is **21 files / 29.3 GB / 69,689,262 records / 3,202 samples**
  on every chromosome, plus chr21's pre-existing file, so all 22 are covered; verifier
  21 of 21 with the **reference build checked against our own `deconstruct` output** rather
  than against typed-in constants.
  ⚠⚠ **THE HEADLINE IS THE DILUTION: `mean beta_v` per position is 4.5178e-08 genome-wide
  against chr21-alone's 3.27e-06 — 72x smaller.** (4.4830e-08 before decision K, which
  removed 0.772% of the chain and so raised the per-position budget by 0.778%.) This is `B-THY-toy-genome`'s budget-dilution
  effect measured on the real axis instead of projected, and **every chr21-scope per-position
  utility figure in this file shrinks by that factor at genome scope.** It sharpens R5 rather
  than resolving it: the whole-chromosome endpoint was already 1:900 against a single draw's
  noise at chr21 scope.
  ⚠ A publishable side finding: the declared rate spans **434x** across autosomes (chr21
  7.81% to chr6 0.02%) and the top four are **four of the five acrocentrics**, whose p arms
  the genetic map does not cover. It quantifies how much of each chromosome's transition
  structure rests on a declared constant rather than a measured map.
  ⚠ Two failures got us here, both mine and both archived in `docs/memory.md` section 9:
  a commit that landed MID-RUN cost nine chromosomes and split the damage by process start
  time, and the manifest was found to describe the CONFIG rather than the BUILD — it
  relabelled chr21's 11-day-old coordinate-cut arrays as `chain_rule` = "declared_rate_fill".
  Both are now gated in code.
- ⚠⚠ **NEXT: THE UTILITY'S SCOPE. R5 is a live, quantified design choice, not a dead end**
  (reopened 2026-09-30). `Var_R(u)` is set by how many mosaic SEGMENTS the utility spans, so
  a CONTIGUOUS window of 1,000 positions gives **94.9x** the variance of 1,000 scattered
  ones at identical `sum beta`. And the ceiling is exact:
  `gain/noise per draw = eta_tau * sqrt(Var) <= eta_tau/2`. Our configuration needs
  **6,552,960** (the DIPLOID figure; 3,276,480 is the haploid one — corrected 2026-10-01) draws to show the tilt at 2 sigma; a contiguous window of 100 needs about
  **1,600**; the best possible at `tau` = 0.5 is **53**. **The question for Dylan and his PI
  is what the release is FOR** — whole-chromosome fidelity, which is unmeasurable at any
  feasible n, or regional fidelity, which is measurable this week. Pre-registration item
  under constraint 6, so it must be settled BEFORE `A-LCK-preregister` and before any
  evaluation. See R5.
- ⚠ **DO NOT RUN THE EVALUATION FIRST.** At the current utility it would pass trivially and
  prove nothing: the release carries almost no target information, so `attacker_accuracy`
  would sit at chance and two-target TV near 0 — which standing constraint 7 says is not
  evidence of anything. And `utility_retained` is **85.0% of `u_ceiling` at `tau` = 0 too**
  (0.6931 against 0.6932 at `tau` = 0.9), so quoting it as the mechanism's utility would be
  exactly the free lunch `target_fidelity` warns about — now measured on `u_path` ITSELF,
  not merely on the cohort-level metrics that warning was written about.
- ⚠⚠ **THE FIRST RanPanMixer NUMBER: THE TILT RESHUFFLES THE PATH AND BUYS NO
  UTILITY.** HG00096, chr21, `beta` scope `chain`. Between `tau` = 0 and `tau` = 0.9 the
  released path MOVES — 8 of 10 disjoint-stream draws differ, by a mean of **417 state
  positions** and **54 allele positions** per haplotype — while expected utility moves
  **2.94e-06** across the whole `tau` grid. `utility_retained` is **0.6931** against a
  `u_ceiling` of 0.8155 at every `tau` (0.6932 at `tau` = 0.9). ⚠ **[CORRECTION -
  2026-10-10]** This bullet read **0.6940**, the superseded v1 value from the
  correlated-seed run, which this same file already retracts two sections down
  (0.6940 -> 0.6931, at `tau` = 0 through 0.50, n = 10 disjoint-stream draws, over
  253,846 called chain positions) and which line 78 above already quotes correctly. A retraction that does not reach the board's present tense has not
  landed; `docs/plan.md` line 525 still carries it inside a historical entry, where
  it belongs. Not a bug: `log_z_p` tracks `eta_tau` at a
  constant 69.5% and the sampler is proven exact on the toy. Measured sampling-free from
  the derivatives of `log Z_p`: `Var_Q(u)` = **4.046e-06**, `sd(u)` = **2.011e-03**, so
  the gain is about **1:900** against a single draw's noise.
  ⚠ **AND REWEIGHTING `beta_v` CANNOT FIX IT — measured 2026-09-30.** Concentrating the
  budget 3,000-fold raises `sum beta^2` 2,200x but `Var_Q(u)` only **16.5x**, still 55x
  below the noise, while collapsing `E_Q[u]` from 0.695 to 0.049. **94.8%** of `Var_R(u)`
  is the mosaic's between-position covariance, which no `beta_v` touches. So R3 is not the
  answer to R5. See `R5` and the corrections in `docs/memory.md` 2026-09-30 — the first
  version of this bullet claimed the path was bit-identical, which was true only of the
  one seed it checked.
- ⚠ **WE NOW HAVE OUR OWN chr21 VCF, AND IT IS NOT PANMIXER'S.** `vg deconstruct`
  at vg 1.68 gives **340,849 records** against the published 340,824, and **1,189,258**
  missing matrix cells against 1,278,653 — vg 1.68 calls 89,395 more genotypes and
  emits **zero** bare-`.` GT fields where 1.36 emitted 237,594. **Every figure in
  these docs must now say which VCF it belongs to.** The table in `docs/memory.md`
  2026-09-29 is the mapping. ⚠ **THE DERIVED FIGURES ARE SUPERSEDED BY DECISION C.**
  Under the coordinate cut chr21 had T = 305,886, dropping 34,371 below the span and
  592 above. **Under the declared-rate fill chr21's T is 340,848** — nothing is dropped
  for want of a cM, and only ONE position leaves the chain (`support` == 0). The three
  cM-source classes reconcile with the old cut exactly: 26,044 declared below the map
  + 592 declared above it + 8,327 interpolated plateau = 34,963 = what the cut dropped.
- **⚠ THE SCOPE CORRECTION OF 2026-09-28.** The v1 scope decision — *"maximise reuse
  of PanMixer, change exactly ONE thing"* — imported PanMixer's threat model, in
  which the target is a cohort member, into a design whose target is EXTERNAL to
  `G`, `D`, every panel and the attack database. An audit of the six docs plus the
  primer filed 91 findings; adversarial verification upheld 76 and found 51 more.
  Dylan's words: *"Boundaries aren't a weak function of the target's genotypes
  because the target is completely outside of any other dataset."* Corrected in
  commit `1f8ea3f` (`CLAUDE.md`, `terms.md`, `plan.md`, `data.md`). **When in doubt,
  the question to ask is: is this sentence about PanMixer, or about us?**
- **The exit-gate list itself had a hole, and it was in the headline clause.** Exit gate
  (b) is `tv_empirical <= tau` by exact enumeration. It holds and it is nearly useless:
  the worst exact TV on the toy reaches only **27.4% to 43.1%** of `tau`, so the bound is
  never approached. A mutation test settles what that costs — apply the tilt with
  `epsilon_tau` in place of `eta_tau`, twice the intended exponent, and the TV gate
  PASSES at every `tau` while `log_z_p` passes below `tau` = 0.5. A sixth clause,
  `ratio_spread <= exp(epsilon_tau)`, is now a PRIMARY gate: the normaliser cancels out
  of it, a constructed input pair SATURATES it exactly at every `tau`, and it is the only
  clause that fires. Defined in `docs/terms.md`.
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

## THE DECISIONS OF 2026-10-08 — the scope change, settled by Dylan

Nine decisions opened Phase B. They were taken after a seven-surface audit
(71 agents) of what the change implies; the audit's residual list is below.
**Decisions A-I are settled and are not to be re-litigated.**

| | decision |
|---|---|
| **A** | Primary pre-registered endpoint = the **genome-wide** utility, with expectations reported **EXACTLY** (derivatives of `log Z_p` / state marginals, no sampling). A contiguous-window utility is a clearly-labelled **secondary** arm. |
| **B** | `beta_v` sums to 1 **genome-wide over the 22 autosomes**; each chromosome takes a SHARE. Not per chromosome. This is what the claim rests on. |
| **C** | Regions where the map gives **no cM at all** are no longer excluded — they are sampled under a **declared constant recombination rate**. ⚠ **Interior zero-cM plateaus inside the mapped span are NOT overridden** (see the register). Supersedes R2. |
| **D** | **Autosomes only.** No chrX, no chrY. |
| **E** | Arm 1 uses **per-target leave-one-out, K = 86**, removing BOTH of the target's haplotypes. Applied by **masking two rows in memory** at release time, not by materialising 44 cohort trees. |
| **F** | PanMixer's starting data is **rebuilt on our GRCh38 build**. Both arms run on ONE consistent build — the thing their own release got wrong. |
| **G** | **3 seeds on a single configuration** as a tilt diagnostic only, not expected in the paper. The paper's attack numbers rest on **R = 1 release per (target, `tau`)**, which is what constraint 5 requires. |
| **H** | Arm 1 adopts **`utility_loss_panmixer` only** from PanMixer's utility suite. `af_loss`/`ld_loss` are not adopted — but as a **scope choice**, not because they are undefined for us. In Arm 1 they ARE defined, and `docs/terms.md` has been corrected accordingly. |
| **J** | The **`beta_share` allocation rule** (proposed by Dylan, 2026-10-10): a chromosome's share is its **variant count** over the genome-wide total, `T_c / T_total`. This is definitional rather than merely defensible — it makes mean `beta_v` per position `1/T_total` on EVERY chromosome, i.e. it is exactly what "normalise `beta_v` over the whole genome in one go" means, with the per-chromosome share as bookkeeping rather than a second reweighting. ⚠ `T_c` is defined on the FULL 88-haplotype cohort, fixed before any target is selected, so Arm 1's leave-one-out cannot move it: under LOO up to 5 chr21 positions of 305,886 (0.0016%) would drop to support 0, which would otherwise make `beta_share` weakly target-dependent. ⚠ Consequence: the shares cannot be written as 22 numbers until preprocessing runs, so the pre-registration pins the RULE. |
| **I** | Arm 1's linkage attack scores "self" from the target's **own assembly-derived path** (defined for all 44), not its 1000G panel row (defined for only 39 of 44). PanMixer's panel-row semantics run as a labelled parity column. |

### ✅ DECISION L (Dylan, 2026-10-11) — the exact KL is its OWN quantity, not an A2 substitute

**The question.** `scripts/exact_kl.py` computes `kl_exact` exactly and derives
`acc_ceiling` from it by Pinsker. Was that a substitution for pre-registered attacker
**A2**, which says "exact enumeration of TV and `ratio_spread` over all admissible input
pairs"?

**The answer: no, and A2's own positive control says so.** It reads "already passes" —
A2 was written as the TOY check, where all 2,304 supported paths can be enumerated, and
it stays satisfied there. At chr21 scale the support cannot be enumerated at all, so
nothing replaces A2. `kl_exact` and `acc_ceiling` are new names for new quantities
(`docs/terms.md`), reported alongside the pre-registered attackers rather than in place
of one.

⚠ **Why introducing it after the signature is not selective reporting.** It is
DETERMINISTIC given the model and the input pair — exactly one value, nothing to choose
among after seeing an outcome. A measured attacker's accuracy does not have that
property, which is why A1 had to be pre-registered and this does not. **No addendum is
required**, and this entry is the dated record of the decision.

⚠⚠ **PINSKER IS LOOSE AND THAT MUST BE SAID EVERY TIME.** `TV <= sqrt(KL/2)` is an
inequality, so the true ceiling is closer to 0.5 than `acc_ceiling` prints. The
looseness runs in the direction that FLATTERS US, which is why the script prints the
caveat on every run and `docs/terms.md` repeats it in both entries.

**What each piece now carries, which is the real gain:** A1 is a measured attacker at
chance (a LOWER bound on distinguishability, constraint 7); A2 is exact enumeration on
the toy, where it is possible; `kl_exact` / `acc_ceiling` is a theory ceiling at genome
scale covering attackers nobody has written. Three different kinds of evidence rather
than one stretched to cover all three.

⚠ A new GATE comes with it: `acc_ceiling` must not exceed `p_succ_bound`. If the
computed ceiling is looser than what Theorem 1 promises, the COMPUTATION is wrong — and
it is the computation to suspect, not the theorem.

---

### ✅ DECISION K (Dylan, 2026-10-10) — the untraversable blocks are EXCLUDED

**Was a blocker; resolved the same day.** The forward pass raised on the rebuilt chr21
at position 10,543,014: with ZERO genetic distance between two chain positions the
transition is exactly the identity, so a donor that becomes unavailable has nowhere to
go, and when all surviving mass sat on those donors the path dead-ends.

**THE RULE.** A maximal run of chain positions with zero genetic distance between
neighbours is a block inside which no path can switch donors. If no single cohort donor
has a called allele at EVERY position in the block, no path traverses it, and **the
whole block leaves the chain.** It reads only the genetic map and `support(v)`, both
functions of `D`, so it is fixed before any target is opened and standing constraint 1
is untouched.

**MEASURED, all 22 autosomes: 20 blocks, 172,142 positions, 0.772% of the chain.**
14 of the 20 overlap a centromere. chr9 loses the most (74,292 positions, 7.15% of its
chain), then chr20 (18,360, 3.58%) and chr21 (8,327, 2.44%); chr10, chr14 and chr22 lose
nothing. Of 52,870 zero-distance runs genome-wide only these 20 are untraversable — the
other 52,850 are traversed by at least one donor and are kept, emitted as verbatim donor
runs exactly as `docs/terms.md` says they should be.

⚠ **THE WHOLE BLOCK GOES, NOT THE MINIMUM NEEDED.** The cheaper alternative — keep the
sub-run the best-covering donor spans, which would have retained about 96% of these
positions — was REJECTED: in that sub-run a single donor is forced, so the mechanism has
no choice at all there, and those positions would enter the denominator of every
per-position utility figure while carrying no variability.

**Consequences, all propagated 2026-10-10:** genome-wide T is **22,134,562**, `mean
beta_v` is **4.5178e-08** (up 0.778%), and because decision J allocates `beta_share` by
variant count every chromosome's share moved — chr9 by **-6.43%**, chr20 by -2.83%.

---

### ⚠ STILL OPEN — do not write these as decided

| | question | blocks |
|---|---|---|
| 1 | ✅ **SETTLED 2026-10-10 (Dylan): ONE GLOBAL 1.2707 cM/Mb**, the autosomal mean mapped rate (3,543.8 cM over 2,788.8 Mb). Declared in `configs/inputs.chr21.yaml` under `declared_rate`, with a test asserting it stays equal to the measured value so the config and the maps cannot drift apart. Global rather than per-chromosome because the measurement reversed the earlier recommendation: the five acrocentrics hold 96% of the map-silent sequence and have the genome's highest mapped rates, measured entirely on their q arms. ⚠ A DECLARATION, not an estimate of the region it applies to. | — |
| 2 | Whether **interior zero-cM plateaus** are eventually overridden too. Deferred deliberately: `docs/terms.md` states P(switch) = 0 there is CORRECT, and `tests/test_toy_enumeration.py` enforces it as a Theorem-1 **support** property. Overriding it reverses a tested position, and is the most rework-costly item on the board. ✅ **NOW SIZED** (2026-10-10): **93.9 Mb / 3.27%** of autosomal bp, LARGER than the class the declared rate covers; chr9 alone 25.61 Mb, chr16 10.69 Mb; and **8,327 chr21 positions**, all with `support > 0`. | nothing yet |
| 2b | ⚠ **The N-GAP POLICY — UNDECLARED, no longer blocked.** The gap table is acquired (`data/raw/ucsc_gap/gap.txt.gz`, 12 KB). ✅ **MEASURED: 79.8% of the map-silent sequence is gap** (68.8 of 86.2 Mb), so accumulating the declared rate across it asserts 109.5 cM of which **87.4 cM sits on placeholder bp** — 2.5% of the autosomal genetic length invented from spacers. ✅ **And NO variant of ours falls inside any gap** (chr21: 0 of 340,849), so gaps hold nothing and skipping their spans changes **no position's chain membership** — it only fixes `Delta_x` between variants that straddle one, dropping the map-silent distance 109.5 -> 22.1 cM. Recommendation: SKIP gap spans. ⚠ Dylan's call; it is a pre-registration item. | preprocessing |
| 2c | ✅ **SETTLED 2026-10-10 as decision J** — `beta_share` is allocated by variant count, `T_c / T_total`, with `T_c` on the full 88-haplotype cohort. Supersedes the genetic-length recommendation made the same day, which was wrong: it would have made a position's weight depend on which chromosome it sat on. | — |
| 3 | ✅ **CLOSED 2026-10-10.** The **`tau` grid** is written into `paper/prereg_2026-10-10.md` section 1: 0, 0.05, 0.10, 0.25, 0.50, 0.90, 0.99, each with its role declared. Unblocked by `B-EVL-exact-endpoint`. ✅ **SIGNED 2026-10-10** (Dylan), `src/store.py:git_rev` b5a220f122ca, so the grid is FROZEN under constraint 6. | — |
| 4 | ✅ **CLOSED 2026-10-10.** **Arm 2's `beta_scope`** is pinned to the value "chain" in the pre-registration section 2, which closes R1 by DECLINING its fix, on R4's grounds. | — |
| 5 | ✅ **CLOSED 2026-10-10.** **R3** — `w_v = 1 / support(v)`, pinned and declared as INHERITED rather than argued for, which is the honest form. | — |
| 6 | ✅ **DROPPED 2026-10-10.** The **secondary window** is not pursued: the pre-registration section 0 drops the contiguous-window second arm. | — |
| 7 | ✅ **CLOSED 2026-10-10.** **Arm 2's attacker set** is pinned (A1-A4, each with a mandatory positive control) and the **target count** is all 10 — chr21 for the powered attack, genome-wide for exact utility. ⚠ **[CORRECTION - 2026-10-10]** This row said the repo has "no implemented attacker of any kind" and that `attacker_accuracy` "names a module that does not exist": `src/attacks.py` was built 2026-10-10 with the Neyman-Pearson optimal statistic, the exact ceiling, a BALANCED `src/attacks.py:empirical_accuracy`, the arm guard and the pre-registration lock, with `tests/test_attacks.py` over it. | — |

## ⚠ HANDOFF — read this first after a compaction (written 2026-10-12)

**State:** chunks 0-7 done. Working tree clean, 335 tests pass. Next is chunk 8.

**The three things a fresh session will not infer from the code:**

1. ⚠⚠ **PRIVACY IS SETTLED; UTILITY IS THE OPEN QUESTION.** See the Current state bullet.
   `acc_ceiling` = 0.501555 at `tau` = 0.99 (bound 0.9950) is exact and strong. Utility
   is now MEASURED genome-wide and exactly: `e_q_u` moves **+8.345e-08** across `tau` = 0
   to 0.99 (0.708967646 -> 0.708967729), `e_q_u_over_ceiling` 0.858617102 -> 0.858617203,
   HG00096, ARM 2, 22 autosomes, n = 1. **The MOVEMENT is the figure.** The LEVEL — 0.8586
   genome-wide, 0.6931 on chr21 — is reached at `tau` = 0, where the release carries zero
   target information, so it is NOT attributable to the mechanism. Never quote the level
   as utility. ⚠ And do not re-derive the movement from budget dilution: `mean beta_v`
   scaling predicts 72x, sum(`beta_v`^2) predicts 128.6x, measured 35.2x — both failed, so
   `Var_Q(u)` from PAIR marginals is the open measurement. Arm 1 utility is still
   UNMEASURED at any scope (chunk 8).
2. ⚠ **ARM 1'S CHAIN IS PER-TARGET** (consequence of decision K). Leave-one-out removes
   both of the target's haplotypes, which lowers `support(v)` and can EMPTY a
   zero-distance block's donor intersection that was non-empty for arm 2. So chunk 8
   needs 44 chains, 44 `T` values, 44 `beta_share` allocations, and the
   `untraversable_block` set recomputed per target. Reusing arm 2's mask reintroduces the
   dead end decision K removed.
3. ⚠ **ONE CHUNK-9 DECISION IS OPEN:** PanMixer applies leave-one-out at the SAMPLING
   step but not at PMI block scoring (`hmm.py:122` uses all 88 haplotypes while
   `hmm.py:42-43` holds the target out). Arm 1 must match it or diverge DELIBERATELY and
   say which; the two give different numbers and the head-to-head is only direct if the
   choice is declared. Recommendation on record: hold out at both steps, state the
   difference.

**Running jobs at handoff:** `scripts/exact_utility.py` in tmux session `eut`, logging
to `logs/B_exact_utility.md` and writing a status sentinel named exactutil.status under
logs/ ON COMPLETION (so that file does not exist while the job runs — which is why it is
named in prose here rather than as a path, since the contract check resolves paths and
would fail on a file that is correctly absent). Genome-wide exact `E_Q[u]` at `tau` = 0
and 0.99 for one target. It is the number the utility bullet is waiting on; chr1 took
51 s, so the 44 chromosome-passes are roughly half an hour.

**Standing hazards this session hit repeatedly, all in `docs/bugs.md`:** common random
numbers (THREE occurrences — a per-position tilt of ~1e-06 means two runs sharing a seed
ARE the same draw); `pgrep -f` matching your own command line (five occurrences); a
number read without checking what it is SCOPED to (CPU vs elapsed, system vs process);
and a self-computed checksum verifying storage rather than transfer.

---

## ⚠ THE CHUNK MAP — the review-gated plan, and where we are

⚠ **[CORRECTION - 2026-10-11] THERE ARE TEN WORKING CHUNKS, 1 through 10, plus chunk 0
as pre-work. An earlier version of this map had THIRTEEN, which I invented.** The error
was treating two stretches of forced work as NEW chunks and renumbering the planned ones
behind them. They were not new chunks — they were work INSIDE chunks 5 and 6:

- the `tau`-definition audit and the pre-registration SIGNATURE are **chunk 5's
  completion**, not a chunk 7. Chunk 5 was "decide the pre-registration's eight fields";
  signing it is what finishes that, and the audit reopened and re-closed field 6.
- the release writer and **decision K** are **chunk 6(d)**, not a chunk 8. Chunk 6
  always included "route releases through provenance"; decision K was forced by actually
  running it.

Renumbering the planned work is the worse half of the mistake: it moved the Arm 2 attack
from chunk 7 to chunk 10 in a document people read to know what is next.

One chunk = one review gate: nothing downstream starts until the chunk is reviewed.

| chunk | what it covers | state |
|---|---|---|
| 0 | Pre-work: clear the known-wrong facts out of the authority docs before phase B code reads them | ✅ done |
| 1 | `B-THY-diploid-gate` — a diploid, two-segment falsification gate | ✅ done |
| 2 | `B-IMP-genomewide` — `beta_v` sums to 1 genome-wide, per-chromosome share | ✅ done |
| 3 | `B-DAT-map-coverage` — how much of each autosome the genetic map speaks to | ✅ done |
| 4 | The declared rate settled; the exact endpoint at genome shares; gaps skipped; the arm guard; the toy genome | ✅ done |
| 5 | **The pre-registration**: eight fields decided with Dylan, audited, and SIGNED 2026-10-10 | ✅ done |
| 6 | **Genome-scale data, and releases that can be quoted.** (a) 22-autosome cohort build (b) 1000G 30x panel (c) genome-wide TARGET preprocessing (d) releases through provenance, which forced decision K | ✅ **4 of 4 done 2026-10-11.** 220/220 target pairs; chr10's panel was found CORRUPT by using it and refetched |
| 7 | **The Arm 2 privacy result**: A1's sampled attack AND the exact `acc_ceiling` | ✅ **done 2026-10-12.** A1 a clean null (17/18 cells below bound, none significantly above); `acc_ceiling` = **0.501555** at `tau` = 0.99 against a bound of 0.9950 |
| 8 | **Arm 1 leave-one-out releases** — 44 cohort members, diploid, 22 autosomes, K = 86 | ⬜ needs 6(c). ✅ **ALIGNMENT CONFIRMED 2026-10-11**: PanMixer obfuscates ONE individual per experiment, independently, against the ORIGINAL panel, and already does leave-one-out (`hmm.py:42-43`) — so our design matches theirs and is NOT order-dependent. ⚠ One difference to declare: their LOO applies at sampling but not at PMI block scoring (`hmm.py:122` uses all 88 haplotypes). ⚠ And see the per-target chain risk below |
| 9 | **The PanMixer head-to-head** — their own evaluations on our build, so the comparison is direct | ⬜ needs 8 |
| 10 | **The writeup** — two sections, one per arm | ⬜ |

⚠ **TWO STRETCHES OF WORK IN CHUNKS 5 AND 6 WERE UNPLANNED, and that is the structure
working rather than failing.** The audit began with a question about one sentence in the
proposal's introduction and ended by reopening a signed-off field; decision K began by
running a release and hitting a dead end no test could see. A plan that cannot absorb
those is a plan that would have shipped them.

⚠⚠ **A RISK DECISION K CREATED FOR CHUNK 8, recorded now so it is not discovered then.**
Arm 1 removes BOTH of the target's haplotypes from the donor panel, which lowers
`support(v)` and can EMPTY a zero-distance block's donor intersection that was non-empty
for arm 2. So **arm 1's chain is per-target**: 44 chains, 44 `T` values, 44 `beta_share`
allocations, and the `untraversable_block` set recomputed per target. Reusing arm 2's
mask there would reintroduce exactly the dead end decision K removed. See `in_chain` and
`untraversable_block` in `docs/terms.md`.

⚠⚠ **AND AN OPEN QUESTION ABOUT THE SIGNED PRE-REGISTRATION, raised 2026-10-11.**
Decision K (21:59) POST-DATES the signature, and it changed the released object: 172,142
positions left the chain, so `T` and every `beta_share` moved. Two consequences:
section 7 item 6 of the signed file says the interior zero-cM plateaus are "emitted as
verbatim donor runs", which is now only partly true — 20 blocks of them are EXCLUDED;
and `beta_v`'s normalising set, a pre-registered field, has different contents under the
same name (`chain`). ✅ **RESOLVED 2026-10-11 (Dylan): a dated ADDENDUM**, `paper/addendum_2026-10-11_decision_K.md`,
SIGNED, with the signed pre-registration left untouched. It states both corrections and
records that the frozen fields — `tau` grid, attacker set, endpoints, attack database,
decoy design, `delta` — are unaffected and that no attack outcome had been inspected.
⚠ It is deliberately NOT named `prereg_*.md`, because `require_signed_prereg` globs that
pattern and an addendum must never be able to become the document a result is stamped
with. Three tests pin that.

---

## ⚠ PHASE B WORK ORDER — what must happen before the first number

Ordered, and the order is load-bearing.

⚠ **THE CHUNK MAP ABOVE IS THE LIVE PLAN (reinstated 2026-10-11 at Dylan's request).**
This work order is the ORDERING CONSTRAINT inside it, not a second plan — every item here
belongs to a chunk, and the chunk is what carries the review gate.
⚠ **[CORRECTION - 2026-10-11]** A note here said the chunk labels were "now HISTORICAL
... and they are not extended", reasoning that two numbering systems where one is stale is
worse than one. That was the wrong fix for a real problem: the labels HAD gone stale —
chunk 5 was never archived under one, and the 2026-10-10 audit and signing carried no
number — but the answer is to maintain them, not to retire the frame Dylan asked for.
They are maintained again and the gaps are filled in the map above.

1. ✅ **`B-THY-diploid-gate` — DONE 2026-10-10.** `tests/test_diploid_gate.py`,
   64 clauses, mutation-tested against the real `ffbs.py` with each revert verified
   byte-identical. The budget invariant `sum beta == ploidy` is the clause that
   catches the genome-wide composition error, and **it is forward-looking**:
   `utility.beta_v` still hard-demands `sum == 1`, so the error cannot yet be
   injected into the mechanism. The gate will bite when step 5 implements the
   share — which is the point of building it first.
2. ✅ **`B-DAT-map-coverage` — DONE 2026-10-10.** `scripts/map_coverage.py`,
   output at `logs/B_map_coverage.md`. Map-silent 3.00% of autosomal bp, interior
   zero-cM 3.27%, autosomal mean mapped rate 1.2707 cM/Mb. ⚠ **The N-gap half is
   BLOCKED** — it needs the GRCh38 assembly's gap coordinates and no FASTA, AGP or
   assembly report is on disk. ⚠ Position counts exist for chr21 only; the other
   21 arrive with preprocessing, because they need each chromosome's own VCF.
3. ✅ **`B-EVL-exact-endpoint` — DONE 2026-10-10.** `logs/B_exact_endpoint.md`.
   The endpoint is monotone and exactly computable at every candidate share, so
   **the `tau` grid is unblocked.** ⚠⚠ But the route matters: the finite-difference
   derivative route returns a NEGATIVE variance at a genome-wide share, so the
   pre-registration must pin the endpoint as E[differ] via the **marginals**, and
   `Var_Q(u)` may no longer be quoted at genome scope from derivatives.
4. ✅ **`B-LCK-preregister` — DONE. ✅ **SIGNED 2026-10-10** (Dylan), `src/store.py:git_rev` b5a220f122ca.**
   `paper/prereg_2026-10-10.md`. Settled and written in: the declared rate (1.2707
   global, gaps skipped), the share allocation rule, the primary endpoint and the
   route it must use, the attacker set with mandatory positive controls, the arm
   guard, Arm 1's "self" source, and the framing (the toy is the correctness
   demonstration; the headline rests on the proof plus genome-scale measurement).
   ✅ **The secondary windowed-utility arm is DROPPED** — the toy genome
   demonstration does its job without a second release of the same person, which
   would have composed.
   ✅ **ALL EIGHT FIELDS FILLED AND ✅ **SIGNED 2026-10-10** (Dylan), `src/store.py:git_rev` b5a220f122ca.** Field 6 was reopened on 2026-10-10 —
   the pinned attack design was ONE-SIDED and could not measure the quantity
   `p_succ_bound` bounds — and re-closed the same day at the BALANCED 25-draw design
   (`delta` = 0.056 unchanged, 7,000 releases, ~18.3 h). ⚠⚠ The attacker set, the `tau`
   grid, the endpoints, the attack database and the decoy design are now FROZEN. A
   change means a NEW dated file; this one is never edited. ⚠ **[CORRECTION - 2026-10-10]** This line read "8 FIELDS REMAIN OPEN",
   which was true when written and was false from the moment the fields were filled.
   The load-bearing one is **the minimum detectable advantage the experiment is
   powered for**: the repository forecasts a chance-level null, and without a power
   statement a null is not evidence (standing constraint 7).
5. ⚠ **`B-IMP-genomewide` — ONE RESIDUAL LEFT: THERE IS NO RELEASE WRITER.**
   `src/mechanism/ffbs.py` does not import `src/store.py` at all, so a release lives
   only in memory; `Release.record()`'s docstring says "the arrays go to disk through
   `store.save()`", which describes what a CALLER is expected to do and not what any
   code does. Until that exists no genome-wide release is quotable, because an
   artifact with no provenance sidecar is not quotable by this project's own rule —
   and that lands squarely before item 7.
   ✅ The other residuals this item carried are closed: the per-chromosome SHARE is
   implemented and gated, the driver and the 22 manifests exist (item 6), and the
   `log_c` accumulator is hardened with exact summation.
   ⚠ **[CORRECTION - 2026-10-10]** this line read "the genome-wide mechanism (W1-W14
   below)". **There is no W1-W14 list anywhere in the docs** — a search of every
   markdown file finds the reference and nothing else, so it was lost in an earlier
   rewrite and the item has pointed at nothing since. A cross-reference to a list
   that does not exist is worse than no cross-reference: it reads as though the work
   is specified somewhere.
6. ✅ **`B-DAT-rebuild` — DONE 2026-10-10.** Decision F's consistent GRCh38 rebuild,
   all 22 autosomes on one axis: 22,307,009 records, **22,134,562** chain positions,
   driver exit 0 and verifier 22 of 22. ⚠ The chain figure is post-decision-K; the
   build itself produced 22,306,704 and decision K then removed 172,142.
   ⚠ **NAMING DRIFT, MINE, 2026-10-10:** I logged this as `B-DAT-genomewide-build`
   when `B-DAT-rebuild` already named it here. `CLAUDE.md` forbids inventing a synonym
   for something already named, and two names for one work item is how a ledger stops
   being searchable. The ledger row keeps its original text under the canonical name,
   cross-referenced.
7. Releases, then attacks.

**Needs no decision and reduces risk now:** harden the `log_c` accumulator
(`math.fsum`/Kahan — the naive running-float form breaches the `log_z_p`
tolerance only at genome scale, after passing every chr21 test); re-gate float64
against genome scale (the per-position tilt falls BELOW float32 eps under a
genome-wide `beta`, so a dtype downgrade would silently produce an untilted
draw); delete the decompressed GFA after deconstruct; fix the dead scratchpad
path in `scripts/churn_across_tau.py`.

## Headline result

### The mechanism, chr21, HG00096 — v1 diagnostic, 2026-09-30

Metrics are rows, `tau` values are columns. Emitted by
`scripts/run_mechanism_chr21.py`; full table at `logs/A_HG00096_chr21_v1.md`.
⚠ **A v1 NUMBER. NOT PRE-REGISTERED.** `A-LCK-preregister` has not run and R1-R5 are
unsettled.

| metric | `tau`=0 (floor) | 0.05 | 0.10 | 0.25 | 0.50 | 0.90 |
|---|---|---|---|---|---|---|
| `utility_retained` * | 0.6931 | 0.6931 | 0.6931 | 0.6931 | 0.6931 | **0.6932** |
| `u_ceiling` * | 0.8155 | 0.8155 | 0.8155 | 0.8155 | 0.8155 | 0.8155 |
| `utility_retained` / `u_ceiling`, at each `tau` | 85.0% | 85.0% | 85.0% | 85.0% | 85.0% | 85.0% |
| `attacker_accuracy` * | n/a — no attacker implemented; pre-registration must fix the attacker set first, and a weak attacker would prove nothing (constraint 7) | — | — | — | — | — |
| `p_succ_bound` * | 0.5000 | 0.5250 | 0.5500 | 0.6250 | 0.7500 | 0.9500 |
| `log_z_p` * | 0.0000 | 0.0348 | 0.0697 | 0.1774 | 0.3815 | 1.0226 |
| `log_z_p` as a share of `eta_tau` | 0 | 69.5% | 69.5% | 69.5% | 69.5% | 69.5% |
| `two_target_indistinguishability` * | n/a at this scale — plug-in TV over a support this size is biased UPWARD and cannot falsify anything. Exact on the toy; worst pair 43.1% of `tau` at `tau` = 0.9 | — | — | — | — | — |
| `target_fidelity` * | n/a — needs its own ceiling on this graph (~9.9% of a target's non-reference calls have no record in `G`, carrying 20.5% of its `-log f`). `u_path` is fidelity-on-what-the-graph-can-see and is NOT a substitute | — | — | — | — | — |
| released state positions differing from `tau`=0, per haplotype, mean of 10 seeds | 0 | — | — | — | — | **417.2** (8 of 10 draws move; range 0–1,641) |
| `donor_switches`, both haplotypes, mean of 3 draws | 13908 | 13908 | 13908 | 13908 | 13908 | **13909** |
| `eta_tau` | 0.000000 | 0.050042 | 0.100335 | 0.255413 | 0.549306 | 1.472219 |
| `epsilon_tau` | 0.000000 | 0.100083 | 0.200671 | 0.510826 | 1.098612 | 2.944439 |
| `af_loss` / `ld_loss` | n/a BY DEFINITION — cohort-distortion metrics, undefined for an external-target release that never edits the cohort | — | — | — | — | — |
| seconds per release, diploid | 9.4 | 9.4 | 9.4 | 9.5 | 9.4 | 9.4 |

**Provenance.** `logs/A_HG00096_chr21_v2.md` · `chr21/hprc-v1.0-pggb/K88/T305886` ·
44 individuals / 88 haplotypes · `indicator/w=inv_support/beta=chain`, `beta` over
305,886/305,886 chain positions · 3 draws per `tau`, seeds 1000-1002 with
`SeedSequence.spawn` streams · `baseline_log_norm` -584.8135 · `chain_span` PROVISIONAL
(12,968,320 < POS <= 46,680,243) · R1-R5 UNSETTLED.
⚠ **SUPERSEDES `logs/A_HG00096_chr21_v1.md`**, whose three "independent" draws per `tau`
were four CORRELATED haplotype draws — the diploid seeds were `seed` and `seed + 1`, so one
release's haplotype 1 was bit-identical to the next release's haplotype 0. Re-run
2026-09-30 after the fix: `utility_retained` 0.6940 -> 0.6931, `donor_switches`
14,169 -> 13,908, and `sd(u)` at `tau` = 0.5 falls from 2.389e-03 to **1.173e-03** against
a diploid prediction of `sd(u_hap)/sqrt(2)` = 1.422e-03 — so the spread now corroborates the
exact value instead of coincidentally matching the haploid one. `log_z_p` and `u_ceiling` are
unchanged, as they do not depend on the draw.

**COLUMNS** · *`tau`=0* is the floor: an untilted draw from `baseline_model`, a real
cohort haplotype mosaic carrying zero target information — which is why
`target_fidelity` has to be its own axis. · *0.05-0.25* is the defensible range;
`tau` = 0.10 already bounds equal-prior identification at 0.55. · *0.90* is past any
claim we would make and is included because a bound must be checked where it is loose.

**ROWS** · *`utility_retained`* moves in the 4th decimal only, and only at `tau` = 0.9
(0.6931 -> 0.6932); the exact expected-utility movement across the whole grid is 2.94e-06. ·
*`u_ceiling`* is 0.8155, not 1: deferred decision R1 measured rather than predicted. ·
*`log_z_p`* DOES move, exactly proportionally to `eta_tau` at a constant 69.5%, which is
the first-order identity `log Z_p ~ eta_tau * E_R[u]` with `E_R[u]` = 0.694 — so the tilt
is applied at the right magnitude. · *released positions differing* is the whole finding
in one row. · *seconds per release* 9.4 s for a diploid release over T = 305,886 at
K = 88, so compute is not the constraint and never was.

**SYNTHESIS.** The mechanism is correct, the guarantee holds, and the tilt is doing
visible work that buys nothing. It reweights 305,886 positions by about one part in a
million each (1.80e-06 at mean `beta`, `tau` = 0.5), and its entire budget summed over
those positions is `sum_v eta_tau * beta_v = eta_tau` = 0.549. That bounds the number of
INDEPENDENT flips at order `eta_tau` — but one flip re-seeds the whole downstream path
through the backward recursion, so the released path moves by a mean of 417 state
positions while expected utility moves 2.94e-06. The achievable gain is
`eta_tau * Var_Q(u)` = 2.222e-06 against a per-draw `sd(u)` of 2.011e-03.
⚠ **`beta_v` is not the lever either**, which is the correction of 2026-09-30:
concentrating it 3,000-fold buys 16.5x the variance and still falls 55x short, because
94.8% of `Var_R(u)` is mosaic covariance set by the chromosome's genetic length. What
remains is structural: `phi_v` must decompose per position for `ffbs` to be EXACT, and
that is what forces `u` to be a near-uniform average over T positions, which is what makes
`Var_R(u)` tiny. A utility with real variance would not factor, and an approximate draw
carries no bound — the reason the k-mer Jaccard was ruled out. **That tension is the thing
to take to the PI.** See R5.

### The ingest measurement that gates use of PanMixer's data

Metrics are rows; panels are columns.

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

`CLAUDE.md` carries the invariant; this is the full checklist it points to.

Two results are comparable only if their config stamps match on:
**graph+cohort id · arm · `tau` · utility id · `loo_id`**. `loo_id` is REQUIRED
throughout Arm 1 (decision E makes per-target K = 86 the build for a whole paper
section) and is literally `none` in Arm 2, where the target is external and a
genuine release holds nobody out.

1. Metrics are rows, configurations are columns. Report the full metric set every
   time — `n/a` with a reason, never an omitted row.
2. Every number carries its arm, its `tau`, its denominator and its `n`.
3. Row labels are the exact `docs/terms.md` identifier; gate metrics marked `*`.
4. Order rows: primary -> gates -> structure -> context -> demoted.
5. Every table carries a provenance line and is followed by one bullet per column,
   one per row, and a prose synthesis. Every number in prose traces to a cell.
6. Emit tables from the scoring script; never hand-assemble one.

**Floor and ceiling are columns, not prose asides.** The floor is `tau` = 0, an
untilted draw from `baseline_model`.
⚠ **The ceiling is `u_ceiling` — per target, per `beta` scope, per `loo_id`, and
below 1.** The earlier wording here, *"the ceiling is the target's own path,
`u_path` = 1"*, was CORRECTED in `docs/terms.md` on 2026-09-30 and survived here
as a stale copy; it holds only for a target with a call at every weighted
position, which no real target has. Reporting `utility_retained` against an
implicit 1 is what produced the original misreading. In Arm 1 the ceiling is a
function of the leave-one-out, so it differs per target.
⚠ **The two arms' ceilings are not comparable to each other** — a cohort member's
alleles are largely carried by the other 86 haplotypes, an external target's are
not. Compare within an arm, across mechanisms; never across arms.

⚠ **`eps_pmi` and `tau` are incomparable, and so are `capacity` and `tau`** (see
`docs/terms.md`). **The head-to-head aligns the arms on `target_fidelity`**, which
is pre-registered as the matching axis because it is computable with no attack
run. Empirical attack success is the OUTCOME, never the alignment.
⚠ **CORRECTED 2026-10-08.** This section previously read *"aligns the arms on
MEASURED empirical attack success only"*. That selects the comparison point after
inspecting attack outcomes, which standing constraint 6 forbids — the plan itself
carried the violation it warns about everywhere else. Both full frontiers are
still reported.

## ⚠⚠ DEFERRED BY DECISION — REVISIT BEFORE ANY RESULT IS PUBLISHED ⚠⚠

Three things are knowingly left wrong-or-suboptimal so v1 can run (Dylan, 2026-09-30).
None is a bug; each is a decision to carry a known cost now and pay it later. **A number
produced before these are settled is a v1 number and must say so.**

### R1. `beta_v` is NOT normalised over `readable` — 18.5% of the tilt budget is inert

**Carrying:** `beta_v` sums to 1 over all 305,886 chain positions, but only **253,846**
(83.0%) are readable, so `u` can never exceed **0.815**. We pay the full `tau` for 81.5% of
the leverage, and the released path at the other 17% is a pure prior draw with no reference
to the target at all.
**The fix, when we take it:** normalise over `in_chain & readable` — max `u` = 1.000, nothing
invented, and `readable.npy` is already emitted and verified target-independent.
**⚠ The catch that makes this urgent rather than cosmetic:** per standing constraint 6 the
utility function must be pre-registered BEFORE any attack outcome is inspected. Changing the
normaliser after seeing results is not available. **Decide before `A-LCK-preregister`.**
**⚠ And it does not generalise as written** — see R4.

### R2. ✅ SUPERSEDED 2026-10-08 BY DECISION C — `chain_span` is no longer a cut

**Settled.** Map-silent regions are no longer excluded: they are sampled under a
**declared constant recombination rate**, so the coordinate cut that this item was
carrying is retired. Two things survive it and must not be lost:
- ⚠ **The rate's VALUE is not chosen**, and it cannot be until `B-DAT-map-coverage`
  runs. Do not invent one. The reach to the acrocentric p-arms also means the
  declared rate multiplies **fabricated base-pair distance** in GRCh38 placeholder
  sequence; say so rather than letting it pass as a measurement.
- ⚠ **Interior zero-cM plateaus inside the mapped span are NOT overridden.** That
  is open question 2, deliberately deferred, and reversing it is the most
  rework-costly item on the board.
- The `support == 0` exclusion is unaffected and still cohort-only.

**The original item, preserved:**

**Carrying:** 89.74% of chr21; the p arm and 592 telomeric variants are not released at all.
**The fix, when we take it:** nothing cheap. Every map that covers the p arm interpolates
across a void (pyrho spends ~10 cM on two markers), the graph carries the sequence only on a
CHM13 backbone at 43.9% missingness, and T2T maps do not lift back monotonically. The full
option set with measurements is in `chain_span` in `docs/terms.md`.
**Who decides:** Dylan with his PI.

### R3. `w_v = 1/support(v)` is inherited, not argued for

**Carrying:** `support(v)` counts donors with a CALLED allele, so this upweights
**poorly-called** positions — which is not the same as **informative** ones. A site where 87
of 88 donors are called but the target carries a singleton is highly informative and gets low
weight. It came from PanMixer, where it was a *cost* weight doing a different job.
**Why it is survivable for now, measured:** on this cohort it is nearly uniform. Median
support is **88**; **98.6%** of chain positions have support >= 80 and hold **96.8%** of the
weight; the top 1% of positions by weight hold only 2.8%. The 88x ratio exists in principle
and applies to almost nothing.
**The fix, when we take it:** decide what `beta_v` should actually emphasise — probably
rarity (`-log f`) rather than call rate. Also a pre-registration item.

### R4. ⚠ NEW — `readable` is target-independent ONLY because our targets share one panel

Raised 2026-09-30. `readable.npy` is computed from the target VCF's **site list**, and it is
identical across our ten targets **only because all ten were cut from the same 1000G panel
and therefore share one site list**. For a genuinely novel genome — the real setting — each
target arrives from its own calling pipeline with its own site list, so `readable` would
**differ per target**, and normalising `beta_v` over it would be target-dependent after all:
the exact constraint-1 violation R1 is meant to avoid.
**So R1's fix works in the test harness and does not, as written, survive contact with a real
target.** Before adopting it, `readable` must be redefined as a function of something fixed
and public — `G` itself, or a named reference callset pinned in advance — never of the
target's own VCF.
**This is the clearest example so far of something that works in the test setup and breaks in
the deployed one.** Do not let it through on the strength of the ten-target check.

### R5. ⚠⚠ NEW AND THE MOST CONSEQUENTIAL — the tilt moves nothing at chr21 scale

Raised 2026-09-30, from the first mechanism run. **Measured: the released path is
bit-identical at `tau` = 0, 0.1, 0.5 and 0.9** for a fixed seed — 0 of 611,772 positions
differ — and `utility_retained` is 0.6940 at every `tau`.

**Carrying:** nothing, yet. This is not a decision that has been deferred; it is a finding
that has to be understood before `A-LCK-preregister` can mean anything.

⚠ **CORRECTED 2026-09-30, SAME DAY.** The paragraph above originally opened by claiming the
released path is "bit-identical at `tau` = 0, 0.1, 0.5 and 0.9 — 0 of 611,772 positions
differ". **That was false**, and true only of the single seed it checked: at ten
disjoint-stream draws, **8 of 10 move**, by a mean of 417 state positions and 54 allele
positions per haplotype. It was also measured under COMMON RANDOM NUMBERS, so agreement was
the expected outcome rather than a finding. **The finding's substance is unchanged and now
rests only on the exact derivative measurement below.** The honest statement is that the
tilt VISIBLY RESHUFFLES the released path and the reshuffling is UTILITY-NEUTRAL.

**It is not a bug, and four independent things say so.** `log_z_p` tracks `eta_tau` at a
constant **69.5%**, which is the first-order identity `log Z_p ~ eta_tau * E_R[u]` with
`E_R[u]` = 0.694. `ffbs` is proven to draw exactly from `Q^tau_p` on the toy, where the
per-position tilt is 1.28 at `tau` = 0.9 rather than 1.000005. Two independent routes to
`log_z_p` agree to 1e-12. And the observed seed-to-seed spread matches the `sd(u)`
computed from the second derivative of `log Z_p`.

**The arithmetic, exactly.** `sum_v beta_v = 1`, so the total tilt pressure summed over
positions is `sum_v eta_tau * beta_v = eta_tau` — 0.549 at `tau` = 0.5. Spread over
305,886 positions that is a per-position reweight of **1.80e-06**, and the expected number
of positions at which it changes the draw is of order `eta_tau`, which **cannot grow with
T**. Sampling-free from the derivatives of `log Z_p` at `tau` = 0.5:
`Var_Q(u)` = **4.046e-06**, `sd(u)` = **2.011e-03**, achievable gain
`eta_tau * Var` = **2.222e-06**. Signal to noise about **1:900**; resolving it needs
~**8.3e+05** draws at 9.4 s each, about 90 days for one target at one `tau`.

⚠⚠ **`beta_v` IS NOT THE LEVER EITHER — MEASURED, AND THIS KILLS THE OBVIOUS FIX.** The
first version of this block asserted "THE LEVER IS `sum beta^2`, NOT `tau`" and quoted
`sum_v beta_v^2` = 3.27e-06. Both were wrong. The real `sum_v beta_v^2` is **4.556723e-06**;
3.269192e-06 is `1/T` = `mean(beta_v)`, the value `sum beta^2` takes only for UNIFORM
weights — two statistics that coincide in the uniform case were conflated. And
`Var_Q(u)` measured from `d2 log Z_p / d eta2` under `beta` held uniform over the top-n
positions by inverse support:

| `beta_v` scheme | sum of squared weights | `Var_Q(u)` | gain at `tau`=0.5 | `E_Q[u]` |
|---|---|---|---|---|
| inv_support (current) | 4.557e-06 | 4.046e-06 | 2.22e-06 | 0.6951 |
| uniform over all T | 3.269e-06 | 4.203e-06 | 2.31e-06 | 0.7075 |
| top 100,000 by inverse support | 1.000e-05 | 1.140e-05 | 6.26e-06 | 0.6610 |
| top 10,000 | 1.000e-04 | 1.114e-05 | 6.12e-06 | 0.2926 |
| top 1,000 | 1.000e-03 | 1.113e-05 | 6.12e-06 | 0.0873 |
| top 100 | 1.000e-02 | 6.676e-05 | 3.67e-05 | 0.0494 |

Concentrating the budget **3,000-fold** raises `sum beta^2` **2,200x** and `Var_Q(u)` only
**16.5x** — the gain stays **55x below** the per-draw `sd(u)` of 2.011e-03 — while `E_Q[u]`
collapses from 0.695 to 0.049, because a utility defined on 100 positions is a different and
far weaker quantity. The reason is that **94.8%** of `Var_R(u)` is between-position
COVARIANCE from the mosaic's segment structure (`sum_v beta_v^2 p_v(1-p_v)` = 2.115e-07
against `Var_Q(u)` = 4.046e-06), and the segment count is set by the chromosome's genetic
length and `Ne*r`, which no `beta_v` touches.

**So R3 is NOT the answer to R5.** R3 is still worth settling on its own merits — call rate
is not informativeness — but it will not make the mechanism measurable.

⚠⚠ **WHAT IS ACTUALLY LEFT IS STRUCTURAL, AND IT IS THE THING FOR THE PI.** `phi_v` must
decompose along chain positions for `ffbs` to draw from `Q^tau_p` EXACTLY. That requirement
is what forces `u` to be a near-uniform average over T positions; that is what makes
`Var_R(u)` tiny; that is what makes the tilt powerless. A utility with real variance would
not factor — and an approximate draw carries no bound at all, which is exactly why the k-mer
Jaccard was ruled out as the in-sampler potential. **Exact sampling and a high-variance
utility are in tension, and the mechanism currently sits at the exact-sampling end.**

**The option set, with the first entry now MEASURED AND CLOSED:**
- ~~**Concentrate `beta_v`** on SCATTERED positions~~ — **MEASURED AND IT DOES NOT WORK.**
  16.5x the variance at 3,000x the concentration. ⚠ But see the next entry: the table above
  tested only scattered sets, and that confound was the whole answer.
- ✅ **CONCENTRATE `beta_v` ON A CONTIGUOUS REGION — MEASURED 2026-09-30 AND IT DOES WORK.**
  Same n, same `sum beta` = 1, only the layout differs: at n = 1,000 a contiguous window
  gives `Var_Q(u)` = 1.0565e-03 against 1.1132e-05 scattered, **94.9x**; at n = 100 it is
  **125.1x**. `Var_R(u)` is governed by how many independent mosaic SEGMENTS the utility
  spans, not by `sum beta_v^2` — scattered positions each sit in their own segment, so
  concentrating onto them changes nothing. **And there is an exact ceiling:**
  `gain/noise per draw = eta_tau * sqrt(Var) <= eta_tau/2`, since `u` in [0,1] bounds `Var`
  at 1/4. So the effect is ALWAYS sub-noise in one release — that is what a `tau` bound
  MEANS — and the only question is how many draws it takes to see. Our configuration needs
  **6,552,960** (the DIPLOID figure; 3,276,480 is the haploid one — corrected 2026-10-01); a contiguous window of 100 needs about **1,600** (4 hours); the theoretical
  best at `tau` = 0.5 is **53**. We sit 250x below the ceiling, not at it.
  ⚠ **The cost is what the release is FOR:** a windowed utility preserves target information
  in the window and is a pure prior draw everywhere else. That is a different product from
  whole-chromosome fidelity, and it is the PI's question, not `beta_v`'s.
- ⚠ **OPTION 2 BELOW WAS ALSO WRONG** and is corrected here. Releasing k windows each with
  its own `tau` composes and is forbidden. Defining ONE release whose utility sums to 1 over
  ONE window does NOT compose: the chain still runs genome-wide, one path is released, and
  the positions outside the window are a prior draw. A narrower UTILITY is not a second
  release. See `docs/memory.md` 2026-09-30.
- **Shrink the scope `tau` covers.** A `tau` per megabase rather than per chromosome
  raises every `beta_v` by the ratio of scopes — but that is k independent releases of
  one genome and they COMPOSE, which standing constraint 5 forbids. It is the same scope
  algebra as the genome-wide question already recorded, pointed the other way, and it is
  not free.
- **Accept it and say so.** The honest reading of the current configuration is that at
  `tau <= 0.9` the release IS a prior draw, the bound holds trivially, and the paper's
  contribution is the guarantee rather than the utility. That is a defensible position and
  it must be a deliberate one.
- **Model B as the fallback** gives T = 100,757 rather than 305,886, so a 3x larger
  per-position tilt. Three-fold against a shortfall of 900 — not a fix.

**Who decides:** Dylan with his PI, and BEFORE `A-LCK-preregister`. ⚠ Whether `Var_R(u)`
itself shrinks with T is NOT established and must not be claimed either way: the mosaic's
segment count is set by the chromosome's genetic length, not by T.

---

## Ledger — Phase A

| ID | Intervention | Endpoint | n | Result | Verdict | memory |
|---|---|---|---|---|---|---|
| B-ATK-arm-guard | Make an Arm 1 attacker unable to produce a number that reads as an Arm 2 privacy result | the ill-posed combinations must RAISE, not return | 23 clauses | **BUILT.** `src/attacks.py` refuses database-search and membership attacks in arm 2 — the first has no correct answer (target outside the DB by construction, so 0% is the experimental design), the second inverts (a donor mosaic identifies DONORS). Arm is part of the spec's identity; a stated positive control is a MANDATORY field, putting constraint 7 in the type. Also the Neyman-Pearson optimum from the proposal's own algebra, asserted to ATTAIN the `(1+TV)/2` ceiling exactly. ⚠⚠ **And `p_succ_bound` CANNOT detect a miscalibration** — a doubled exponent reaches 0.8412 against a 0.9500 bound at `tau` = 0.9, because TV reaches only 27-43% of `tau`. Third instance of a gate against a never-approached bound having no power | done | 2026-10-10 |
| B-DAT-gap-skip | Stop the declared rate manufacturing distance from placeholder base pairs | gap-free spans only, and no position's chain membership may change | chr21 + 4 edge cases | **SKIPPED, per Dylan.** the gap reader and the gap-free declared-distance helpers in `src/cohort/genmap.py`. chr21's 10.33 Mb map-silent lead is 3.87 Mb gap-free (37.4%), so declared distance drops 13.122 -> 4.913 cM; autosome-wide 109.5 -> 22.1 cM, of which 87.4 cM was placeholder. Safe because NO variant of ours is inside any gap (chr21: 0 of 340,849), so nothing enters or leaves the chain | done | 2026-10-10 |
| B-THY-toy-genome | Demonstrate what the mechanism DOES, exactly, on a 2-chromosome toy genome | the tilt visible, and the sampler matching the enumeration | 2 chr, 7 `tau`, 40,000 draws | **THE TILT IS VISIBLE.** Shares 7/9 and 2/9 by decision J; mean `beta_v` 0.111 against chr21's 3.27e-06. TV(Q, R_D) reaches **0.3090**, `E_Q[u]` moves **0.4502 -> 0.6477** against chr21's 2.94e-06, `Var_Q(u)` 7.3e-02 against 4.05e-06. `ratio_spread` equals `exp(epsilon_tau)` to the last digit at every `tau`. Same `ffbs` path, no special case — **so it is budget dilution that makes the chr21 effect invisible, not the mechanism.** ⚠ A correlated-seed bug in the demo was caught by its own joint assertion: every marginal was right and only the joint was wrong | done | 2026-10-10 |
| B-EVL-exact-endpoint | Re-measure the PRIMARY endpoint at a genome-wide share, since every figure on record used the whole budget on chr21 | monotone and exactly computable at the share that will be used | 3 shares x 7 `tau`, two independent routes | **MARGINALS ROUTE WORKS, DERIVATIVE ROUTE BREAKS.** Movement scales EXACTLY linearly in the share (5 s.f. at two shares); monotone at every share; 3.45e+09x float64 headroom at the smallest. ⚠⚠ But `Var_Q(u)` from the finite-difference second derivative goes **NEGATIVE** at chr21's genetic-length share and scatters over 5x at 1/22 — and every recorded `Var_Q(u)` came from that route. The pre-registered endpoint must name the marginals route; the derivative route survives only as a chr21-scale diagnostic. ⚠ New residual: the share ALLOCATION RULE is undecided and the candidates span 2.8x | done | 2026-10-10 |
| B-DAT-map-coverage | Measure how much of each autosome the map speaks to, before the declared rate is chosen | the three classes reported separately, and chr21's split reconciling with the record | 22 autosomes, bp; chr21 positions | **MEASURED.** Map-silent **3.00%** of autosomal bp (86.2 Mb), **96% of it the five acrocentric p arms**; interior zero-cM **3.27%** (93.9 Mb), LARGER than the class the rate covers, chr9 alone 25.61 Mb. Autosomal mean mapped rate **1.2707 cM/Mb**, per-chromosome spread **1.84x**. chr21: 26,044 map-silent + 8,327 in-span plateau = 34,371, RECONCILES. ⚠⚠ **A per-chromosome rate is LESS defensible than a global one** — the acrocentrics' rates are the genome's highest and are measured entirely on the q arm, so they would extrapolate a hot q-arm rate onto heterochromatin. ⚠ N-gap policy BLOCKED: no reference assembly on disk | done | 2026-10-10 |
| B-EVL-exact-utility | Measure the PRIMARY utility endpoint genome-wide and exactly, since Dylan asked whether utility had been measured at all and every figure on record was chr21 at `beta_share` = 1 | the movement, not the level; and the completeness of the fan-out asserted | 44 exact passes; 22 autosomes x 2 `tau`; n = 1 target | **THE MOVEMENT IS +8.345e-08 ACROSS THE WHOLE `tau` RANGE.** `e_q_u` 0.708967646 -> 0.708967729 from `tau` = 0 to 0.99, `u_ceiling` 0.825708741, `e_q_u_over_ceiling` 0.858617102 -> 0.858617203. Zero Monte-Carlo error, no finite difference — the marginals route, because the derivative route breaks at genome scale. ⚠⚠ **TWO FIRST-ORDER PREDICTIONS OF THE MOVEMENT BOTH FAILED, IN OPPOSITE DIRECTIONS**: `mean beta_v` scaling said 72x, sum(`beta_v`^2) scaling said 128.6x, measured **35.2x**. So budget dilution alone does not set the movement; positive correlation along the copying path is the untested explanation and needs `Var_Q(u)` from PAIR marginals. ⚠ **The level is NOT the result** — 0.8586 is reached at `tau` = 0, where the release carries zero target information. ⚠ I introduced a definition drift and caught it: the ratio row was labelled `utility_retained`, which `docs/terms.md` defines as the SAMPLED `u_path`; renamed `e_q_u_over_ceiling`, and `--from-json` added so the relabelled table regenerated from saved results rather than 30 min of recomputation. ⚠ Two stale "not yet implemented" claims about the genome-wide normaliser retired — decision J built it on 2026-10-11 | done | 2026-10-12 |
| B-EVL-arm2-ceiling | Put an exact ceiling on every Arm 2 attacker, since sampling provably cannot | the ceiling must respect `p_succ_bound`, and the computation must reproduce its own second-order scaling | 70 cells; 10 targets | **0.501555 AT `tau` = 0.99 AGAINST A BOUND OF 0.9950.** No attacker, written or not, exceeds **50.16%** at the loosest grid point; 50.01-50.03% across the operating range. ⚠⚠ Validated by a check I did NOT build in: `kl_exact / eta_tau^2` constant to **0.31%** over a 53x range of `eta_tau`, the second-order behaviour theory requires. Implied Var(u_p - u_q) = 5.4664e-06, consistent with the recorded 4.05e-06 at a different scope. Gate passes 70/70. ⚠ The theorem is loose by a factor of ~300 at genome scale — R5 stated exactly rather than by failing to detect it. ⚠ Pinsker flatters us; quote the caveat | done | 2026-10-12 |
| B-DAT-targets | Map every prepared target onto every autosome's site axis, and give the target VCFs a producer | 220/220 pairs with provenance, and the paths must be distinct data rather than copies | 220 pairs; 10 targets | **CHUNK 6 IS 4 OF 4.** 187,268,720 called cells over a genome-wide `in_chain` of 22,134,562, mean fraction **0.846047**. ⚠ The target VCFs had NO PRODUCER: chr21's were hand-made without a sidecar, so the only target set that existed was irreproducible. ⚠⚠ **A CORRUPT PANEL SURFACED BY USING THE DATA, NOT BY CHECKING IT** — chr10 arrived with a CRC32 mismatch at offset 393,732,438 that passed every check in `scripts/verify_panel.py` INCLUDING the sha256, because the sidecar is computed FROM the download; the fresh copy is the SAME SIZE with a different hash. Fixed with `bgzip -t` in both the fetcher and the verifier. ⚠ Two false readings of mine in one pass: a monitor reported a MISSING BINARY as a corrupt file, and I explained a uniformity that 3-decimal formatting had invented | done | 2026-10-11 |
| B-ATK-arm2-run | Run the pre-registered Arm 2 attack and establish whether its null is trustworthy | the exact route must agree with the sampler, and every cell must have its own stream | 1,750 trials; 5 hypotheses; 3,320 validation draws | ⚠⚠ **FIRST RUN INVALID — COMMON RANDOM NUMBERS, THIRD OCCURRENCE.** Seeded from the draw index alone, so 10 targets x 7 `tau` shared 25 seeds; at a per-position tilt of **7.7e-07** the released path is **bit-identical** at `tau` = 0, 0.05 and 0.25, so the seven columns were ONE experiment. Clustered SEs move z from -3.81 to **-1.99**. ✅ The IMPLEMENTATION is correct: exact `E_Q[Lambda]` = **+1.784e-07**, positive for all ten targets, sampler at **z = -0.16** over 3,320 draws, `log_z_p` bit-equal. ⚠ Three of my own claims were wrong, including a first-order `eta_tau` "fingerprint" that was an artefact and that I used to narrow the search. ⚠ Decoy selection also ran through Python's salted `hash()` and was NOT reproducible. ⚠ Cost **110 CPU-hours** against the prereg's 18.3 — it omitted the 6 forward passes per trial. ✅ Both fixed, re-run launched with 1,750 distinct streams | rerunning | 2026-10-11 |
| B-FIX-untraversable | Clear the blocker: make the rebuilt chain releasable | the release must RUN, and the rule must be recomputable independently of the code that applied it | 20 blocks; 15 tests; 3 mutations | **DECISION K, AND THE RELEASE NOW RUNS.** chr21 at T = 332,521, baseline log_norm -879.4501, 9.1 s per diploid draw. Of **52,870** zero-distance runs genome-wide only **20** are untraversable — **172,142 positions, 0.772%**, 14 of 20 overlapping a centromere. ⚠ Not a pairwise property: every adjacent pair can share a donor while no donor spans the run, and THREE donors are the minimum to construct that case. ⚠ Propagated: T 22,306,704 -> **22,134,562**, `mean beta_v` -> **4.5178e-08**, and every `beta_share` moved (chr9 -6.43%). ⚠⚠ Two holes in my own checks, both found by the opposite-direction test: the verifier asserted `in_chain` was exactly `support > 0`, and its replacement only recomputed the blocks when something was ALREADY excluded — so a chain never rebuilt PASSED, the one case it exists for | done | 2026-10-10 |
| B-IMP-release-writer | Make a release quotable: route it through the provenance chokepoint, and refuse every wrong write | removing any gate must fail its own test, and the path must work end-to-end on real data | 20 tests; 7 mutations | **BUILT, GATED, AND IT CLOSES ITEM 5's LAST RESIDUAL.** `src/releases.py`. 7 of 7 gates fire when removed, control clean; verified end-to-end on real chr21 arrays with the axis digest, `src/store.py:git_rev` and five input digests reaching every sidecar. ⚠ The gate worth having: **genome scope is checked against the 22-autosome set**, because a `GenomeRelease` of one chromosome with share 1.0 satisfies every in-memory check — shares sum to 1, no duplicates, both bounds hold — so nothing before the write could tell it from a release of a genome. ⚠ My own fixture was caught by `GenomeRelease` doing its job: a flat per-chromosome `log_z_p` over 22 parts summed to 1.05 against an `eta_tau` of 0.1003. ⚠⚠ **AND RUNNING IT FOUND A BLOCKER**: the rebuilt chr21 dead-ends the forward pass on an interior zero-cM plateau, 9,971 such steps genome-wide. Item 7 is blocked on that decision, not on this | done | 2026-10-10 |
| B-LCK-preregister | Freeze `tau`, the utility, the endpoints, the attacker set, the attack database and the decoy design before any outcome is inspected | the lock must OPEN, and reverting the signature must fail a named test | 8 fields; 4 mutations | ✅ **SIGNED 2026-10-10** (Dylan), `src/store.py:git_rev` b5a220f122ca. All eight fields filled, zero `[OPEN]` markers, eight ticked boxes; `attacks.require_signed_prereg` now returns rather than raises. The live lock test's expectation FLIPPED in the signing commit, which is what its own docstring demanded. ⚠⚠ **FIELD 6 WAS REOPENED AND RE-CLOSED THE SAME DAY**: the drafted design drew only from `Q_p`, so it could not measure the equal-prior quantity `p_succ_bound` bounds — it would have reported a VIOLATION with the theorem intact. Re-closed BALANCED at 25 draws/target per `tau` per hypothesis, `delta` unchanged at 0.056 (the balanced estimator is more precise per observation), 7,000 releases, ~18.3 h. ⚠ Signing would also have frozen three stale copies of the superseded 300-draw design, including two ADJACENT TICKED BOXES pinning 300 and 50 | done | 2026-10-10 |
| B-FIX-tau-definition-audit | Audit every statement of the `tau` guarantee against the proposal's abstract, after Dylan found the introduction defining it differently | each finding must survive two independent adversarial refutations | 6 sources, 36 candidates, 78 agents | **20 CONFIRMED, AND ONE WAS A LIVE FALSE-POSITIVE RISK ON A PRIVACY CLAIM.** ✅ The DEFINITION is fine: the abstract says TV <= `tau` and attacker <= (1+`tau`)/2, `docs/terms.md` matches it, and no file states the loose form. ⚠ The proposal's INTRODUCTION is wrong as written — "distinguished with the probability of no more than `tau`" is impossible for `tau` < 0.5, since a coin gets 0.5; it means ADVANTAGE, which equals TV. Needs a one-line fix before submission. ⚠⚠ `scripts/arm2_attack.py` compared a ONE-SIDED rate to a bound that only covers the balanced quantity — it would have reported a violation with the theorem intact. ⚠⚠ The pre-registration's field 6 is REOPENED: the pinned design draws only from `Q_p` and so cannot measure what it is compared to. ⚠ Three stale copies of the superseded 300-draw design, including TWO ADJACENT TICKED BOXES in the signature checklist pinning 300 and 50. ⚠ 12 stale doc claims, and the contract's absence-marker skip had no converse — a known hole recorded 2026-10-08 and left open | done | 2026-10-10 |
| B-DAT-panel | Acquire the 1000G 30x panel for all 22 chromosomes and prove it is on OUR reference build | the build must be checked against our own deconstruct output, and a GRCh37 length must be refused | 21 chromosomes; 9 tests | **COMPLETE AND BUILD-VERIFIED.** 21 files, 29.3 GB, **69,689,262 records**, **3,202 samples on every chromosome**, zero failures, no leftover partials; plus chr21's pre-existing panel, so all 22 are covered. ⚠⚠ **The build check compares the panel's contig lengths against OUR OWN pangenome's**, not against constants I typed — a wrong-build panel produces a low match rate that reads as biology, which is exactly how decision F's hazard would present. All 21 verified. ⚠ Truncation is caught by SEEKING to the far end, because a cut-short download still opens and still indexes against its own partial content. Deleting the build comparison fails exactly its two tests. ⚠ Checksums are RECORDED but not re-verified (`--deep` re-reads 29 GB) | done | 2026-10-10 |
| B-DAT-rebuild (logged 2026-10-10 as `B-DAT-genomewide-build`, a synonym I should not have invented) | Build all 22 autosomes on one consistent axis, and make a manifest unable to claim a rule its artifacts do not support | driver exit 0 over all 22, verifier exit 0, and chr21 reconciling with the superseded cut | 22 chromosomes; 8 checks each | **ALL 22 BUILT AND PASSING.** 22,307,009 records, **22,306,704 chain positions**, only **305 excluded genome-wide**. ⚠⚠ **`mean beta_v` = 4.4830e-08, 73x smaller than chr21-alone's 3.27e-06** — the dilution is now measured, not projected, and every chr21-scope per-position figure shrinks by it. Declared rate spans **434x**, top four are four of five acrocentrics. chr21 reconciles with the old cut **to the position** (26,044 + 592 declared + 8,327 interpolated = 34,963 = exactly what the cut dropped). ⚠ **The manifest described the CONFIG, not the BUILD** — chr21's 2026-09-29 arrays were relabelled `chain_rule` = "declared_rate_fill" because every resume sentinel was satisfied and only the manifest re-ran. Now the claim is CHECKED against `cm_source.npy`, whose absence the coordinate cut cannot fake; fired on the live case. ⚠ chr21's old tree moved to `DEPRECATED_chr21_coordinate_cut` — its LD blocks carried the OLD axis digest | done | 2026-10-10 |
| B-FIX-resume-deconstruct | Make the expensive step resumable, and stop the verifier failing a correct chromosome | recovery must reproduce the originally recorded digest exactly, and a tie must not be reported | 3 real chromosomes; 14 tests | **BOTH FIXED, AND ONE WAS A FALSE POSITIVE OF MINE.** (1) `deconstruct` was missing from the resume map, so any late failure cost a whole chromosome rebuild (measured 437.9 s, ~2 CPU-hours); it is now resumable AND the skip branch repairs missing provenance rather than carrying the gap forward. chr10 recovered in **9 s** against a measured 437.9 s full rebuild, and the recomputed graph digest matches the originally recorded one **byte for byte on chr9, chr16 and chr19** — so the recovery is identical provenance, not an approximation. (2) ⚠ The verifier's position check demanded STRICT increase and failed chr9, which has one tied position where a parent bubble (`lv` 1) and a nested 1 bp deletion (`lv` 2) share a start. Nesting is routine at 10-18% of records per autosome, so the tie is legal output; the check now forbids only a DECREASE and reports ties as a column | done | 2026-10-10 |
| B-FIX-deconstruct-prov | Stop the deconstruct step deleting the graph it is about to checksum | reinstating the original ordering must fail a named test | 6 tests, 3 fail on the reverted bug | **FIXED, AND THE CAUSE WAS MY OWN MID-RUN COMMIT.** `8d56346` inserted the GFA cleanup ABOVE the provenance call, so `src/store.py:sha256_file` was handed a path that had just been unlinked — after vg had spent 1-2 h. Six chromosomes lost (chr10-14, chr17), three more in flight. ⚠⚠ **THE COMMIT LANDED WHILE THE BUILD WAS RUNNING**, and each chromosome imports `vcf.py` in a fresh subprocess, so the split is by PROCESS START TIME: chr9's manifest is stamped 32 min AFTER the commit and is fine, while everything started after it failed. ⚠ No existing test could see it — they all passed a plain `.gfa`, so the unlink never fired. ⚠ And the traceback named the provenance call, not the cleanup three lines above. Recovery is cheap because vg writes the VCF BEFORE the crash: all 9 survived, and `src/cohort/vcf.py:attach_deconstruct_prov` plus `src/store.py:sha256_gz_plain` rebuild the sidecar without re-running vg, saving ~12 h | done | 2026-10-10 |
| B-DAT-cohort-verify | Check the built cohort tree is CORRECT, not merely complete | every injected corruption must be caught, and the unmutated control clean | 10 mutations; 13 tests; 5 built chromosomes | **BUILT, AND IT HAS POWER — 10 of 10 caught, control clean.** The driver asserts COMPLETENESS (every chromosome ran, every file present and non-empty); a row-misaligned per-site array passes all of it and then scores the wrong variants, because the arrays are positional and a shift yields numbers rather than an error. The check with the real power is `support` against BOTH the called-genotype count and `INFO/AN` — two independent routes to one quantity, so they can disagree. chr1/2/4/16/19 pass every check. ⚠ **chr21's tree is PRE-DECISION-C** — no `cm_source.npy`, `chain_rule` absent — so it is not comparable with the other 21 and must not be quoted genome-wide; `_build_one` never skips, so its batch rebuilds it. ⚠ Two defects in the checker itself, both found by the mutation test and not by use: a missing artifact raised instead of reporting, suppressing all 22; and the length check fell through into positional checks that raise on broadcast | done | 2026-10-10 |
| B-IMP-genomewide (part) | Make `beta_v` sum to 1 genome-wide with a per-chromosome SHARE, and make the error that permits detectable | 5 injected errors must each fire a named clause | 25 clauses; backward compat on real chr21 | **SHARE IMPLEMENTED, EXACTLY BACKWARD COMPATIBLE.** Weights at a share are bitwise that share times the whole-budget weights; at `tau` = 0.5 the ratio of the two `log_z_p` values equals the share exactly, and `log_z_p` is 69.46% of its own cap under both, reproducing chr21's recorded 69.5%. `baseline_log_norm` and `u_ceiling` match the record to the decimals they were quoted at; exact summation moves the normaliser by 1.1e-13. ⚠⚠ **3 of 5 mutations survived the first 20-clause suite** — loosening an assertion changes no value, so a test that independently checks the same bound cannot see it. Also the aggregate's share check was hidden behind its own normaliser check. Both fixed, all 5 now fire. ⚠ No driver, no manifests, releases still bypass provenance | partial | 2026-10-10 |
| B-THY-diploid-gate | Falsify the DIPLOID and multi-segment release, which every prior gate missed | injected calibration errors must fire a named clause at every `tau` | 64 clauses, 3 mutations, 6 `tau`; two-segment sweep EXHAUSTIVE over 1,024 pairs | **BUILT AND HAS POWER.** The composite primitives `exact.py` lacked, plus one invariant that catches both composition errors: `sum beta == ploidy`. ⚠ The first version checked `exact.py` against itself, passed 57 clauses and would have caught NEITHER doubled-exponent mutation — both live in `ffbs.py`. The clauses with the power are the two tying reference to mechanism, and the strong one needs no draws: `log_z_p` is deterministic. ⚠ TV fires on nothing; the bracket alone catches a doubled exponent only at `tau` = 0.9 | done | 2026-10-10 |
| B-FIX-chunk0-stale-facts | Clear wrong facts from the authority docs before phase B code reads them | each distinctive stale number greps clean; historical copies carry a correction pair | 5 files of live text, 4 figures | **4 corrections + 1 new measurement.** 1.61e-06/13.5x used the site axis not the chain (-> 1.80e-06/15.1x); 21 TB is 21.1 GB, wrong by 1000x in 5 places; 15.03% is one draw, exact is 14.8077%. ⚠ The contract's correction-format check matched only column 0, so it had skipped every correction written inside a bullet — strengthened, and it immediately caught 3 | done | 2026-10-10 |
| A-EVL-churn-across-tau | Does raising `tau` move a few positions or reshuffle many? | the mismatch SET, not just its size, at a fixed seed | 2 models, 7 `tau`, 3 seeds, common axis | **Both, and still tiny.** At `tau`=0.99 model C moves 40.7 positions to net −12.7 and model B 84.7 to net +23.3 — about 3x more movement than the net. Below `tau`=0.25 model C's churn is EXACTLY ZERO. ⚠ The changes CLUSTER: 12 of C's 23 and 51 of B's 64 fall in one 400-variant window, because a flip swaps a whole donor segment. Sampled nets are noise; exact is −3.30 (C) and −4.72 (B) | done | 2026-10-08 |
| A-IMP-model-b | Build the LD-block chain so switching happens only between blocks | model C reproduced exactly by identity blocks; blocks from a panel that excludes the target | chr21, 8 capacities swept on the other arm | **BUILT.** 59,472 units, `expected_switches` 3,394.1 against 3,532 predicted independently; floor 14.5699% vs C's 14.8077%, slope −1.7738 vs −1.2381 — model B is a MILD contiguity lever, 1.43x, against 16x from re-scoping the utility. Cost: 1,222 positions dropped, 1,212-variant verbatim donor runs. ⚠ One bug, caught only by the prediction: transitions were start-to-start, not inter-unit gaps | done | 2026-10-07 |
| A-EVL-capacity-sweep | Does PanMixer have a knob, and what does it reach? | its reachable `target_fidelity` range | 8 capacities, HG00438, seed 123 | **YES — `capacity`, and it SATURATES.** 1.20% at 0.01 rising to **16.41% at 0.2-0.3**, then four settings byte-identical. Our 14.81% is reached at **capacity ~= 0.17**, correcting the 0.45 interpolated from two points | done | 2026-10-07 |
| A-FIX-mechanism-audit | Independent audit of the same day's mechanism build; docstrings declared non-evidence | every gap between docs and code, adversarially verified | 127 agents, 6 dimensions, 3 refuters per gap | **6 real defects the 88-test suite passed over**, 4 confirmed by mutation and all fixed; 30 gaps survived of 40 filed. Overturned the day's headline row and closed R5's leading option | done | 2026-09-30 |
| A-THY-toy-enumeration | Brute-force correctness gate on a hand-written graph, built BEFORE the sampler | sampler matches exact enumeration; exact TV <= `tau`; `log_z_p` bracketed; `tau`=0 reproduces `R_D` | K=4, T=7; 96 output paths; 3,571,128 input pairs per `tau` | all clauses PASS. ⚠ AND A SIXTH WAS NEEDED: exact TV reaches only 27.4-43.1% of `tau`, so a doubled tilt exponent passes it at every `tau`; `ratio_spread` is sharp and is the only gate that fires | done | 2026-09-30 |
| A-IMP-cohort-hmm | `baseline_model` from the cohort alone, constants re-derived from the Methods | the initial distribution and the `stay_switch` transitions built from `G` and `D` only, with no target reachable | chr21, K=88, T=305,886 | built; no target parameter, asserted from the AST; E[switches] **7,021.0** vs the **7,155** model C predicted; O(TK) agrees with dense O(TK^2) to **1.11e-16** | done | 2026-09-30 |
| A-IMP-utility | Bounded utility: `phi_v` an indicator, `beta_v` normalised inverse-support weights | `u_path` in [0,1] and `beta_v` summing to 1, both ASSERTED at runtime | chr21 | built; `u_ceiling` = **0.815455** over 253,846/305,886 positions, independently reproducing R1's predicted 0.815 | done | 2026-09-30 |
| A-IMP-ffbs-sampler | Algorithm 1: stored alphas, exact backward sample, no argmax | draws from `Q^tau_p`; `log_z_p` range-asserted; `A` never materialised | 40,000 draws/`tau` on the toy; 18 releases on chr21 | built; output distribution AND state-path joint both within Monte-Carlo error; 9.4 s per diploid chr21 release. ⚠ `log_z_p` needs the `baseline_log_norm` subtraction — see `docs/bugs.md` | done | 2026-09-30 |
| A-EVL-thin-tilt | Does the tilt actually change the released path at chr21 scale? | a measured utility difference between `tau` = 0 and `tau` = 0.9 | 1 target, 6 `tau`, 10 disjoint-stream seeds + exact derivatives | **The path MOVES (8/10 draws, mean 417 state positions per haplotype) and the movement is UTILITY-NEUTRAL**: gain `eta_tau * Var_Q(u)` = 2.222e-06 against `sd(u)` = 2.011e-03, about 1:900. Not a bug. ⚠ Reweighting `beta_v` cannot fix it — 3,000x concentration buys 16.5x the variance, 94.8% of which is mosaic covariance. ⚠ The first verdict here read "NO — 0 of 611,772 positions differ" at n=3 under common random numbers; corrected same day | done | 2026-09-30 |
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

**Provenance.** The five 2026-09-30 rows are `tests/test_toy_enumeration.py`,
`tests/test_mechanism.py` and `logs/A_HG00096_chr21_v1.md`, all reproducible from a clean
checkout (the toy rows need no `data/` at all). The ingest row's evidence is workflow
`wf_a30688a6-3d7` plus hand verification recorded in `docs/memory.md`.

## In Progress

⚠ **Both sections below are now SPECIFICATIONS OF BUILT THINGS, not work in flight.**
They are kept here rather than compressed into the ledger because each carries decisions
and measurements that the ledger's one-line rows cannot hold, and because the deferred
items R1-R5 above point into them. The live work is `[A-LCK-preregister]`, which is
BLOCKED on R5.

### [A-IMP-preprocess] THE PREPROCESSING PIPELINE — specified 2026-09-28, BUILT 2026-09-29

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
   the same person and would compose — exactly what standing constraint 5 forbids. This is
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

#### ⚠ TWO THINGS THE FIRST RUN SURFACED THAT NEED A DECISION (2026-09-29)

**0. THE TOOL IS INDEPENDENT OF PANMIXER (Dylan, 2026-09-29).** PanMixer is in this
repository as a published example of a problem in the same space — a different threat model,
a different release object, a different attack vector — and must not be overweighted in what
we build. Concretely, as of 2026-09-29: our own conda env (`environment.yml`, name
`ranpanmixer`) with nothing taken from their `environment.yaml`; `src/` importing only
stdlib, numpy and pyyaml; our own `vg deconstruct` producing our own VCF; and the
development fixture sliced from **our** VCF rather than theirs. **No PanMixer-derived
artifact is an input to anything we build.** Their allele indices differ from ours anyway
(see `docs/memory.md` 2026-09-29, the ALT-reordering finding), so a cross-build join would
be wrong as well as unwanted.

**1. `absent_policy` — 16.4% of the tilt budget is inert (was 19.0% before the split-ALT
fix of item 1d), and there is a better third option.** The mechanics, stated precisely because the phrasing matters:

- Our site axis is the **pangenome's** 340,849 records. The target is read from the 1000G
  30x panel, which has its own site list. **24.1% of our records have no 1000G record at
  all** — graph-derived structural and nested variants that a short-read SNV/INDEL/SV
  callset never called. The direction is graph → panel, not panel → graph.
- At those positions we do not know the target's allele, so `path.npy` holds `-1`. `phi_v`
  is then undefined, the position contributes nothing to `u`, and the sampler emits a pure
  prior draw there. **The release carries a cohort allele at those positions with no
  reference to the target whatsoever.**
- Because `beta_v` is normalised over **all** chain positions, those inert positions sit in
  the denominator. So `u` can never exceed **0.836**, and the tilt tops out at
  `exp(0.836 · eta_tau)` instead of `exp(eta_tau)` — we pay the full `tau` for 83.6% of the
  leverage. Measured: 260,308 of 305,886 chain positions have a target call, identical for
  all ten targets because it is a property of the two site sets, not of any individual.
  Note the inert positions are **disproportionately weighted** — 16.4% of the weight against
  14.9% of the positions — because graph-only variants sit in harder regions with lower
  support and therefore higher `w_v = 1/support(v)`.

| option | effect |
|---|---|
| "missing" (current) | honest, but max `u` = 0.836; 16.4% of the budget cannot move |
| "reference" | fills allele 0. Recovers the budget but **invents a call** — absence from a callset means the caller reported no variant there, NOT that this person is reference. Worse than inert: `phi_v` would then reward donors matching a fabricated value |
| **normalise `beta_v` over the panel-matched positions** | max `u` = 1, nothing invented. **Legitimate, because the matched set is a property of the PANEL and is fixed before any target arrives** — identical for all ten, so it is target-independent and can be declared as part of the utility definition |

⚠ **CORRECTION 2026-09-29 (independence audit).** The justification below was attached to
the WRONG set. `in_chain & (target_reason == 0)` — positions where a particular target has a
call — is **NOT** invariant: measured, **4 distinct sets across the ten targets**, because a
contradiction at a poly-A repeat drops a row for one target and not another. Normalising
`beta_v` over it would be exactly the constraint-1 violation this section warns about.
The set that IS invariant is now emitted as **`readable.npy`** — a row is readable if the
panel has a usable record for it, by strict match or full multi-allelic coverage, and both
tests read only the panel's **site list**, which is the same file for every target.
**Verified: 1 distinct set across all ten**, and locked by a test. **Normalise `beta_v` over
`in_chain & readable`**, which gives max `u` = 1.000 over **253,846 of 305,886** chain
positions (83.0%). Normalising over all chain positions caps `u` at **0.815**.

**Measured 2026-09-29, and it settles the target-dependence worry.** With `target_reason.npy`
in place the 19% splits as: **45,578 positions code 1** (no record in the target callset,
16.4% of the weight) and **0 positions code 2** (record present, no call). A phased reference
panel has no missing genotypes within its own site list, so every target is fully called
there — identical for all ten. **The restricting set therefore provably contains no
per-person variation**, which is exactly the property the third option needs.

⚠ The third option is **not** the same as "normalise over positions where this target has a
call" — that would make `beta` target-dependent, which is exactly the constraint-1 violation
that voids Theorem 1. The distinction is whether the restricting set is the panel's or the
person's. **Recommended, not yet decided.** Also note this is partly an artifact of using a
panel as the target source; a real external genome called against GRCh38 with full coverage
would have a different unmatched fraction.

**1b. Three distinct losses, not one.** They are easy to conflate and the numbers differ:

| | direction | measured |
|---|---|---|
| target variant with **no graph record** — unrepresentable at any `tau` | target → graph | **10.2%** of non-ref calls, **21.5%** of information |
| target variant where the graph has the POS but a different REF/ALT spelling | target → graph | ~9.9%, joinable only under a relaxed key we deliberately do not use |
| **graph record with no target record** — we cannot read the target there | graph → panel | **14.9%** of chain positions, **16.4%** of `beta_v` weight (was 17.6% / 19.0% before the split-ALT fix) |

The first two are ceilings on fidelity. The third is a ceiling on the *budget*. Reporting
one when you mean another inverts the conclusion.

**1d. ⚠ PART OF THE INERT BUDGET WAS A BUG, NOW FIXED.** The 1000G panel is decomposed into
biallelic records while 7.9% of our records are multi-allelic, so a strict full-ALT join
could never match those — 8,377 records, 8,283 of them in the chain, whose target alleles we
actually held. `split_alt_join` recovers them by matching each of our alts against its own
panel record. Chain positions with a target call went **82.4% -> 85.1%**, and representable
non-ref calls **80.0% -> 83.5%**. So the remaining inert fraction is **16.4%** of the weight
(over **14.9%** of chain positions), not 19.0% over 17.6%. ⚠ 3 recoveries in 10 targets hit a **poly-A repeat** where the panel's records are
not mutually exclusive; those rows are left unresolved rather than guessed, and the
contradiction counter is what makes that visible.

**1c. `support(v) = 0` is now excluded from the chain (2026-09-29), not patched.** Measured: of the
**305,886** chain positions, only **252,025** have a target allele at all, because 24.1% of
pangenome records have no matching record in the 1000G panel the target comes from. Under
the default `absent_policy = "missing"` those positions have no `phi_v`, so they can never
tilt — they contribute to the `beta_v` normaliser and never to utility. **19.0% of the
budget is spent on positions that cannot move.** Identical across all ten targets, because
it is a property of the two site sets rather than of any individual.
The alternative, `absent_policy = "reference"`, fills them with allele 0 and recovers the
full budget — but it INVENTS a call: absence from a panel means the site was not in the
callset, not that the sample is reference there. Both are implemented and recorded in each
target's provenance. **Not decided.** The honest framing is that "missing" understates
achievable utility and "reference" overstates confidence in the target's path.

There is exactly one such site on chr21, at **pos 8,460,681**. It was outside `chain_span`
by luck rather than design, and the cut is provisional. `chain_mask` now excludes
`support == 0` outright.
⚠ **Why exclusion rather than `w_v = 1/(support(v) + 1)`**, which was the obvious patch:
`support(v) = 0` means **no donor haplotype has a called allele there**, so the sampler has
nothing to emit and RENORMALISE has an empty pool. The position is unusable, not merely
awkward. `1/(support+1)` would keep it in the chain and hand it a finite `beta_v` share it
can never earn — re-creating the inert-budget problem of item 1 one position at a time, and
hiding a position that cannot produce output. The exclusion is a function of the cohort
alone, so it is target-independent.

**3. The target now has its own reason array (2026-09-29).** A `-1` in `path.npy` had two
meanings that `phi_v` cannot distinguish, because both simply mean "no tilt here":
`no_record_in_target_callset` (our join found nothing — a limitation of the join, identical
for every target from the same panel, and no statement about the person) and
`record_present_but_no_call` (the callset covered the site and reported no genotype for this
individual — a real failure to call them). They behave identically in the mechanism and
differently in the reporting: a position the panel never covered arguably should not count
against `target_fidelity`, while one it covered and could not call arguably should. The
cohort's codes cannot express either, since those are about graph nesting. Emitted as
`target_reason.npy` beside `path.npy`.

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

### [A-IMP-v1-sampler] v1 SAMPLER SPECIFICATION — pinned 2026-09-25, NARROWED 2026-09-28, ✅ BUILT 2026-09-30

**Scope decision (Dylan, 2026-09-25; NARROWED 2026-09-28).** The original wording was
*"maximise reuse of PanMixer and change exactly ONE thing, the sampler — keep its blocks,
its data model, its donor panel and its whole evaluation suite."* That went too far: it
imported PanMixer's threat model, in which the target is a cohort member, into a design
whose target is external to `G`, `D`, every panel and the attack database. **Narrowed:**
we reuse PanMixer's DATA and its forward recursion. Blocks as a chain unit, the
donor-panel-minus-target, and the cohort-distortion metrics are COMPARISON-ARM choices,
not properties of our design. **The chromosome-wide chain is no longer deferred — model C
is the architecture (decided 2026-09-28). ✅ MODEL B IS NOW BUILT TOO (2026-10-07):
`hmm.from_arrays(..., blocks=...)`, with model C asserted identical to its identity case,
so the fallback is one argument away rather than a rewrite.**

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
per-position tilt at tau = 0.5 is **1.80e-06 under C**, only **15.1x** the noise floor,
against 5.45e-06 and **45.7x** under a block partition.
⚠ **CORRECTED 2026-10-10.** This read **1.61e-06** and **13.5x**, computed over the
340,849-record SITE AXIS rather than the chain's T = 305,886. The chain is the
denominator — `eta_tau * mean(beta_v) = eta_tau / T` — so the correct figures are
0.549306 / 305,886 = 1.80e-06 and 15.1x. The old pair understated the tilt's headroom
by 10%. ⚠ **And under decision B this margin goes the other way:** with `beta_v`
summing to 1 genome-wide the per-position tilt falls BELOW float32 eps entirely, so
float64 stops being load-bearing-by-10% and becomes the only thing preventing an
untilted draw. See `docs/bugs.md`. Rounding would consume a meaningful
fraction of the target's entire influence. Alpha storage is **240 MB** under C and **71 MB**
under the model-B fallback, both at K = 88; a K = 86 head-to-head figure is ~2.3% smaller
and must say so. ⚠ Break the chain at genuine discontinuities — outside `chain_span`, the
centromere, long assembly gaps — never at block boundaries.

**EXIT GATES** — ✅ **ALL MET 2026-09-30.**
1. ✅ `A-THY-toy-enumeration` passes: sampler matches exact enumeration (output
   distribution AND state-path joint); exact TV <= `tau` across a grid of `tau` and all
   3,571,128 input pairs per `tau`; `log_z_p` within [0, `eta_tau`] on every draw;
   `tau` = 0 reproduces `baseline_model` to 1e-15.
   ⚠ **AND A SIXTH CLAUSE, because clause (b) above turned out to be the weakest member
   of this list.** Worst exact TV reaches only 27.4-43.1% of `tau`, so a tilt applied
   with `epsilon_tau` in place of `eta_tau` passes it at every `tau`.
   **(f) `ratio_spread <= exp(epsilon_tau)`**, which a constructed pair saturates
   exactly and which is the only clause that catches that error.
2. ✅ One code path, and under model C there is no special case to test: every chain
   position is a variant, so no singleton/T = 1 branch exists at all. What the gate still
   tests is that **no allele-frequency branch exists in the source** — now a test that
   scans every module in `src/mechanism/`. (The original
   wording was written against the block-structured spec, where singletons were the
   special case; retained here only for the model-B fallback.)
3. ✅ **REWRITTEN AND MET 2026-09-30.** The replacement is `donor_switches` against the
   prior expectation and the analytic segment length: measured **7,021.0** expected
   switches per haplotype on chr21 at K = 88, against the **7,155** the model-C decision
   predicted from an independent calculation, and 14,169 realised for a diploid release.
   The original wording is preserved below because it is a good example of a gate that
   would have rejected a correct model.
   ⚠ **MIS-SPECIFIED — must be rewritten before it is used as a gate** (see
   `docs/memory.md` 2026-09-27, A-DAT-genetic-map). As written it demands switch mass of
   order 0.1 to 1 between adjacent chain positions, but those sit about 12 bp apart, where
   the biologically correct value is the measured 0.0026 — the gate would reject a correct
   model. Its second clause, "a sampled block is demonstrably NOT a single donor copied
   verbatim", fails for about 91% of blocks even with the correct map, because an LD block
   is by definition a stretch with little recombination and v1 chains within one block.
   Replacement candidates, all measured 2026-09-27: expected donor switches per chromosome
   (3,622.6, against PanMixer's ~0); median P(mosaic) among blocks with >= 20 variants
   (0.2747); or measured segment length against the analytic n/(Ne*r) = 0.0068 cM.
4. ✅ The missing-data policy is chosen, implemented, and covered by tests that fail if
   the behaviour changes — including one asserting a dead-ended chain raises rather than
   returning a path, which is what forced RENORMALISE to be a potential rather than a
   renormalised transition row.
5. ✅ Every chain position carries a genuine interpolated cM — no clamped values survive
   into the chain, and `switch_mass_from_cm` refuses a NaN or a non-monotonic map with a
   named failure. Interpolate with `left=nan, right=nan` and assert no NaN remains
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

### [A-THY-toy-enumeration] Brute-force correctness gate — DONE 2026-09-30
- **Why:** the ONLY test that can falsify our implementation of Theorem 1.
- **Built:** `src/mechanism/toy.py` (K = 4, T = 7, reads nothing from disk) +
  `src/mechanism/exact.py` (the O(K^T) reference, which imports neither `hmm.py` nor
  `ffbs.py` — asserted, because a check sharing code with what it checks cannot fail)
  + `tests/test_toy_enumeration.py`, 15 tests.
- **Built BEFORE the sampler**, which is what the prerequisite list got wrong: clauses
  (b), (c) and (d) need only the DEFINITIONS of `R_D` and `u_path`, so the gate
  falsifies Theorem 1 with no sampler in existence. Only clause (a) needs `ffbs`.
- **Result:** all clauses pass. 16,384 state paths -> 2,304 with mass -> 96 output
  paths. `ffbs` matches enumeration on the output distribution AND the state-path
  joint. Exact TV <= `tau` over all 3,571,128 input pairs per `tau`.
- ⚠ **A SIXTH CLAUSE WAS NEEDED, and it replaces (b) as the headline.** Worst exact TV
  reaches only **27.4-43.1%** of `tau`, so the bound is never approached and cannot
  detect a mis-scaled tilt. Mutation-tested: applying the tilt with `epsilon_tau` in
  place of `eta_tau` PASSES clause (b) at every `tau` and passes (c) below `tau` = 0.5.
  **(f) `ratio_spread <= exp(epsilon_tau)`** — the normaliser cancels out of it, a
  constructed pair saturates it exactly at every `tau`, and it is the only clause that
  fires. See `ratio_spread` in `docs/terms.md`.

### [A-IMP-cohort-hmm] The target-independent prior — DONE 2026-09-30
- **Built:** `src/mechanism/hmm.py:CohortHMM`. No function in the module takes a
  target-shaped argument and it imports neither `utility` nor `ffbs`; a test reads the
  AST and fails on either.
- **Constants re-derived from the Methods** (Ne = 10,000, r = 1.26, `Delta_x` in cM,
  `d = Delta_x * Ne * r`), never from PanMixer's code. The rank-one decomposition
  `(1-s) + s/K` is cross-checked against the published closed form
  `e^(-d/n) + (1-e^(-d/n))/n` at construction, to 1e-15, at K = 4/44/86/88.
- **The old gate wording is retired, as flagged.** "Switch mass O(0.1-1) at realistic
  separations" would have REJECTED a correct model — adjacent chr21 variants sit about
  12 bp apart. Replaced by `donor_switches` against the prior expectation and the
  analytic segment length: measured **7,021.0** expected switches per haplotype at
  K = 88, against the **7,155** the model-C decision predicted independently.
- ⚠ **RENORMALISE is a 0/1 POTENTIAL, not a renormalised transition row.** Across a
  `Delta_x` = 0 step the transition is the identity, so a donor that becomes
  unavailable has a masked row summing to exactly zero — no chain to renormalise.
  **1,058 zero-distance steps sit inside chr21 `chain_span`**, so this is real. `R_D`
  is therefore the Li-Stephens chain CONDITIONED on availability, which is also the
  form this spec already pre-registered. Consequence: `baseline_log_norm`.

### [A-IMP-utility] Bounded utility — DONE 2026-09-30
- **Built:** `src/mechanism/utility.py`. `phi_v` an indicator, `beta_v` normalised
  inverse-support weights, `u_path` and the `beta_v` sum both ASSERTED at runtime
  rather than documented. `diploid_utility` is the MEAN, with a test that catches a SUM.
- **Support is the FULL cohort**, K = 88 — our target is external, so nothing is held
  out. `support(v)` is read, not computed, and cross-checked against `INFO/AN`.
- ⚠ **`u_ceiling` is NEW and is now PRIMARY.** `utility_retained` was defined against a
  ceiling of 1. On a real external target under `beta` scope `chain` it is **0.815455**
  over 253,846/305,886 chain positions — deferred decision R1 measured, not predicted.
  Report both or neither; see `docs/bugs.md`.
- **The `beta_v` scope is an explicit argument** with no silent default: `chain` (the
  current, undecided-safe choice) and `chain_readable` (R1's fix, with R4's warning in
  the docstring). It is stamped on every release as part of `utility_id`.

### [A-IMP-ffbs-sampler] Algorithm 1 — DONE 2026-09-30
- **Built:** `src/mechanism/ffbs.py`. Every alpha column stored in float64 (215 MB at
  chr21 scale), tilt potentials multiplied in, per-position normalisation, `log_z_p`
  range-asserted, no argmax anywhere — the last checked by walking the AST, because the
  first version of that test grepped the source text and fired on the docstrings
  explaining why there is no argmax.
- **`A` is never materialised.** Verified by peak-allocation measurement at K = 3,000,
  not by assertion: a correctness test cannot see this, since a materialised matrix
  would give identical numbers and simply not run at scale.
- **O(TK) agrees with the O(TK^2) dense reference to 1.11e-16** on the toy and 1e-12 at
  K = 60, T = 40 on random data.
- ⚠ **`log_z_p` = the tilted `sum_t log c_t` MINUS `baseline_log_norm`.** Skipping the
  subtraction puts it at **-584.81** on chr21 against a bracket of [0, 0.5493], and the
  tempting fix of widening the bracket would have retired the cheapest correctness
  assertion in the project while leaving it in place as decoration. `docs/bugs.md`.
- **A diploid release is ONE release**: the mean utility factorises exactly, so the two
  chains are independent at HALF the exponent and `log Z_p` is their sum, landing in
  [0, `eta_tau`]. One private input, one `tau`, one bracket.
- `new_haplotypes` for the comparison arm is NOT emitted yet; it was never part of this
  gate and belongs to `[A-EVL-headtohead]`.

### [B-THY-diploid-gate] A diploid, TWO-SEGMENT falsification gate — ✅ DONE 2026-10-10

⚠ **[CORRECTION - 2026-10-10]** This section's heading read "— THE TOP PRIORITY" and
its body claimed "no gate in the repo can currently detect" the full-`eta_tau`-per-
segment composition error. Both are now false: `tests/test_diploid_gate.py` has 64
clauses, is mutation-tested against the real `ffbs.py`, and the invariant
`sum beta == ploidy` catches BOTH composition errors. The body below is retained as
the record of what the gate was built to do.

⚠ **RE-SCOPED 2026-10-08.** Decision E makes Arm 1 diploid AND decision B makes it
multi-chromosome with a per-chromosome budget share, so a diploid gate on a single
chain with `beta` summing to 1 no longer covers the build. The gate needs **two toy
segments with declared shares `w_1 + w_2 = 1`, both haplotypes** — the existing toy
has 96 output paths, so the product support stays enumerable. It must also catch the
genome-wide composition error: running each segment at the FULL `eta_tau` rather
than its share, which is the same shape as the diploid defect below and which **no
gate in the repo can currently detect.** Added to `[A-EVL-headtohead]`'s
prerequisites, where it was missing. It is the cheapest blocking item on the
board: it reads nothing from disk and needs no rebuild, download or rate decision.
- **Why:** ⚠ **every falsification gate in this project is HAPLOID, and the diploid path is
  the only one that has ever been run.** All 18 releases of the first chr21 run carry
  `ploidy: 2`. Exact TV, `ratio_spread`, the sharp pair, the state-path joint, clause (f)
  and the KILL CRITERION all run on a single chain. Named by the audit's completeness critic
  as the single most important thing still missing, and it is right.
- **Why it is not just a missing test:** `src/mechanism/exact.py` has NO diploid primitive
  — no product-support enumeration, no diploid `utility_vector`, no diploid
  `release_distribution`. The gate cannot be applied to the diploid path because there is
  nothing to apply it with. This is missing REFERENCE CODE, not a missing assertion.
- **What is already patched:** the diploid `log_z_p` bracket, which was absent — a doubled
  diploid exponent was ACCEPTED at every `tau`, returning a record whose published
  `log_z_p` violated its own stated invariant. Fixed and regression-tested by mutation.
  ⚠ **A patched assertion is not a gate.** The bracket catches one error class; the toy
  gate catches the ones nobody predicted.
- **Prerequisite:** none, and ⚠ **MEASURED 2026-09-30: the toy does NOT need shrinking** —
  this entry first claimed the product support is K^(2T) and would need K = 3, T = 4. Wrong:
  the support is pairs of OUTPUT paths, and the toy has 96 of those, so **9,216 ordered
  pairs** — trivially enumerable on the existing fixture.
  **The hard part is the INPUT space, not the support.** 2,673 distinct haploid inputs give
  7,144,929 diploid inputs and 5.1e+13 pairs of them, so the exhaustive all-pairs sweep the
  haploid gate runs is not available. The gate needs CONSTRUCTED extremal pairs instead —
  `toy.sharp_pair` generalises, since a diploid pair differing on both haplotypes everywhere
  still drives `u_p - u_q` to +/-1 — plus a seeded random sample, with the sample size and
  the seed recorded so the coverage is stated rather than implied.
- **Exit gate:** exact TV <= `tau` and `ratio_spread <= exp(epsilon_tau)` over pairs of
  DIPLOID inputs; `log_z_p` in [0, `eta_tau`] for the pair; `tau` = 0 reproduces the product
  baseline; and the mutation (full `eta_tau` per haplotype rather than half) fails at least
  one clause at EVERY `tau` on the grid, not only above 0.25.

### [A-FIX-audit-gaps] What the 2026-09-30 audit left open — needs Dylan's decisions
- **Why:** each is a real gap the audit confirmed, and each needs a choice rather than code.
- **The four that change behaviour:**
  1. ⚠ **`missing_policy`'s guard deviates from the recorded decision.** `docs/terms.md`
     says "Declare the fallback (uniform over all donors), count how often it fires, and
     report it." The code RAISES instead, on the reasoning that `support == 0` is excluded
     upstream so it cannot fire and a loud failure beats a silent fallback. **That is the
     implementer overriding a written decision and it needs Dylan's sign-off either way.**
  2. **Duplicate-haplotype collapse is not implemented** — no collapse, no policy, no
     assertion that duplicates are absent. ✅ **MEASURED 2026-09-30 AND IT IS BENIGN HERE:
     all 88 haplotypes are distinct** over the 305,886 chain positions, and the closest pair
     differs at **32,741** of them (10.7%). So every K = 88 this project quotes is correct;
     what was missing is that nothing said so and nothing checked. Downgraded from a
     correctness gap to a missing ASSERTION plus a label — but the assertion is worth adding,
     because both the uniform initial distribution and the per-state switch probability
     `s_t` divided by K are both functions of K, so a duplicate would put twice the prior mass on one
     haplotype, and a different cohort or a shorter `chain_span` could produce one.
  3. **Standing constraint 5 is unenforced.** No release cache, and `docs/data.md`'s
     "(target_id, tau, utility_id, seed) ... collisions must raise, never overwrite" exists
     in no code. `--seeds` is labelled an evaluation device but nothing structurally stops
     it being used as a release procedure.
  4. **`graph_id` defaults to the empty string** and is hand-assembled by the caller, so the
     only cohort stamp on a release can be silently absent or wrong. It should be derived
     from the cohort manifest's `axis_digest`, which already exists.
- **The rest, all cheap:** `hmm.expected_switches` has no test and no caller, so the 7,021.0
  figure quoted in five places is produced by no code in the repo; `ratio_spread` is pinned
  as a gate to be marked `*` in tables and appears in no table; `src/cli.py` has no
  mechanism subcommand, so the mechanism is reachable only through a script its own
  docstring calls "operational, not mechanism"; `scripts/toy_chain_models.py` is an orphan
  model-B relic referenced by nothing, and is where the wrong `sum beta^2` figure came from;
  and the post-split-ALT figures in this file's item 1 / 1b / 1d were reported inconsistent
  with the code's output and are UNVERIFIED by me.
- **Scope limits to disclose in any writeup:** one target of ten prepared, one chromosome
  of 22, and no attacker of any kind.

### [A-LCK-preregister] Freeze the evaluation before looking at it — NEXT, BUT BLOCKED
- **Why:** the proposal requires the utility function be chosen before any attack
  outcome is inspected; otherwise every later privacy claim is unfalsifiable.
- **Prerequisite:** ✅ `A-THY-toy-enumeration` passing (2026-09-30).
- ✅ **UNBLOCKED 2026-10-08 by decisions A and B**, which settle R5's question: the
  genome-wide utility is the PRIMARY endpoint, and expectations are reported EXACTLY
  from the derivatives of `log Z_p` rather than estimated from draws. R1's scope
  question survives for Arm 2 only (open question 4) and R3 survives untouched
  (open question 5) — both are still pre-registration items that cannot be
  revisited afterwards.
  ⚠ **The sentence this block used to open with was RETRACTED on the day it was
  written** — see the correction at R5 below. The released path is NOT bit-identical
  across `tau`; 8 of 10 disjoint-stream draws move. The item's substance (do not
  pre-register a grid the mechanism provably does not vary over) is unchanged and
  rests on the exact derivative measurement instead.
- ⚠ **NOW BLOCKED ON SEQUENCING, NOT ON A DECISION.** Decision B gives each
  chromosome a SHARE of the budget, so the endpoint's movement is far smaller than
  the chr21-normalised figure on record. `B-EVL-exact-endpoint` must run first, or
  the grid is pinned in ignorance of whether the endpoint moves over it at all.
- **Renamed `B-LCK-preregister`.** ONE dated file at `paper/prereg_<date>.md`
  covering both arms — ✅ **WRITTEN 2026-10-10 as `paper/prereg_2026-10-10.md`**,
  all eight fields filled, ✅ **SIGNED 2026-10-10** (Dylan), `src/store.py:git_rev` b5a220f122ca. ⚠ **[CORRECTION -
  2026-10-10]** This read "(PLANNED — it does not exist yet)" while the same file's
  work order already named the drafted document.
- **Exit gate:** a committed, dated file pinning, **per arm**: the `tau` grid, the
  utility id and `beta` scope, the primary endpoint, the attacker set, **which
  attack database**, and — per decision C — the declared recombination rate.

### [A-EVL-headtohead] The comparison on PanMixer's data
- ✅ **FIRST STEP IDENTIFIED 2026-10-03, and it is cheap.** `target_fidelity` — the share
  of the target's called alleles the release alters — is computable on either arm's output
  and is the axis to align on. Measured: ours is **14.8077%** EXACTLY — 75,177.401 of
  507,692 target-known cells, zero Monte-Carlo error — and effectively FLAT across the
  whole `tau` grid (the exact movement from `tau` = 0 to 0.99 is **3.28 cells**).
  ⚠ **CORRECTED 2026-10-10.** This read **15.03%**, which is 76,283/507,692 from a
  SINGLE DRAW at seed 1000 — the same statistic on the same denominator, but sampled,
  and it ran 1,106 cells high. `docs/data.md` records that the exact value supersedes
  every sampling estimate of this quantity. The `capacity ~= 0.17` crossing below is
  unaffected: the capacity sweep already used the exact 14.81%, so only this line was
  out of step.
  PanMixer's is **10.50%** at capacity 0.1 and **16.41%** at capacity 0.5. **Their capacity
  dial spans our fixed point**, so a matched-fidelity comparison exists at **capacity ~= 0.17** (measured 2026-10-07 by sweeping all eight capacities; the earlier 0.45 interpolated across a FLAT region from only two points) — which is the alignment this item needs and did not previously have.
  ⚠ Three obstacles, all quantified in `docs/memory.md` 2026-10-03: 10.25% of their records
  sit outside our `chain_span`; the two VCFs agree on only 92.7% of sites by
  (POS, REF, ALT), so the count must run on the biallelic-SNP intersection (249,750
  records); and their obfuscated output was not kept, so their run must be repeated on OUR
  VCF.
- **Why:** the direct comparison this project is for.
- **Prerequisite:** `B-THY-diploid-gate`, `B-DAT-rebuild`, `B-LCK-preregister`.
- **Exit gate — REWRITTEN 2026-10-08 for decisions E, G and H.** For each of the 44
  subjects, a `loo_cohort` at K = 86 with the donor panel and support counts
  recomputed without BOTH of the target's haplotypes, applied to BOTH mechanisms,
  with PanMixer's as-published numbers also reported so the cost of the correction
  is visible. Both emit `new_haplotypes`; the shared evaluators run unchanged.
  **R = 1 release per (subject, operating point)** — the paper's `n` is the subject
  count — plus **3 seeds at ONE declared configuration** as the tilt diagnostic,
  labelled not-for-publication (decision G). The EXACT expectation is printed
  beside every sampled mean, so the small R is a diagnostic and never the evidence.
  Attacks run PER operating point, never pooled, because a frontier of k points is
  k correlated releases and our guarantee is one-release-per-genome.
  ⚠ **Three clauses of the original gate were false or unbuildable** and are
  removed: *R >= 10 seeds per (subject, level)* (superseded by decision G, and it
  demanded ~150x the compute Dylan has agreed to spend); and *allele frequencies*
  and *block assignment* recomputed per subject — **model C has neither**, no
  allele-frequency table and no blocks, so those two could never have been met.
- **Utility axis (decision H):** `utility_loss_panmixer` only. `af_loss` and
  `ld_loss` ARE defined in Arm 1 — the release set is cohort-shaped — but are not
  adopted. That is a scope choice, not an impossibility, and `docs/terms.md` has
  been corrected to say so.
- **Attack axis (decision I):** "self" is scored from the target's own
  assembly-derived path, defined for all 44. PanMixer's panel-row semantics —
  defined for only the 39 members present in the 30x panel — run as a labelled
  parity column.
- **Must include:** a `target_fidelity` axis, or our `tau` = 0 point looks like a
  free lunch — a pure-prior draw is a real cohort haplotype and scores near-perfectly
  on every cohort-level utility metric while carrying zero target information.

### [B-ATK-lrt-arm2] The Arm 2 attacker — RE-FILED 2026-10-10

⚠ **IT WAS BLOCKED ON A DEPENDENCY IT DOES NOT HAVE.** Filed as
`[A-ATK-matched-pair]` with prerequisite "`A-EVL-headtohead` plumbing" — the Arm 1
`new_haplotypes` emitter — which the Arm 2 likelihood-ratio attacker does not
consume at all. Arm 2 has no leave-one-out, no K = 86 panel and no PanMixer
evaluator in the loop, so a whole paper section was waiting on the other arm for no
reason. **Prerequisite: none.**

✅ **THE ATTACKER DOES NOT HAVE TO BE FOUND IN THE LITERATURE.** The proposal
specifies no attacker — seven pages of theory, no evaluation section, the word
"adversar" does not appear — but its algebra hands over the provably optimal one.
With `Q^tau_p(y) = R_D(y) exp(eta_tau u(p,y)) / Z_p`, the Neyman-Pearson statistic
between two input paths is `eta_tau [u(p,y) - u(q,y)] - [log Z_p - log Z_q]`, with
`R_D` cancelling. Its advantage EQUALS `TV(Q_p, Q_q)` — the quantity Theorem 1
bounds — so no attack can beat it. Ingredients all exist: `calibration.eta_tau`,
`utility.u_of_output`, `ffbs.log_z_p`, `toy.sharp_pair`.

⚠⚠ **THE BINDING RISK IS POWER, NOT STRENGTH.** Even the optimal attacker returns a
null at any feasible `n` (R5: signal-to-noise about 1:900). A null without a power
statement is what standing constraint 7 forbids calling evidence; a null WITH one is
the privacy evidence. So every attacker needs a POSITIVE CONTROL — a stated
condition under which it must demonstrably SUCCEED — and the pre-registration must
pin the minimum detectable advantage the experiment is powered for.

⚠⚠ **THE STANDARD GENOMIC ATTACKS ARE MIS-SPECIFIED AGAINST THIS RELEASE, AND ONE
INVERTS.** The release contains NONE of the target's own alleles: every released
allele is some cohort donor's, and the only target-dependent quantity in the object
is the tilt. So a database-search linkage attack has no correct answer in Arm 2 (the
target is outside the database by construction, which makes 0% the experimental
design rather than privacy), and a membership attack INVERTS — the release is a
donor mosaic, so a high match fraction identifies DONORS, not the target.
**An Arm 1 attacker carried into Arm 2 unexamined would produce a number that looks
like a privacy result.** Homer/Sankararaman assume a pool the target contributed to;
the Beacon LRT assumes a presence oracle; we publish neither.

### [A-ATK-matched-pair] The discriminating experiment — SUPERSEDED by the above
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
