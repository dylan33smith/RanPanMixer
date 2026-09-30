# memory.md — the linear ledger

**DO NOT read this file whole.** `grep` it.

```
grep -n "^## 20" docs/memory.md          # list entries
grep -n "INCORRECT\|CORRECTION" docs/memory.md   # what we got wrong
```

`docs/plan.md` exists so a new session does not have to come here.

## Rules

1. **Append only.** Never delete, never overwrite. New entries go directly
   BELOW the APPEND marker, so everything after the marker reads
   **newest-first**. (The rule previously read "newest at the bottom", which
   never matched the file; the file is the source of truth for history.)
2. **In-place correction.** Prepend `[INCORRECT] - ` to the wrong line, preserving
   its text verbatim; insert `[CORRECTION - YYYY-MM-DD]: ` directly below it. The
   wrong version stays — it is what makes the reasoning legible later.
3. **Provenance or it didn't happen.**
4. Metric names come from `docs/terms.md`.

---

## UNEDITED ARCHIVE — prior history

There is none. The repository's first commit (`30b3a1a`, 2026-08-29) contained two
PDFs and nothing else: the proposal and the documentation framework. No prior
notes, logs or results exist to import. This banner records the absence
deliberately, so that a later reader does not go looking for a history that was
never written.

---

## 2026-08-29 — A-FIX-docs-system: adopted the six-file documentation system

**Goal.** Stand up the documentation structure of
`archive_docs/The_Six_File_Lab_Record.pdf` on a greenfield repository, with the
content derived from `paper/Private_Genome_Path_Release.pdf`.

**Method.** Read both source PDFs. Wrote `docs/terms.md` and `docs/data.md` first
(the framework's advice: they carry the anti-drift load and cannot be generated
from existing material). Moved the two PDFs out of `docs/` so that directory holds
exactly the six files — the proposal to `paper/`, the framework to `archive_docs/`,
both with `git mv` so history follows. Wrote the verifier before any code exists.

**Result.** Six files plus `tests/test_docs_contract.py`. `CLAUDE.md` is 137 lines
against a 150-line budget and carries no findings. Every path in `docs/data.md`
except the documentation set is `PLANNED`, which is the honest state: this project
has no code, no graph, no cohort and no target genome.

**Deviations from the framework, and why.**
- The framework's Stage A / Stage B measurement split (selection vs
  characterisation) does not apply to a release mechanism. It is replaced in
  `docs/terms.md` by the denominator trap this project actually has: utility over
  ALL blocks versus over SWITCHED blocks only. Same failure mode, different axis.
- The reference verifier asserts that `memory.md` contains at least one
  `[INCORRECT]` line. On a new project that fails for the wrong reason, so that
  check skips until the first retraction exists and fails loudly thereafter.
- Two checks in the reference implementation had no analogue here (a pinned
  regression baseline, and a contested-values registry). They are replaced by a
  check that recomputes the `tau` calibration table from its own definitions —
  the closest equivalent for a project whose primary artifact is a proof.

**Source:** `paper/Private_Genome_Path_Release.pdf`,
`archive_docs/The_Six_File_Lab_Record.pdf`.

### Decision — `viterbi_release` retired before use

MAP/Viterbi decoding of the tilted HMM is not a draw from `release_distribution`
and carries no privacy guarantee, while producing output indistinguishable in
appearance from a valid release. It is the natural thing to reach for in an HMM
codebase. Retired in `docs/terms.md` so that the name resolves to its refutation
rather than to nothing. **This is a decision from the proposal's text, not a
measurement.**

### Finding — two typesetting defects in the proposal PDF

Verified by rendering page 3 at 300 dpi, not by text extraction (which mangles
superscripts and would not be evidence).

- **Theorem 1**, the pointwise likelihood-ratio bound, is typeset as
  `e^{-eps} <= Q_p(y)/Q_q(y) ^{eps}` — the upper half of the inequality has lost
  its `<= e`, leaving `eps` floating as a superscript beside the fraction. The
  intended statement is `e^{-eps} <= Q_p(y)/Q_q(y) <= e^{eps}`.
- **The proof's first line** has the same defect: `1 <= exp(eta*u(p,y))^{eta}`
  where `1 <= exp(eta*u(p,y)) <= e^{eta}` is intended.

Both are the same LaTeX slip losing `\le e` before a superscript. **The
mathematics is unaffected:** the chain still yields `TV <= tanh(eta_tau) = tau`,
and `epsilon_tau = 2*arctanh(tau)` is correct as printed elsewhere.

⚠ **Why this is recorded rather than ignored.** A one-sided bound transcribed
literally into an assertion would pass on every input, because the missing half is
the half that can fail. `docs/terms.md` states both halves explicitly and the
`log_z_p` bracket `[0, eta_tau]` is written as a two-sided assertion for the same
reason. Fix the manuscript when `paper/` is next edited.

**What this does NOT establish:** nothing about the correctness of the mechanism.
It is a defect in the document, found by reading it, before any code was written.

---

<!-- APPEND NEW ENTRIES BELOW THIS LINE -->

## 2026-09-29 (later) — A-IMP-preprocess: independence from PanMixer, and a real join bug

**Goal.** Dylan's direction: the tool must be COMPLETELY INDEPENDENT of PanMixer, which is
here only as a published example of a problem in the same space. Plus re-run chr21 through
the committed code path rather than by hand, and recheck everything.

### Independence — done and verified

- **Our own conda env** `ranpanmixer` (`environment.yml`, plus `requirements.txt` for a
  venv). Nothing from their `environment.yaml`. `src/` imports only stdlib, numpy and
  pyyaml; scipy is present for the sampler. **pytest is in our env**, so both suites run
  normally (it is absent from `panmixer`, which is why `test_preprocess.py` is self-running).
- **The dev fixture is rebuilt from OUR VCF.** It had been sliced from PanMixer's.
- **`vg` relocated** to `data/tools/vg`, pinned by sha256 in the manifest, overridable via
  `$RANPANMIXER_VG`. Verified by re-running chrY deconstruct through it: 34,014 records.
- **The 1000G panel staged** to `data/raw/`, recorded as evaluation-only.
- `grep` finds **no** reference to `external/PanMixer` in `src/`, `configs/` or `tests/`, and
  tracing every provenance sidecar shows PanMixer's chr21 VCF sha256 `5b880a3a…` appears in
  **no** artifact's inputs.

### chr21 re-run through the committed path

8:33 / 9.475 GB, **340,849 records, 45 samples** — reproducing the by-hand run exactly
(8:22 / 9.477 GB / same counts). The CLI performed its own decompression, confirming the gz
fix. So the artifacts are now provably the output of committed code.

### ⚠ A REAL BUG: the strict join could never match a multi-allelic record

**The 1000G panel is entirely decomposed into biallelic records** — every one of its
1,002,752 chr21 records declares exactly one ALT — while **7.9% of ours (27,021) are
multi-allelic**. A strict (POS, REF, full-ALT-string) join therefore **cannot succeed on any
of them**, for a purely notational reason, and the first implementation silently discarded
target information we actually held.

Recoverable by matching each of OUR alts against its own panel record and mapping back to
OUR allele numbering:

| | records |
|---|---|
| strict match | 258,610 |
| **recoverable via split alts** | **8,377** (all multi-allelic; 8,283 inside `chain_span`) |
| genuinely absent from the panel | 73,862 |

Implemented as `split_alt_join`, on by default. Effects, measured across all ten targets:
sites matched **258,610 -> 266,987**; chain positions with a target call **82.4% -> 85.1%**;
and once the accounting was made consistent, representable non-ref calls **80.0% -> 83.5%**
(a mean of 2,269 calls per target were being miscounted as unrepresentable).

⚠ **The decomposition assumption is not universally safe, and the counter proves it.**
**3 contradictions** across the ten targets (0.04% of recoveries) — a haplotype reading as
carrying two different alts. Examined: chr21 **28,186,675**, REF `C`, our alts
`CA,CAAA,CAA,CAAAA` — a **poly-A repeat** where the panel's biallelic records are NOT
mutually exclusive, and haplotype 0 is `1` at both `C->CA` and `C->CAA`. The code counts it
and **leaves the row unresolved rather than guessing**, and those calls fall back to
`pos_only`. That is the correct behaviour and the contradiction counter is what makes it
visible.

**Second-order bug this exposed:** the recovery initially improved `path.npy` while
`representability.json` still classified those calls as unrepresentable — the two artifacts
disagreed about the same fact. Fixed by deferring the classification until the reconstruction
resolves.

### Three other items closed

- **`support(v) == 0` is excluded from the chain**, not patched. `w_v = 1/(support+1)` was
  the tempting fix and is wrong: at support 0 no donor has an allele, so the sampler has
  nothing to emit and RENORMALISE has an empty pool — the position is unusable, and a finite
  weight it can never earn would just hide that. Removes 0 positions today (the one such
  site is already outside the span), so the guard no longer rests on that luck.
- **`target_reason.npy`** distinguishes `no_record_in_target_callset` from
  `record_present_but_no_call`. Measured: the inert 19% is **53,861 code 1 and ZERO code 2**
  — a phased panel has no missing genotypes inside its own site list. So the unmatched set
  is identical for all ten targets and provably carries no per-person variation, which is
  the property needed to normalise `beta_v` over the panel-matched positions legitimately.
- **A fabricated checksum, caught by the check written to catch it.** A `vg` sha256 was
  written into the manifest with its tail invented from a truncated terminal display. It
  looked plausible. Corrected, and a test now compares every declared tool digest against
  the binary on disk — verified to fail against a corrupted one. A checksum that was never
  computed is the worst violation of the traceability rule, since being checked is its only
  purpose.

**Provenance.** 40 tests pass in the `ranpanmixer` env. Audit workflow `wf_814ac6b2-5f2`
launched four auditors; all four were cut off before returning, but one flagged the
multi-allelic lead that this entry's main finding came from — chased and confirmed by hand.

---


## 2026-09-29 — A-IMP-preprocess: the preprocessing pipeline, built and run end to end

**Goal.** Implement preprocessing for the cohort and one target, verify every step, and
audit the result. First code this project has that is its own rather than PanMixer's.

**Result: it runs, and every artifact cross-checks.** 17 unit tests, 28 consistency checks
on the real output, all passing. Still no mechanism and no RanPanMixer number.

### Our own chr21 build, and how it differs from PanMixer's

`vg deconstruct` on the chr21 PGGB GFA, **8:22 wall / 9.5 GB peak at -t 32** — far cheaper
than the 1-4 h estimated. ⚠ **Sample grouping was the risk and it is settled**: paths are
PanSN `SAMPLE#HAP#CONTIG` with each assembly split across many contigs, and `-H` is
deprecated in vg 1.68, but a chrY smoke test (2:01, 34,014 records) produced 17 clean sample
columns. chr21 produced **45** (44 donors + chm13), as it should.

**Our VCF is NOT the published one, and every derived number must say which it is.**

| quantity | OURS (vg 1.68) | PanMixer's (vg 1.36) |
|---|---|---|
| chr21 records | **340,849** | 340,824 (+25) |
| sample columns before chm13 | 45 | 45 |
| missing matrix cells | **1,189,258** (3.965%) | 1,278,653 (4.26%) |
| bare `.` GT fields | **0** | 237,594 |
| declared alleles | **767,376** | 767,324 |
| T inside `chain_span` | **305,886** | 305,887 |
| dropped below span / above | **34,371** / 592 | 34,345 / 592 |
| sites at support 0 | 1 | 1 |
| max alleles at one record | 90 | 90 |

**⚠ AND THE ALLELE INDICES DIFFER — the most important finding of the run.** "340,849 vs
340,824" understates it badly. On the strict (POS, REF, ALT) key only **94.9%** of records
match: 17,454 are ours-only and 17,429 theirs-only. But comparing on (POS, REF) and then
asking whether the ALT *set* matches shows the cause:

| | records |
|---|---|
| same ALT set, order may differ | **339,358** |
| genuinely different ALT set | 386 |
| (POS, REF) only in ours | 1,105 |
| (POS, REF) only in theirs | 1,080 |

So the 5.1% mismatch is **almost entirely ALT REORDERING**: `T -> TAA,TA,TAAA` in ours
against `T -> TAA,TAAA,TA` in theirs, the same variant with the alleles listed in a
different order. **Allele index 1 in our VCF is not allele index 1 in theirs** at those
sites. Any PanMixer-derived artifact joined positionally to ours — the allele matrix, a
released `new_haplotypes`, the attack database masks — would be silently wrong there, with
no error and no shape mismatch. This is precisely the hazard `site_axis_digest` exists to
catch, and it is now a measured fact rather than a worry: the two VCFs' axes genuinely are
different objects. **Never mix artifacts across the two builds.**

⚠ **vg 1.68 calls 89,395 more genotypes than 1.36 on the same graph** — and emits **zero**
bare-`.` GT fields where 1.36 emitted 237,594, writing `.|.` instead. Both counts were
verified independently with bcftools against both VCFs. Two consequences: the "haploid
genotype" cause of a `-1` does not merely happen to be always `.` on our data, it **does not
occur at all**; and the chrY divergence measured on 2026-09-28 (92.7% of tuples matching)
was atypical — chr21 differs by 25 records in 340,824, i.e. 0.007%.

### Verified, not assumed

- **`INFO/AN` == the non-missing count, exactly, at all 340,849 sites.** The prior
  measurement was on PanMixer's VCF; it holds on ours too. That makes the cross-check a
  free row-alignment canary, and it is not vacuous — the tiny fixture carries one
  deliberately wrong `AN` and the check fires on exactly that record.
- **Reason codes hand-verified.** Every count on the fixture matches hand calculation,
  including the precedence collision at pos 720 where S3 is both conflict-named and
  would-be-inherited: conflict wins, as pinned.
- **Code 4 is not "unexplained".** 571,182 uncategorised cells, of which **99.3% sit in a
  run of >= 10 consecutive missing** — they are assembly gaps, which deliberately get no
  per-cell code and live in `missing_runs.tsv` as intervals instead. The bucket is behaving
  exactly as designed.
- **The coordinate cut is doing its job.** 8,327 records sit in [10,326,676, 12,968,320] —
  inside the genetic map but outside `chain_span` — and every one carries cM **0.584144
  exactly**. That is the centromeric plateau the cut exists to exclude, observed directly.
- **Representability on all 10 targets** (in-region, so a partial fixture cannot report its
  own window as a failure): mean **80.0% exact / 10.2% invisible / 21.5% of information
  invisible**, against the 2026-09-27 figures of 80.4% / 9.9% / 20.5% measured on six
  samples and the published VCF. Independent reproduction.

**⚠ NEW FINDING — representability is ancestry-dependent.** The two AFR targets are the
worst on both axes: NA19350 and NA20294 carry the most non-reference calls (78,595 and
79,060 against ~60,000 for EUR/EAS) and have the lowest exact-match rate (**78.9%** against
81.1% for HG01468, AMR). The cohort represents a EUR target better than an AFR one. This is
exactly what the deliberately multi-continental target set was chosen to expose, and it is
a utility-ceiling effect that a single-population test set would have hidden entirely.

### Two real bugs, both found by running on real data

1. **`_csv.Error: field larger than field limit`** — `load_sites` could not read its own
   output. The hypervariable record's ALT string is **23.2 MB** and Python's csv module
   refuses fields over 128 KB. Fixed with `csv.field_size_limit`. Found in seconds because
   the dev slice was chosen to contain that record.
2. **`sites.tsv` had no provenance sidecar** — the one artifact that defines what every row
   means was the only one written outside the `store.save()` chokepoint, because
   `write_sites_tsv` streams it to compute the digest in the same pass. Added
   `store.attach_prov()` as the narrow, deliberate exception.

Both now have regression tests. A third apparent bug was **my own audit check** being
wrong, not the code: "all cM outside the span are NaN" is false by design, because the
centromeric plateau is inside the map and outside the span.

**One fixture error of my own**, caught by the canary it was meant to test: record 522 had
`AN=4` where the truth is 5. Corrected, so only the deliberate error at 1300 remains.

### Provenance

`data/cohort/chr21/` 240 MB, 12 artifacts each with a sidecar, one axis digest
`251b554f...` across all of them and all ten target paths. `data/targets/<id>/chr21/` holds
`path.npy` + `representability.json` per target. Dev fixture at `data/dev/` is the identical
code path via `--root`. `pyyaml` was added to the `panmixer` conda env (absent upstream).

---


## 2026-09-28 — A-FIX-panmixer-framing + A-DAT-genetic-map-parm + A-DAT-input-set + A-DAT-evaluation-fit

**Goal.** Start the build. Create a clean branch, decide the chain architecture, settle the
preprocessing input and output sets, and correct the documentation where the PanMixer
comparison had been forced onto a problem it does not describe.

**Method.** Two adversarially-verified workflows (10 and 12 agents) plus hand verification.
Every number below is measured in this project on chr21 unless labelled otherwise.

### Decisions (Dylan)

1. **Model C is the chain architecture** — one chain over the chromosome, units = individual
   VARIANTS. Model B (units = LD blocks) is the retained fallback; **model A is ruled out**.
   This answers open question 8, raised 2026-09-22 and deferred to the PI.
2. **The scope decision is narrowed.** *"Maximise reuse of PanMixer, change exactly ONE
   thing"* imported PanMixer's threat model — target-as-cohort-member — into a design whose
   target is EXTERNAL to everything. Dylan: *"Boundaries aren't a weak function of the
   target's genotypes because the target is completely outside of any other dataset."*
3. **Write fresh, do not port.** Copy functions from PanMixer where they are right, but the
   file structure and pipeline are ours, so nothing is inherited without a decision.
4. **The pipeline starts from our own inputs and produces all outputs.** `vg deconstruct` is
   pipeline step 1. Accepted knowingly: our VCF will not match PanMixer's.
5. **`missing_policy` = RENORMALISE.** The cost — penalising poorly-assembled donors for a
   defect that is not theirs — was accepted deliberately, to be noted in the writeup rather
   than modelled around at this stage.
6. **Target mapping is a separate stage and a separate command**, so a new target does not
   re-run the cohort stage.
7. **`chain_span` for chr21 starts at 12,968,320**, flagged for revisit.
8. **The k-mer Jaccard is deferred** to a post-build evaluation/reporting tool.

### The framing audit

91 findings filed across `CLAUDE.md`, the five docs and the three primer sections;
adversarial verification **upheld 76, rejected 15, and found 51 further instances** the
auditors missed. The rejections were almost all the predicted failure mode — a correct
description of PanMixer read as a claim about us.

Highest-leverage: **`CLAUDE.md` constraint 1** ended with an unscoped *"Recompute
leave-one-out"* inside the hard rules, in the file auto-loaded every session — so the error
regenerated every session. Also `terms.md`'s `loo_cohort`, which defined OUR cohort as
`HPRC \ {target}` at Status PRIMARY; `plan.md`'s exit gate 5, which made PanMixer benchmark
compatibility a blocking gate of the sampler; and `03-panmixer.md:2363`, which asserted
block boundaries are a function of the target's genotypes — flatly false of our setting and
internally contradicted by line 1220 of the same file.

### The genetic map, and the p arm

**The defect found is worse than the one we were looking for.** `np.interp` clamps by
default, so all 26,018 variants below the map's start (10,326,676) take its first value —
and the map's own first interval (10,326,676 -> 12,968,320, spanning the centromere) carries
**dcM = 0.000000**. Together: **34,345 consecutive positions with identical cM**, so
`Delta_x = 0`, `P(switch) = 0`, and the sampler emits **one donor's real haplotype verbatim
across 7.25 Mb** — departure (b) at ~200x the scale of the worst single block. A privacy
failure, not a utility one, and silent.

**No map fixes it.** deCODE 2019 (native hg38, highest resolution) starts at **14,143,522**,
worse — it would leave 47,805 variants (14.03%) unmapped against 26,018 (7.63%). The 1000G
hg38 map has a single value-0 interval over chr21:0–10,326,676. Bherer 2017 is GRCh37 and
starts later still. pyrho hg38 (Spence & Song 2019) DOES cover from **5,088,754 bp** and
would give every variant a distinct cM — but it spends **9.98 cM (ACB) / 9.15 (CEU) across a
single 3.4 Mb interval defined by only two markers**, attributing 23–26% of chr21's genetic
length to a region with zero observed crossovers. Interpolation across a void.

**Why the region is like that — the p arm.** Every chromosome has a centromere dividing it
into a short arm (**p**, French *petit*) and a long arm (**q**). Five human chromosomes —
13, 14, 15, 21, 22 — are **acrocentric**: the centromere sits near one end, so the p arm is
small and is almost entirely repetitive DNA plus tandem arrays of **rDNA** genes. Those
arrays are near-identical BETWEEN the five, so assemblers cannot place them uniquely and
GRCh38 fills them with N. chr21's p arm is 0–~12 Mb.

Five independent legs support excluding it: GRCh38 chr21 5,010,000–10,814,560 is placeholder
model sequence (25 contigs, fabricated N-gaps, 22 of exactly 50 kb); 1000G's strict
accessibility mask rules **26,014 of 26,018** (99.98%) inaccessible; our cohort matrix is
**36.76%** missing there against 6.57% elsewhere; deCODE's pedigree map assigns **0.0 cM/Mb
to every 1 Mb bin from 0–13 Mb**; and a four-generation CEPH pedigree recorded *"Not a
single allelic recombination was observed on the p-arm"* across 107 transmissions.

⚠ **Two smaller zero-distance defects remain and the cut does not touch them**: 592 variants
above the map end give a 593-position zero-distance run at the telomere (UNDECIDED), and 618
interior map intervals with dcM = 0 hold ~1,650 kept variants (correct to leave — real map
plateaus). ⚠ **What the cut costs is unmeasured** — the share of a target's `-log f`
discarded has not been computed.

### The input set

- **PanGenie drops out entirely.** Exactly two live consumers, both deleted by our design;
  no evaluator reads it. Saves a 5.06 GB download.
- **The external panel drops out of the mechanism** under model C — no blocks, no PLINK.
  It survives only as the attack database.
- **`INFO/AN` is bit-identical to `support(v)`.** `np.array_equal` is True across all
  340,824 chr21 records. Support needs no computation, and the comparison doubles as a free
  row-alignment check.
- **Do NOT compute LD blocks from the 44-sample cohort** (for the B fallback): 7,858 blocks
  covering 50.4% of variants against 14,137 / 74.2% from the panel; largest block 2,562
  against 618. Only 44% of SNPs clear MAF 0.05 in 88 haplotypes.
- **The haploid-genotype cause of `-1` does not exist on chr21.** All 237,594 non-piped GT
  fields have exactly one distinct value, `.`. Four causes, not five.
- **The graph is published only as GFA** for PGGB freeze1 — no GBZ, no snarls. chr21 =
  434,772,428 B; whole genome 15.64 GB compressed / 86.5 GB uncompressed.
- **`vg` 1.68 does not reproduce `vg` 1.36.** On chrY it recovers 99.55% of the published
  VCF's sites by POS but only **92.7% by (POS, REF, ALT)**, concentrated at multi-allelic
  sites. Also: the published VCF carries 45 samples, our local chr21 carries 44 — `chm13`
  removed by PanMixer's own step — so a rebuild does not drop in as a replacement.
- **The graph buys little the VCF lacks**: 98.37% of consecutive top-level snarls share a
  boundary node, so donor switching is provably a valid graph walk from `AT` alone; all
  340,824 alleles are explicit sequences. The one irreplaceable use is mapping a target's
  reads — and HPRC marks short-read mapping **"untested"** for PGGB.
- **The graph does carry the p arm**: the CHM13-backbone VCF has **185,414 records** in
  `chm13#chr21:1-10 Mb` against 22,746 for GRCh38. But that recovers sequence, not genetic
  distance, and missingness there is 43.9%.

### The evaluation

**Nothing in PanMixer's suite tests the claim we make.** Read every evaluator:
- `af_loss` and `ld_loss` reconstruct an edited cohort (`new_pangenome[subject_index] =
  released_array`). We never edit the cohort, so they are not defined for our release
  object. ⚠ A `--dont_replace` flag exists that skips the substitution, so this is
  conditional rather than structural — but the decision stands regardless: Dylan,
  *"we don't use af_loss for our project at all."*
- **The gap score is the one that transfers**, after three changes: "self" becomes the
  external target's true path; the attack DB excludes the target; the AF source is disjoint
  from the cohort. ⚠ Not reusable verbatim — `score_genotypes_remove_missing` writes
  `test_hist.png` to the working directory on every call, 11x per row.
- **Their shape assertion is vacuous.** Both sides equal `site_mask.sum()` by construction,
  so a differently-ordered release is silently mis-scored rather than erroring. This is the
  row-order hazard, located.
- **Two of their metrics are broken**: `ld_loss`'s window is **5 bp, not 5 kb** (7,629 pairs
  on chr21 instead of ~10.6M), and the haplotype gap score sets `ref_af = alt_af`, so
  reference matches are weighted by the alt frequency.
- The test we need is **two-target indistinguishability** — new entry in `docs/terms.md`.

### Two implementation facts worth not rediscovering

- **The budget is invariant in T.** `u = sum_v beta_v phi_v` with `sum beta = 1`, so the
  total tilt between a perfectly-matching and a non-matching path is `exp(eta_tau)` however
  the path is partitioned. At tau = 0.5: per-position tilt 1.0000055 at T = 100,757 against
  1.0000016 at T = 340,824, total 1.7321 either way. Finer T does NOT spread the budget
  thinner. ⚠ What it DOES cost is float32 headroom: the per-position tilt under C is only
  **13.5x** the float32 resolution against 45.7x under blocks, so float64 stored alphas are
  load-bearing. And the thin-budget problem is a property of **tau**, not of T.
- **Never materialise `A`.** Li-Stephens `A = (1-s_t)*I + (s_t/K)*J` is rank-one plus
  diagonal, so the forward step is O(K) from one scalar sum and the backward draw likewise.
  A materialised `A` per position would be 340,824 x 88 x 88 x 8 B = **~21 TB**.

**Provenance.** Workflows `wf_1dd6bfa5-5a8` (inputs, 10 agents) and `wf_d0195c3c-52d`
(framing audit, 12 agents), both adversarially verified; `wf_785b6f66-755` applied the primer
corrections. Doc corrections in commit `1f8ea3f`. No mechanism has been run and no
RanPanMixer number exists.

---


## 2026-09-27 — A-THY-chain-architecture: three models, measured, and a toy that shows the tilt vanishing

**Why.** Dylan asked whether "the HMM is per block and the chain never crosses a block
boundary" is the right framing for OUR mechanism, explicitly from first principles rather
than from the docs. It is not: that sentence is a FINDING ABOUT PANMIXER'S CODE which
leaked into our spec as though it were a design principle. Nothing in the proposal or in
Theorem 1 requires it.

**First-principles position.** `baseline_model` is a model of what a plausible haplotype
looks like: a mosaic whose segment boundaries fall where recombination happened, at a rate
set by genetic position in centiMorgans. **LD blocks appear nowhere in that.** They are a
statistical summary of where recombination has historically been rare — downstream of the
same process the transition model already encodes. Using blocks INSTEAD of the map
substitutes a coarse proxy for the real quantity; using both double-counts.

**The three options, measured on chr21** (acquired GRCh38 map, Ne = 10,000, r = 1.26,
n = 86; analytic segment n/(Ne*r) = 0.0068 cM ~ 4.5 kb):

| model | chain | E[switches] | mean segment | alphas |
|---|---|---|---|---|
| A: per-block chains, no inter-block (pinned) | 100,757 separate | 104,379 | 0.39 kb | 69 MB |
| B: one chain, units = blocks | 100,757 | 3,532 | 11.60 kb | 69 MB |
| C: one chain, units = variants | 340,824 | 7,155 | 5.73 kb | 234 MB |

- Columns: one map, one formula, our real positions.
- Rows: A is the pinned spec; B and C are candidates.
- Synthesis: **A is 11x too scrambled and is the worst of the three.** Its fresh-uniform
  reset per block is algebraically the `d` -> infinity limit of the switch formula. Measured
  inter-block switch mass is median **0.0022** against the **0.9884** independence implies,
  and 98.5% of block pairs are below half of it. B discards **51%** of the recombination
  (the within-block half); C is nearest the analytic segment length.

**Supporting measurements.** 89,087 of 100,757 block-dictionary entries (88.4%) are orphan
variants in NO called plink block, so block boundaries are mostly bookkeeping, not
recombination sites. Map resolution is 407 bp median against 60 bp variant spacing, so the
map already resolves finer than blocks do. The largest block is 618 variants over 39.5 kb
which C recombines through 73% of the time and B emits as one donor verbatim.

**Numerical findings.** An unnormalised forward pass underflows float64 after about **183
positions**; per-column normalisation (`log_z_p` = sum of log `c_t`) removes it entirely and
is preferred over log-space (cheaper, and `alpha_t` is already the distribution backward
sampling needs). **Storage must be float64**: float32 resolution 1.2e-07 against a
per-position tilt of 5.5e-06 is only 46x of headroom.

**The toy (4 donors, 7 variants, 3 blocks), model B.** Forward pass reproduces `Z_p` exactly
against enumeration over all 64 paths (1.330113 both ways); `log_z_p` = 0.285264 inside
[0, 0.5493]; E[u] moves 0.5000 -> 0.5386 at tau = 0.5.

**The tilt vanishing, shown directly.** Tiling the same 3-block pattern to longer chains,
holding tau = 0.5:

| T blocks | max psi | E[u] gain |
|---|---|---|
| 3 | 1.2654 | +3.86e-02 |
| 300 | 1.0024 | +3.72e-04 |
| 30,000 | 1.0000235 | +3.73e-06 |
| 100,755 | 1.0000070 | +1.11e-06 |

Exactly 1/T, and **identical under all three architectures** — the chain choice does not
touch the budget. The same toy at T = 3 gives E[u] gains of +0.0402 (A), +0.0386 (B),
+0.0371 (C): A buys marginally the most utility precisely BY discarding the linkage
constraint, which is also why its output is least haplotype-like.

**Status: OPEN, deferred to Dylan's PI**, written up as a decision in `docs/plan.md`. If C
is chosen the factorisation argument becomes unnecessary and the singleton special case
disappears, making C simpler than what is pinned rather than harder.

## 2026-09-27 — A-DAT-genetic-map: map acquired, and exit gate 3 is mis-specified

**Why.** The GRCh38 genetic map was the last unacquired input for `A-IMP-cohort-hmm`.
Earlier this session I reported the SHAPEIT4 b38 maps as returning 302 — that was wrong,
a redirect I failed to follow. Both candidates return 200 and were downloaded.

**Result 1 — the two candidate maps are the same map.** SHAPEIT4 and Beagle both yield
chr21 with **44,618 positions, identical positions, identical cM, maximum absolute
difference 0**. The choice is format and provenance only, not content. Two independent
distributions agreeing exactly is evidence the map is not a single-source artifact. Both
registered in `docs/data.md`; the SHAPEIT4 copy is the working file, being pinnable by
commit SHA rather than served from a lab webserver.

**Result 2 — 7.63% of our chain positions sit outside the map.** Map span is
10,326,676 - 46,680,243 bp; our chr21 variants span 5,713,651 - 46,699,788. So **26,018
variants (7.63%) fall below the map start** and 592 (0.17%) above. That region is the
acrocentric p-arm and heterochromatin, where no recombination map exists for a real
reason. Interpolation clamps, giving flat cM, so **296 multi-variant blocks there have
P(mosaic) = 0 exactly** and can only ever be one donor copied verbatim. This is the hazard
flagged from Pickrell's README, confirmed on our data; it is a property of the region, not
of the map we picked.

**Result 3 — the transitions work, and exit gate 3 as written does not test them.**
Map interpolated onto all 340,824 chr21 positions, published formula (Ne = 10,000,
r = 1.26, d = `Delta_x` * Ne * r, total switch mass = (n-1)(1-exp(-d/n))/n, n = 86):

| quantity | value |
|---|---|
| median `Delta_x` between adjacent chain positions | 1.811e-05 cM (about 12 bp) |
| per-step switch mass, median | 0.0026 |
| per-step switch mass, p75 / p90 | 0.0098 / 0.0304 |
| steps with switch mass in [0.1, 1] | 3.0% |
| expected donor switches across chr21 | **3,622.6** |
| P(block is a mosaic), median over 11,670 multi-variant blocks | 0.0866 |
| P(block is a mosaic), median over the 3,523 blocks with >= 20 variants | 0.2747 |
| PanMixer shipped, same measurement | 5.8e-12; 300/300 blocks single-donor |

- Columns: one map, one formula, our real chain positions.
- Rows: per-step figures are between ADJACENT variants; block figures span a whole block.
- Synthesis: median per-step switch mass is **4.5e8 times** PanMixer's, and the model gives
  mosaic segments of a few kb, consistent with the analytic n/(Ne*r) = 86/12,600 =
  0.0068 cM. **But exit gate 3 asks for "switch mass of order 0.1 to 1 at realistic
  separations", and that target is biologically wrong at this density.** Adjacent chr21
  variants are about 12 bp apart; switch mass near 0.1 there would mean roughly one
  ancestral recombination per 12 bp. The measured 0.0026 is the CORRECT value.

**Exit gate 3's second clause is mis-specified too, and more subtly.** It requires that
"a sampled block is demonstrably NOT a single donor copied verbatim". With the correct map
the median multi-variant block still has only **8.7%** chance of containing a switch, so
about 91% of blocks WILL be one donor copied verbatim. That is NOT PanMixer's failure
repeated: theirs was 100% from an arithmetic bug at 5.8e-12; ours is ~91% because **an LD
block is by definition a stretch with little recombination**, and v1's chain scope is one
block. The mosaic emerges ACROSS blocks — precisely the scope v1 defers. Inside v1 the tilt
acts on donor SELECTION per block, not on within-block recombination.

**OPEN — needs restating before implementation (Dylan).** Exit gate 3 should be rewritten
against quantities meaningful at block scope. Candidates: expected donor switches per
chromosome (3,622.6); median P(mosaic) among blocks with >= 20 variants (0.2747); or
agreement between measured segment length and the analytic n/(Ne*r). Also open: whether to
extrapolate cM below 10.33 Mb or accept 296 permanently single-donor blocks.

## 2026-09-27 — A-DAT-representability-split: 9.9% invisible, carrying 20.5% of the information

**Why.** The earlier entry today put the invisible share at 19.6%. Its median allele
frequency came out at 0.28 — implausible, since a variant that common should appear in 88
haplotypes. Re-measured by splitting on whether the graph has ANY record at the position.

**Method.** Same three 1000G samples outside the HPRC 44 (HG00096, HG00097, HG00099), chr21.
Three categories per target non-reference call: exact (POS, REF, ALT) match into the graph;
POS present but no exact match; POS absent entirely. Allele frequency is leave-one-out over
the 3,202-sample panel with the target's own two haplotypes removed, singleton-floored at
0.5 counts. Information is `sum of -log f`, which is the published linkage score of Eq. 5.

| category | share of sites | median AF | MAF < 1% | share of information | ratio |
|---|---|---|---|---|---|
| exact match — measured | 80.2% | 0.47 | 0.3% | 68.3% | 0.85x |
| POS only, alt differs — measured | 9.9% | 0.35 | 1.1% | 11.1% | 1.13x |
| **no record — INVISIBLE** | **9.9%** | **0.20** | **16.9%** | **20.5%** | **2.07x** |

- Columns: pooled over three targets; per-sample spread was under 2 points in every cell.
- Rows: only the third is invisible. The second HAS a chain position, because a record
  exists at that POS — the mismatch is a REF/ALT spelling difference, typically a
  multi-allelic graph context, and `phi_t` scores it as an ordinary mismatch.
- Synthesis: **one tenth of a real target's variants have no record at all, and they carry
  one fifth of its identifying information** — a 2.07x over-representation, with rare
  variants enriched about 56-fold (16.9% under MAF 1%, against 0.3% among exact matches).

**The framing that matters (Dylan).** This loss is NOT ON THE `tau` DIAL. It is fixed by the
graph's coverage and happens before `psi_t` is evaluated, so no `tau` recovers it: at
`tau` -> 1, with no privacy guarantee at all, the release still cannot express those alleles.
**The utility ceiling is therefore about 79.5% of the target's information, not 100%**, and
every target_fidelity number must be read against it. The same fact read from the privacy
side: 20.5% of the linkage score is annihilated unconditionally, for every target, for free.
Bad for a clinician, good for a privacy claim, and controlled by neither of us.

**Caveat for reporting.** The middle category is measured but unfair to the mechanism: the
target's exact allele is unavailable for a join reason, so `phi_t` records a mismatch the
sampler could not have avoided. Separate it when reporting utility, or the mechanism is
charged for a representation artifact.

## 2026-09-27 — A-DAT-representability: a fifth of a real target never reaches the mechanism

**Raised by Dylan.** A target variant absent from both the graph and every external callset
never enters any VCF. The release is a path through `G`, so that variant cannot be emitted —
but it also has no record, hence no chain position, hence no `phi_t` term and no `beta_t`
weight. It is removed AND unmeasured. Is anything accounting for it?

**Answer: no, and the magnitude is large.** Measured on six 1000 Genomes samples that are
not among the HPRC 44, chr21, using the strict (POS, REF, ALT) mapping:

| metric | per target |
|---|---|
| non-reference calls | 60,561 |
| with an exact record in the graph | 48,678 (80.4%) |
[INCORRECT] - | **invisible: no record, no chain position** | **11,883 (19.6%)** |
[CORRECTION - 2026-09-27]: 19.6% conflates two categories. Only about half of it — **9.9%** —
has no record at all. The other **9.9%** has a record at that POS with a different REF/ALT
spelling, so a chain position DOES exist and `phi_t` does measure it. See the
A-DAT-representability-split entry above. The 80.4% exact-match figure is unchanged.

Per-sample spread was tight: HG00096 80.3%, HG00097 80.1%, HG00099 80.3%, HG00100 80.6%,
HG00101 80.5%, HG00102 80.5%.

- Columns: one typical external target.
- Rows: "invisible" means no exact-match record exists, so no chain position is created.
- Synthesis (⚠ the invisible share is corrected to 9.9% above): **both bounds move the wrong way.** 80.4% is an UPPER bound on retainable
  fidelity, because having a record is necessary but not sufficient — the allele must also
  be carried by a donor, and 8,823 declared alleles have zero cohort support (2026-09-27
  entry above). And the invisible share — 9.9% as corrected, not the 19.6% written here — is
  itself a LOWER bound, because these are 1000G samples
  whose variants are by construction already known; a genuinely novel variant is invisible
  to this measurement too. The measurement is also SNV-biased, since 1000G short-read calls
  under-ascertain exactly the structural variation the graph is richest in.

**The consequence for our utility claim.** `u(p, y) = sum beta_t phi_t` is defined over
chain positions, and chain positions come from VCF rows. So **u can equal 1.0 while the
release retains only about 80% of the target's variation.** Every downstream instrument
inherits the blind spot: AF loss, LD loss and read mapping all operate on the graph's site
set and none can see a variant with no row.

**Not inherited from PanMixer.** Their target is a cohort member, so their haplotypes ARE
the VCF rows by construction and nothing they carry is inexpressible. This gap belongs to
the external-target setting and has no upstream precedent.

**Privacy is untouched.** The invisible variants cannot be emitted, so they are maximally
protected; and the mechanism's behaviour does not branch on whether the target had one, so
nothing leaks. This is a measurement and utility defect, not a privacy defect.

**Fix, and the shape it must take.** Dylan proposed adding a record at those locations.
That cannot be done inside the mechanism: a target-specific record makes `output_support`
target-dependent, which is the §4.4 Alice/Bob counterexample — disjoint supports, TV = 1,
no guarantee. The safe inversion is to record at MAPPING time as an audit artifact that
never touches `output_support`, `baseline_model`, `beta_t`, `phi_t` or the tilt — the same
pattern as the REASON array. Reporting rule now in `target_fidelity`: the quantity has a
ceiling below 1 independent of `tau`, and must never be quoted without it.

## 2026-09-27 — A-DAT-cohort-support: alleles v1 cannot emit, and what v1 still needs externally

**Question (Dylan).** If v1 routes every block through the tilted sampler, what external
data is still required? And since the target is external and never enters the cohort,
how are target-unique variants handled?

**Result 1 — external dependencies shrink to two for the mechanism.**

| dependency | role in v1 | status |
|---|---|---|
| 1000 Genomes panel | `plink --blocks` boundaries only | present (30x GRCh38) |
| GRCh38 genetic map | `Delta_x` in cM for transition constants | **NOT ACQUIRED — blocks exit gate 3** |
| 1000 Genomes panel | evaluation: attack DB, AF/LD site sets, Beagle reference | present |
| PanGenie callset | **nothing in v1**; needed only to run the RESTRICT `missing_policy` baseline | optional |

Dropping anchors removes PanGenie from the mechanism's critical path entirely. Note also
that block boundaries do not affect correctness — the per-block factorisation gives the
same `tau` for any segmentation with `sum beta_b` = 1 — so plink blocks are inherited for
HEAD-TO-HEAD COMPARABILITY with PanMixer, not out of necessity.

**Result 2 — alleles declared in the VCF that no donor carries.** Measured from
`external/PanMixer/starting_data/chr21/pangenome.npy` (44 x 340,824 x 2) against
`num_alleles.npy`, counting support over all 88 haplotypes:

| quantity | count | denominator | share |
|---|---|---|---|
| declared allele slots | 767,324 | — | — |
| declared alleles with ZERO cohort support | 8,823 | 767,324 slots | 1.15% |
| sites with >= 1 unsupported declared allele | 6,809 | 340,824 sites | 2.00% |
| sites where the unsupported allele is REF | 4,885 | 340,824 sites | 1.43% |
| sites monomorphic across all 88 haplotypes | 4,713 | 340,824 sites | 1.38% |
| sites with no called haplotype at all | 1 | 340,824 sites | 0.00% |

- Columns: counts, then the denominator each share is taken over. Slots and sites are
  different denominators and must not be mixed.
- Rows: "declared" means the allele index is below `num_alleles[site]`, i.e. the VCF lists
  it. "Zero cohort support" means no haplotype of the 88 carries it.
- Synthesis: v1's output support is cohort-supported paths, so **the sampler cannot emit
  these 8,823 alleles.** On 1.43% of sites the unreachable allele is REF, so a target
  carrying reference there cannot be matched by any donor state. The consequence is a
  `phi_t` penalty — **pure utility loss, and bounded at about 2% of sites. No privacy
  consequence**: `Y_D` is fixed and target-independent, so common support (condition 1)
  holds regardless. Monomorphic sites are released deterministically, which is also safe
  because the determinism is identical for every target and cancels in the likelihood
  ratio; determinism only breaks the bound when it depends on the input.

**Framing to keep.** In PanMixer the target is DATA — a column inside the matrix every
cohort statistic is computed over, which is how departure (c) arises. In v1 the target is
a PARAMETER: it never enters the matrix, `support_D`, or the state set, and touches the
mechanism only through `phi_t` comparing the target's block against each donor's, inside `psi_t`. Adding it to the cohort would
break the `R_D` cancellation in Step 3 of the Theorem 1 proof — not a weaker bound, no bound.

**Three senses of "target-unique", only one of which is a real gap.** (i) An allele absent
from `G`: unreachable by the mechanism because `p = Map(g, G)` is a traversal of `G`; lost
at mapping time, a limitation of the setting. (ii) A combination unique to the target: the
ordinary case, and exactly what the mosaic produces — conditional on the transition
constants actually recombining (exit gate 3). (iii) An allele in `G` with no cohort
support: the measured 1.15% above.

## 2026-09-27 — A-DAT-lv-nesting: the biallelic SNP mask is blind to snarl level

**Question.** Does `build_biallelic_snp_mask.py` (Step 11) use the `LV` tag? No — and
no PanMixer Python reads the INFO column at all (`grep -rn "LV=\|INFO" --include=*.py`
over the pinned checkout returns nothing). Step 11 re-reads `pangenome.vcf.gz` directly,
so INFO is present in the line it parses; it splits with `maxsplit=5` and tests only REF
and ALT.

**Method.** One pass over `external/PanMixer/starting_data/chr21/pangenome.vcf.gz`,
parsing `LV=` out of INFO per record, cross-tabulated against
`biallelic_snp_mask.npy`. All 340,824 records carry an LV tag; none missing.

| snarl level | all records | share | inside the mask | share of mask |
|---|---|---|---|---|
| LV=0 (top level) | 309,335 | 90.76% | 252,039 | 91.83% |
| LV=1 | 28,029 | 8.22% | 20,153 | 7.34% |
| LV=2 | 3,228 | 0.95% | 2,141 | 0.78% |
| LV=3 | 222 | 0.07% | 119 | 0.04% |
| LV>=4 | 10 | 0.00% | 6 | 0.00% |
| **nested (LV>0)** | **31,489** | **9.24%** | **22,419** | **8.17%** |

- Columns: all chr21 records, versus the 274,458 records the mask marks True.
- Rows: `LV` is the snarl-tree level written by `vg deconstruct -a`, 0 = top level.
- Synthesis: **9.24% of records are nested.** This independently reproduces the 9.2%
  recorded on 2026-09-22 for records whose POS lies inside an earlier record's span —
  two different methods (tag-based here, coordinate-based there) agreeing to rounding,
  which is a useful cross-check on both.

**Why it matters.** The mask is the site set for the gap-score attack (Step 21), which
assigns each site a Hardy-Weinberg genotype probability and sums `-log` across sites,
treating loci as independent evidence. 8.17% of that set is nested, and a nested allele
is conditional on its parent traversal, so those sites are not independent evidence.
The published Methods DO state a conditional-independence approximation for nested
variants — `p_hat(a_nested) ~ p(a_parent) * p(a_nested)` — but on the SAMPLING side
only; nothing equivalent appears in the scoring or attack path. **Magnitude unmeasured:**
whether this materially moves the gap score is not established here, and the
correlation may be weak in practice. Structure only.

**For v1.** `LV` and `PS` are two of the fields the REASON array should carry (see
`missing_policy`): they are what distinguishes "this donor's parent allele does not
contain the child bubble" from the other four causes of a `-1`.

## 2026-09-26 — A-DAT-missingness-scope: what 11.31% and 0.35% are means OVER

**Why.** `docs/terms.md` (`missing_policy`) and `docs/plan.md` quoted "11.31% versus
0.35%" with no measurement entry behind either, and the primer then described 11.31%
as the mean "across every variant" / "once every variant is a chain position". A
second-round primer audit re-measured it; re-measured again here, independently.

**Method.** `external/PanMixer/starting_data/chr21/pangenome.npy` (44 x 340,824 x 2, int16), missingness
per record = share of the 88 haplotype entries equal to `-1`, averaged over the class.
Anchors = the 273,475 keys of `pangenome_to_thousand_g_alignments.pickle` (the PanGenie
overlap, Correction 5). Block membership from `blocks_dict.json` (100,757 entries, every
record covered exactly once). Target included (the matrix as shipped).

| position class | n | share of variants | mean missingness | records > 25% missing |
|---|---|---|---|---|
| all variants | 340,824 | 100.00% | 4.263% | 30,800 |
| anchors, all | 273,475 | 80.24% | 0.351% | 925 |
| non-anchors, all | 67,349 | 19.76% | 20.150% | 29,875 |
| variants in multi-variant blocks | 251,737 | 73.86% | 1.582% | 8,108 |
| non-anchors in multi-variant blocks | 29,966 | 8.79% | **11.306%** | **7,555** |
| anchors in multi-variant blocks | 221,771 | 65.07% | 0.268% | 553 |
| variants in blocks with >= 2 anchors | 243,327 | **71.39%** | 0.501% | 1,274 |
| anchors in blocks with >= 2 anchors | 221,319 | 64.94% | **0.264%** | 538 |
| variants in multi-variant blocks with <= 1 anchor | 8,410 | 2.47% | 32.853% | 6,834 |
| singleton blocks | 89,087 | 26.14% | 11.840% | 22,692 |

- Columns: n and share are record counts; missingness is the class mean of the
  per-record missing share; the last column counts records over a quarter missing.
- Rows: `anchors in blocks with >= 2 anchors` is what PanMixer's HMM actually steps
  over. `all variants` is what v1 steps over, because v1 sends singletons through the
  same sampler as T = 1 chains. `variants in blocks with >= 2 anchors` is what PanMixer
  SAMPLES by donor copying today; everything else is drawn from allele frequencies.
- Synthesis: **11.31% is the mean over the 29,966 non-anchor records inside
  multi-variant blocks** — the positions v1 adds INSIDE blocks PanMixer already chains.
  It is not the mean over every variant (4.26%), and 0.35% is the mean over all
  anchors, not over the anchors PanMixer's HMM steps over (0.26%). The honest
  like-for-like is **v1's positions 4.26% (30,800 over 25%) against PanMixer's HMM
  positions 0.26% (538 over 25%)**, about 16x (derived). The "32x" (11.31 / 0.35) is a
  ratio of two differently scoped means and is quotable only with both scopes named.

**Second finding — the donor-copied share today is 71.4%, not 73.9%.** 73.86% is the
share of variants in multi-variant blocks. But the sampler sends a multi-variant block
with <= 1 anchor to the allele-frequency route too (8,410 variants, 2.47%), so the share
PanMixer samples by copying a donor is **71.39%** — the figure the 2026-09-18 entry
already gives ("10.4% of blocks ... 71.4% of variants"). `docs/plan.md` v1 REPLACE
item 2 said 73.9%; corrected there the same day. The 8,410 are also the worst-behaved
class in the table (32.85% mean missingness), which is presumably why they have <= 1
anchor.

## 2026-09-26 — A-DAT-primer-provenance: four primer numbers that `docs/` never logged

**Why this entry exists.** A three-critic audit of the expanded primer flagged
`274,458`, `177,865` and `3,077.60` as fabricated, and `8,661` as invented precision,
because `grep` over `docs/` returned nothing for any of them. All four were already in
the primer as committed on 2026-09-23. Each was re-measured from its artifact before
anything was deleted. **All four are exact.** The defect was that they were measured
during a session and quoted in teaching material without being logged here. They are
logged now, so they are quotable.

**1. Biallelic SNP records on chr21: 274,458 of 340,824 (80.53%).** Measured by
loading `external/PanMixer/starting_data/chr21/biallelic_snp_mask.npy` (dtype bool,
shape (340824,)) and summing it. The mask is built by
`external/PanMixer/starting_data/scripts/src/build_biallelic_snp_mask.py`: a record
qualifies iff it has exactly one ALT and both REF and ALT are a single base in ACGT.
**Relation to the other counts:** 65,176 records are indel/structural (see the
2026-09-23 conflicts entry), so 275,648 records are single-base; 274,458 of those are
biallelic, so **1,190 are multi-allelic SNP records**. The two numbers are consistent;
one is a subset of the other. `238,052` (69.8%) is a DIFFERENT quantity — the sites
the gap-score attack keeps after the mask AND the 1000G match — and must not be quoted
as the mask count.

**2. The capacity-0.5 obfuscation run, HG00438, chr21, seed 123.** From
`logs/pm_obf.log` (sentinel `logs/pm_obf.status` = 0), verbatim optimizer and stacker
lines:

| metric | capacity 0.1, seed 123 | capacity 0.5, seed 123 |
|---|---|---|
| selected moves | 47,070 | 177,865 |
| utility_loss / max | 928.92 / 9288.85 (10.0%) | 3,077.60 / 9288.85 (33.1%) |
| pmi_gain (nats, as reported) | 115,315.84 | 120,443.42 |
| alleles changed / total | 68,335 / 650,541 (10.50%) | 106,549 / 649,242 (16.41%) |

- Columns: two capacities from the same run script, same target, same seed.
- Rows: `moves` and `utility_loss` are the optimizer's own log line; `changed` is
  the stacker's. The allele denominators differ slightly between runs (650,541 vs
  649,242); the log does not say why, and it is recorded rather than explained.
- Synthesis: raising the budget five-fold (0.1 -> 0.5 of max loss) spends only 33.1%
  and buys 4.4% more reported privacy (120,443.42 / 115,315.84 = 1.0445). The knapsack
  saturates. Seed 456 at capacity 0.1 gave 47,069 moves and 928.88 — seed-stable.

**3. PanMixer's Python: exactly 8,661 lines.** `git ls-files '*.py' | xargs cat | wc -l`
inside the pinned checkout. Every `.py` on disk is tracked, so the untracked count is
the same. `docs/data.md` previously rounded this to ~8,660.

**4. The published Equation (5) carries its minus sign.** Rendered
`paper/s41467-026-77591-0_reference.pdf` **PDF page 10** (not 9) at 300 dpi:
`L(g_i, g*_target) = − Σ log f_j`, the sum running over the loci j in `S(g_i, g*_target)`. The minus is typeset and
unambiguous. The "sign discrepancy" the primer once carried as unresolved was an
artefact of text extraction dropping the glyph (this PDF has no text layer at all), not
a defect in the paper. With the sign, the equation agrees with the paper's own prose
("emphasizes rare allele matches"). Same page states the read-mapping donors as
HG00138, HG00635, HG01112, HG02698 and **NA18853**, and the HMM constants
Ne = 10,000, r = 1.26, Δx in cM, d = Δx·Ne·r, p = e^(−d/n) + (1 − e^(−d/n))/n.

**Lesson.** "Not in `docs/`" means *unlogged*, not *false*. A checker that cannot
find a number must re-measure it from the artifact before recommending deletion; the
audit here would otherwise have deleted four true measurements and published a
retraction of a correct one. See the matching entry in `docs/bugs.md`.

## 2026-09-25 — A-DAT-missingness: `-1` conflates at least four different things

**Goal.** Dylan asked how a donor can fail to traverse a parent bubble, given that a
variant ought to have a value for every haplotype. Establish what `-1` actually means
before choosing a `missing_policy`.

**The answer for TOP-LEVEL variants is that it cannot.** Every haplotype path runs
along the GRCh38 backbone, so a top-level variant has an allele for everyone; the
2.25% missingness there is assembly gaps and conflicts, not structure.

**For NESTED variants it can, and routinely does.** A child bubble sits inside a
PARTICULAR parent allele. A haplotype taking a different parent allele never visits
that part of the graph, so there is no allele to report. It is NOT APPLICABLE, not
missing — and VCF has no way to say so, therefore `.`, therefore `-1`.

Clean measured example: parent row 37624 at POS 13,205,133 (REF 13 bp; five ALTs of
15, 1, 17, 11 and 13 bp) with child row 37625 at POS 13,205,136, three bases inside
it. The 87 haplotypes carrying parent alleles 0, 1, 3, 4 or 5 all have a child
allele. The one haplotype carrying **allele 2 — the 1 bp allele, i.e. the region is
deleted** — has `-1`. The child variant does not exist on that haplotype because the
sequence it sits in was deleted.

**How much of nested missingness this explains**, over 5,735 nested records with any
missing child call:

| cause | share |
|---|---|
| entirely inherited from a missing PARENT call | 12.1% |
| among haplotypes WITH a parent call, perfectly determined by which parent allele they carry | **66.1%** |

**A second, unrelated cause: assembly gaps.** Run lengths of consecutive missing
positions (strand 0, all 44 subjects): 28,052 runs, median 1, mean 20.0, max 23,478.
**54.4% are isolated single positions** (structural / not-applicable), while **349
runs of length >= 100 cover 448,283 entries** — those are contig breaks, not structure.

**A FIFTH cause, from the converter rather than the data.**
`external/PanMixer/tools/common/VCFtoNP.py` maps a HAPLOID genotype (one allele, no
`|`) to `-1` on strand 1. Measured on 150,000 chr21 records: **3.135% of sample GT
fields are haploid**, and 5.581% of the two allele slots are an explicit `.`. So some
of the matrix's missingness is a conversion decision, not a property of the data.

**The causes are distinguishable — until conversion throws them away.** `LV` and `PS`
separate "not applicable" from "inherited"; run-length separates assembly gaps; the
`CONFLICT` tag names conflicts; GT arity identifies the haploid case. But `VCFtoNP`
reads only `fields[1]` (POS) and `fields[9+i].split(':')[0]` (the GT subfield) — the
entire INFO column is discarded. By the time anything reaches the matrix, all five
look identical, and `1/support(v)` and `pmi[hap == -1] = 0` treat them alike.
**Recommendation for v1: emit an auxiliary REASON array during conversion.** Every
cause is recoverable at that moment and none is recoverable afterwards.

**So `-1` means at least five different things**: (1) not applicable, the parent
allele does not contain this bubble; (2) inherited from a missing parent; (3) an
assembly gap; (4) a CONFLICT, where the sample had multiple paths that disagreed.
A single `missing_policy` is being asked to cover all five, and they do not want the
same treatment: "not applicable" argues for renormalise, an assembly gap arguably
argues for something closer to wildcard. **Consider splitting the policy by cause —
the LV/PS tags and run-length both distinguish them, and both are currently discarded
at VCF-to-numpy conversion.**

**Provenance.** chr21, `external/PanMixer` @ `c182c38`, measured 2026-09-25.


## 2026-09-25 — Design note: if the chain comes from the graph, what happens to LD blocks

**Verified fact.** PanMixer never loads a genetic map. Searching the whole codebase,
every `recomb` hit is an output FILENAME for the RECOMB conference; the only genetic-map
references are the unwired `external/PanMixer/starting_data/scripts/get_genetic_maps.sh`
(which has a broken shebang, `!/bin/bash`, and is referenced only by `.gitignore`).
Its transition "distance" is a VCF row index.

**That explains the architecture.** With no recombination map, the LD blocks have to
carry the recombination information structurally: the per-block HMM treats everything
inside a block as one linked unit and everything across blocks as independent. Blocks
are a discretised stand-in for a distance-dependent transition model.

**Consequence for our design.** In a Li-Stephens copying model the linkage is carried
by the DONOR PANEL and the TRANSITION PROBABILITIES, not by where block boundaries
fall — copying a stretch from donor j reproduces whatever allele combinations j
carries, including rare co-occurrences that exist only because they sit on one
ancestral chromosome. So with a genome-wide chain plus distance-dependent transitions,
LD blocks become redundant for the MODEL. The two jobs they do for PanMixer both
disappear for us: they scope its per-block HMM (we have one chain), and they are the
unit of its knapsack (we have no knapsack).

**Practical requirement this creates.** Distance-dependent transitions need a genetic
map in GRCh38 coordinates. The map repo PanMixer's unwired script points at
(the joepickrell 1000-genomes-genetic-maps repo) offers `interpolated_from_hapmap` and
`interpolated_OMNI`, both GRCh37-era. SHAPEIT4 publishes b38 maps (URL responds 302),
which is the likely source. **Not yet acquired or verified — an open task, not a
decision.** Falling back to bp distance with a constant rate is possible but crude.

**What LD blocks would still be for.** Evaluation comparability with PanMixer
(reporting per-block results on their partition), and possibly grouping for the
`beta_t` weights. Neither is a modelling need.


## 2026-09-25 — A-DAT-anchors: why anchors exist, and whether we need them

**Goal.** Dylan asked why the HMM runs only over anchor variants rather than every
variant in a block, and whether we can drop the concept. Establish the rationale
before deciding.

### The paper's rationale is about TOP-LEVEL SNPs, not external matching
Published Methods, on choosing the branch: *"dependent on the number of top-level
SNPs that exist inside the LD block. If 2 or more exist, enough SNPs exist to build
an HMM, and we can synthesize more accurate new haplotype blocks. If only one
top-level SNP exists, we do not have enough useful information to build an HMM."*
Top-level means on the GRCh38 backbone path (LV=0). The >=2 threshold is simply
that a Markov chain needs two positions to have a transition between them.

### Why restrict to top-level SNPs — the data supports it
Missingness per record class on chr21 (fraction of the 88 haplotypes with no called
allele):

| class | n | mean missing | records >50% missing |
|---|---|---|---|
| top-level SNP (the paper's V_SNPs) | 253,000 | **2.25%** | 421 |
| top-level non-SNP | 56,335 | 3.11% | 722 |
| **nested (LV>=1)** | 31,489 | **22.53%** | 3,515 |

Nested variants are ~10x more often uncalled, because a donor whose path does not
traverse the parent bubble has no meaningful allele at the child — and the matrix
cannot distinguish "not applicable" from "missing". Using them as chain positions
would feed the emission model noise. **That is a real justification.**

Two candidate rationales that do NOT hold:
- **Compute.** Anchors 273,475 vs all variants 340,824 — only **1.25x** more
  positions (24.1M vs 30.0M operations at K=88). Negligible.
- **Multi-allelic emissions.** Median ALT count is 1 even for non-SNP records; the
  90-allele tail is real but rare.

### But the CODE does not implement the paper's set
| | records |
|---|---|
| paper's V_SNPs (LV=0 and SNP) | 253,000 |
| code's anchors (PanGenie overlap) | 273,475 |
| overlap | 234,532 |
| **anchors that are NOT top-level SNPs** | **38,943** (14.2% of anchors) — of which **13,436 are NESTED** and 25,507 are non-SNP |
| top-level SNPs that are NOT anchors | 18,468 (7.3%) |

So the implementation admits 13,436 nested variants as chain positions — precisely
the class the paper's own rationale excludes — while dropping 18,468 top-level SNPs
that qualify. The PanGenie overlap is a proxy for "well-behaved position" that is
neither necessary nor sufficient.

### Implication for our design
Drop the PanGenie-overlap criterion: it has no principled basis and it is what makes
rare non-anchor variants invisible to `eps_j` (measured earlier: in 32% of blocks the
ignored non-anchor rarity exceeds `eps_j` itself). Do NOT simply promote every
variant to a chain position — the 22.5% missingness at nested records is a real
obstacle. The principled replacement is to take positions from the **graph
structure** (top-level / snarl-tree) rather than from an external callset, which is
the same answer open question 8 is circling for `T`.

**Provenance.** chr21, `external/PanMixer` @ `c182c38`, published paper Methods,
measured 2026-09-25.


## 2026-09-24 — A-DAT-noop-moves: 31% of selected obfuscation moves change nothing

**Method.** Walked block 1077 end to end for HG00438 on chr21, then swept the whole
capacity-0.1 run comparing each selected move's released alleles to the original.

**Block 1077** (rows 34816-34820, POS 13,000,333-13,000,621): 5 variants, all 5
anchors, so the HMM branch. Target carries `[0,0,0,0,0]` on both strands.
- `eps_j` = **0.0711 nats** per strand, computed from the target's ORIGINAL block.
- support(v) = 87 of 88 for every variant, so 1/support = 0.01149 each.
- The sampled replacement came back `[0,0,0,0,0]` — **identical to the target**,
  because the degenerate transitions copy one donor verbatim and that donor happens
  to carry the same common alleles.
- Changed variants = 0, so **`eta_j` = 0.00000** (it would be 0.05747 if all 5 changed).

So the move offers positive privacy at **zero** utility cost. The LP takes it
unconditionally. Confirmed in the run: `xsol[625]` is 1 on both strands, and the
released alleles at those rows are unchanged.

**Swept across the whole run: of 47,070 selected moves, 14,526 changed NOTHING
(30.9%).** Every one contributed its `eps_j` to the reported `pmi_gain` of
115,315.84 while altering no allele and costing no utility.

### Follow-up — the inflation is 7% of eps, not 31%
The 31% figure counts MOVES; it is not the share of privacy credit. Measured on a
3,000-move sample of the 47,070 selected:

| selected moves | n | median eps_j | mean eps_j |
|---|---|---|---|
| changed something | 2,112 | 2.1904 | 3.3103 |
| **changed nothing** | 888 | **0.3898** | 0.5964 |

**No-op moves carry only 7.0% of the total eps.** They concentrate on LOW-eps blocks
— common patterns that many donors share, so a randomly drawn donor often carries
the target's own alleles. Block 1077 is typical: 82 of 88 donors carry exactly the
target's pattern, giving eps_j = 0.0711, and -log(82/88) = 0.0706 confirms the
reading of eps_j as the surprisal of the target's block under the cohort.

Note also what is NOT inflated: for a no-op, `eta_j` is correctly **0**, because it
sums `1/support(v)` only over variants that differ. So the utility axis is honest;
it is the privacy numerator that gains ~7% for nothing, and the privacy-per-utility
ratio that is consequently overstated.

**Consequence.** `pmi_gain` — the numerator of the reported "Privacy Risk" axis —
is inflated by moves that did not touch the genome. The effect is systematic rather
than incidental: it follows directly from `eta_j` being computed as
`sum(utility_loss * (new != orig))` while `eps_j` is always positive and is computed
from the original block regardless of what was sampled. It compounds the sampler
degeneracy: the more often the sampler returns the target's own alleles, the more
free privacy the LP books.

**Provenance.** chr21, HG00438, capacity 0.1, seed 123, `external/PanMixer` @
`c182c38`, measured 2026-09-24.


## 2026-09-24 — A-DAT-af-table: where f_v comes from, and two corrections

**Goal.** Dylan asked what population the allele frequencies in `eps_j` are computed
from, and whether "no anchors" implies a degenerate frequency.

**The table is built by a three-way priority** per variant
(`external/PanMixer/starting_data/scripts/src/get_af.py`:27-78). Measured on chr21:

| branch | source of counts | target included? | variants |
|---|---|---|---|
| 1. variant matches the 1000G map | **1000G only** | not under Phase 3; **YES under our 30x repin** | 258,610 (75.9%) |
| 2. matches PanGenie but not 1000G | PanGenie **plus all pangenome haplotypes** | yes | 25,325 (7.4%) |
| 3. novel (neither) | pangenome only, then REF padded to `NOVEL_DENOMINATOR` | yes | 56,889 (16.7%) |

`NOVEL_DENOMINATOR = 2 * (N_1000G + N_pangenome)` = 2*(3202+44) = **6,492** under our
repin, so a novel allele carried by one haplotype gets f_v = 1.54e-4, NOT 1. The
padding exists precisely to stop a novel allele looking common because only its own
carriers were counted; the REF allele absorbs it.

**"Anchor" and "has population data" are different things.** The anchor set is the
PanGenie map; the AF priority also consults the 1000G map. **10,460 non-anchor
variants still receive 1000G frequencies.** A block with zero anchors can be fully
populated with real frequencies.

Median −log f_v where HG00438 carries a NON-REF allele: branch 1 = 0.63 nats,
branch 2 = 0.43, branch 3 = **5.41**. At novel sites the target usually carries REF
(median f_v 0.998), so most novel sites contribute ~0 — the signal is in the 11,760
novel sites where it carries a non-REF allele.

### Correction 7 — the `#FIX ME` zeroing does NOT fire on our data
[INCORRECT] - Alleles like these have population frequency 0 in any external callset, so `-log(0)` -> infinity -> zeroed by that `#FIX ME`. **PanMixer assigns the most identifying sites on the chromosome a privacy score of exactly zero.**
[CORRECTION - 2026-09-24]: Measured for HG00438 across all 319,092 chr21 sites where it carries a called allele: **ZERO have f_v = 0**, so the `#FIX ME` branch never executes. The reason is that the target's own allele is counted in every branch — via the pangenome counts in branches 2 and 3, and (under our 30x repin) via the panel itself in branch 1, since 39 of 44 HPRC donors including HG00438 are in the 30x panel. The zeroing is a LATENT defect whose trigger is panel-dependent: under the shipped Phase 3 panel, where 0 of 44 HPRC donors appear, a branch-1 variant whose allele no 1000G sample carries WOULD give f_v = 0 and fire it. The conclusion that the hypervariable site contributes zero privacy is still correct, but the MECHANISM is different — see Correction 8.

### Correction 8 — the hypervariable site scores zero by non-anchor exclusion, not by zeroing
The 257 kb record (row 51791) gets a perfectly good frequency: HG00438 carries
alleles 89 and 88, each f_v = 1.540e-4, i.e. **8.78 nats** of self-information. That
value is never used. The record is **not an anchor** (`r in align` is False), and its
block 1639 has 92 records with 82 anchors, so the block takes the **HMM branch** —
where the forward algorithm scores only the anchors. The most identifying record on
the chromosome contributes nothing because it is not in the PanGenie callset, which
is a completely different failure from the `-log(0)` zeroing.

### New risk — our own repin weakened target-independence
Branch 1 draws counts from the 1000G panel alone, and the code comment asserts "None
of the pangenome subjects appear in 1000g_phased" — true for Phase 3 (0 of 44
overlap, measured) but **false for our 30x substitution (39 of 44 overlap)**. So our
rebuilt `allele_frequencies.npy` includes the target's own alleles in the
frequencies used to score it. For PanMixer that is merely a stronger version of an
existing self-inclusion. [INCORRECT] - **For US it violates constraint 1**, so any leave-one-out
cohort must also drop the target from the 1000G panel, not only from the pangenome.
[CORRECTION - 2026-09-28]: The measurement is right and the implication is wrong. 39 of 44
HPRC donors including HG00438 really are in the 30x panel — that is true of the pipeline as
built, and it matters for the BENCHMARK, where a cohort member stands in as target. It is
not a constraint-1 violation of our MECHANISM, whose target is external to every panel by
construction. Two further reasons this is now moot: under model C the mechanism uses no
external panel at all, and it computes no allele-frequency table (decided 2026-09-28). The
overlap still needs removing from the ATTACK DATABASE — and note the fix is incomplete, because
all 39 overlaps lie in the 698 samples the 30x release ADDED beyond Phase 3 to complete trios,
so dropping the 39 named samples leaves their first-degree relatives in.

**Provenance.** chr21, HG00438, `external/PanMixer` @ `c182c38`, measured 2026-09-24.


## 2026-09-23 — A-DAT-conflicts: what CONFLICT records do downstream, and how blocks handle spanning variants

**Goal.** Dylan asked whether all CONFLICT records behave alike, what they imply for
blocking, and what actually happens in the released output when a long variant and
records inside its span fall in different blocks.

### Finding 1 — a CONFLICT becomes MISSING DATA, and missingness is not neutral
At a CONFLICT record vg writes that sample's genotype as `.|.`, so `VCFtoNP` stores
**-1**. Measured on chr21:

| | value |
|---|---|
| CONFLICT records | 1,855 (0.54%) |
| mean haplotypes missing at a CONFLICT site | **29.4 of 88** |
| mean haplotypes missing elsewhere | 3.6 of 88 |
| matrix entries missing overall | 1,278,653 (4.26%) |
| sites with >=1 missing entry | 52,926 (15.5%) |

Missingness then drives BOTH of PanMixer's per-move scores, in the same direction:
- **Privacy -> 0.** `external/PanMixer/starting_data/scripts/src/get_support_and_pmi.py` zeroes the score wherever the target's
  haplotype is -1 (`pmi_1[haplotypes[:,0] == -1] = 0`).
- **Utility cost -> up.** support(v) counts non-missing entries, so 1/support(v)
  rises. Median utility weight at CONFLICT sites 0.01818 vs 0.01136 elsewhere —
  **1.60x more expensive to change**.

**Consequence: these sites are worth zero privacy and cost more, so the knapsack
will essentially never select them. They are released VERBATIM.** That is a leak
concentrated exactly where the target's structure is unusual.

### Finding 2 — conflicts are NOT uniform; they concentrate in complex regions
| property | CONFLICT records | all records |
|---|---|---|
| nesting LV=0 | 59.1% | 90.8% |
| nesting LV=1 | **37.5%** | 8.2% (a 4.6x enrichment) |
| nesting LV=2 | 3.4% | 0.9% |
| median size (max REF/ALT) | 4 bp (mean 602, max 311,406) | mostly 1 bp |
| median alleles per record | 4 (max 81) | 2 |

They also land in the SIMPLE part of the block structure: the 1,855 records occupy
1,634 blocks, of which **1,469 are singletons** and only **1.8% are HMM-path blocks**
(baseline 10.4%). So they are mostly alone in their own block and mostly handled by
the independent per-variant allele-frequency sampler, not the haplotype model.

### Finding 3 — the read-mapping evaluation never sees any of this
`external/PanMixer/tools/downstream/utility_out/vg_prep.py` builds the evaluation graph with `bcftools view -v snps`, which keeps
288,020 of 340,824 chr21 records (84.5%) and drops the 52,804 indel/structural
records where conflicts and nesting concentrate. So PanMixer's read-mapping utility
metric is computed on a graph that excludes the hardest regions by construction.
`af_loss`, `ld_loss` and the gap score read the matrix/VCF directly and DO see them.

### Finding 4 — what the released output actually is when a variant spans blocks
Read from `external/PanMixer/tools/panmixer/stacker.py`:37-70. The output is `new_haplotypes.npy`,
shape (n_sites, 2) — **one allele index per RECORD**. It is built as
`new_variants = original_haplotypes.copy()`, then for each selected (block, strand)
the record indices in that block are overwritten from the presampled replacement.

So there is no sequence assembly step, and no coherence check. A long parent record
and the records positionally inside its span are decided **independently**: if the
parent's block is not selected the parent keeps the target's TRUE allele, while
child records in a selected block are replaced. The released VCF then asserts both
at once. Whether those assignments jointly correspond to a walk through the graph is
never computed, checked, or asked.

### Finding 5 — the MECHANISM sees everything; only the EVALUATION filters
An earlier note implied indel/structural records are excluded. That is wrong as a
statement about the mechanism. Measured on our HG00438 chr21 run at capacity 0.1:

| stage | records seen | non-SNP | CONFLICT |
|---|---|---|---|
| mechanism (PMI + knapsack + stacker) | 340,824 (100%) | 65,176 | 1,855 |
| `af_loss` / `ld_loss` (strict 1000G mask) | 258,610 (75.9%) | 21,550 | 37 |
| gap score (posref mask) | 268,851 (78.9%) | 29,965 | 110 |
| gap score after the biallelic-SNP mask | 238,052 (69.8%) | **0** | 49 |
| read mapping (`bcftools view -v snps`) | 288,020 (84.5%) | mixed records kept | **1,167** |

So indels and structural variants ARE scored, selected and rewritten. In that run
**15,206 non-SNP records changed — a 23.3% rate, HIGHER than the 14.9% for SNPs** —
and the 257 kb hypervariable record was rewritten on both strands
(alleles [89, 88] -> [35, 64]).

**Sequence actually rewritten: 3,129,810 bp = 3.13 Mb, 6.7% of chr21.** 222 changed
alleles exceed 1 kb; the largest single rewrite is 257,485 bp. PanMixer reports the
same run as "68,335 of 650,541 alleles changed (10.5%)" — because `1/support(v)` is
per RECORD, the 3.13 Mb never enters the accounting.

Correction to an intermediate figure used while investigating: dropping non-SNPs
does NOT remove most conflicts. `bcftools view -v snps` keeps mixed multi-allelic
records, and conflicts are enriched in those (median 4 alleles), so **1,167 of
1,855 conflicts survive** the read-mapping filter — 37% are dropped, not 84%.

**Net effect on the scoring.** The biggest edits by base pairs land on exactly the
records the read-mapping metric cannot evaluate (its graph is SNP-only), while
`af_loss` and `ld_loss` see 75.9% of records and the gap score's biallelic path
sees no non-SNP records at all. No metric in the suite is sensitive to the 3.13 Mb.

**Provenance.** chr21, `external/PanMixer` @ `c182c38`, measured 2026-09-23.


## 2026-09-22 — A-DAT-hypervariable-sites: measured the worst-case site class, flagged as a risk and an opportunity

**Goal.** Dylan asked whether hypervariable structural sites could affect how well
our mechanism performs, and whether they are somewhere we could improve on
PanMixer. Measure the class before reasoning about it.

**Method.** Found the most multi-allelic site on chr21 and characterised it and its
class from the VCF, the graph traversals (`AT` field) and the numpy matrix.

**Result — the worst site on chr21 is `grch38#chr21:14569980`, graph bubble
`>102270111>102277685`, LV=0 (top level).**

| property | value |
|---|---|
| alleles declared | 90 (89 ALT) |
| distinct alleles carried by the cohort | **88, across 88 haplotypes — every haplotype unique** |
| REF length | 257,439 bp |
| ALT lengths | 257,124 – 514,768 bp |
| size of the single VCF line | **27,983,949 bytes (~28 MB)** |
| `AC` for allele 1 | 0 — an orphan left by `view -s ^chm13` without `--trim-alt-alleles` |
| nested child records (LV>=1, PS = this bubble) | **2,464** |
| total matrix rows for this one region | **2,485** (0.7% of chr21's 340,824 rows) |
| graph nodes in the bubble | 7,575 distinct; traversals 4,971–10,008 nodes |
| allele 58 | 10,008 nodes over 5,014 distinct = ~4,994 revisits -> a **tandem duplication** of the whole region |
| inversions | none (0 traversals use reverse orientation) |

Class frequency on chr21: **7** sites with >=88 alleles, **155** with >=44, **618**
with >=20.

**Why `vg deconstruct -a` emits both levels.** The nested decomposition DOES happen
— 2,464 children exist. The top-level record is emitted IN ADDITION because some
alleles cannot be expressed as a combination of the children: allele 58 traverses
the region twice, and no set of independent per-site substitutions can encode "this
whole region is duplicated". VCF alleles are also flat strings with no syntax for
"allele = this combination of child variants", so the top-level traversals must be
spelled out in full. Hence 28 MB.

**Provenance.** `external/PanMixer/starting_data/chr21/pangenome.vcf.gz` ·
`num_alleles.npy` · `pangenome.npy` · measured 2026-09-22.

### Why this matters to US — four consequences, recorded as risks to watch
1. **Double counting.** The same sequence is represented twice: once in the
   top-level record's giant allele strings and again across its 2,464 children. All
   2,485 are rows in the matrix, so a target's variation here contributes to
   `eps_pmi` at BOTH levels. The paper's nested approximation
   p(nested) ~ p(parent) x p(nested) addresses child-vs-parent dependence but not
   the fact that the parent is itself a scored row.
[INCORRECT] - 2. **The most identifying sites score ZERO.** Every allele here has population frequency 0 in any external callset, so `-log(0)` -> infinity -> zeroed by the `#FIX ME` in `get_support_and_pmi.py`. The sites where each haplotype is unique — a perfect fingerprint — contribute nothing to PanMixer's privacy score.
[CORRECTION - 2026-09-26]: The conclusion stands and the mechanism is wrong. The zeroing branch never fires on this data (Correction 7: 0 of 319,092 sites have f_v = 0, because the target's own allele is counted in every frequency branch). These sites score zero because the privacy score reads anchor alleles only and they are not anchors (Correction 8). Marked in place on 2026-09-26; the copy under Correction 7 was made without marking this original.
3. **Utility weighting ignores scale.** `1/support(v)` is per VARIANT. A move here
   rewrites up to half a megabase yet is costed like a one-base SNP (and, since
   support is 1, is among the most expensive — but for the wrong reason).
4. **Array width is set by the worst site.** `allele_frequencies.npy` is
   (340824, 90): every site is padded to the widest one in the chromosome.

### Opportunity for OUR mechanism — to evaluate, not yet a claim
These sites are exactly where a bounded, path-level guarantee should beat a
per-block edit rule, because our `tau` bound holds regardless of how rare an allele
is, whereas `eps_pmi` collapses to 0 precisely here. Three things must be checked
before claiming any advantage:
- whether `output_support` can contain a *valid* traversal here at all (the
  AF-sampling branch draws variants independently and would produce combinations no
  haplotype carries — see the graph-validity open question);
- whether `phi_t` should be length-weighted rather than per-variant, given a single
  block can span 257 kb;
- how to avoid double-counting the top-level record and its 2,464 children in both
  `u_path` and any privacy measurement.
**Nothing here is measured about our mechanism — it does not exist yet.**

### Follow-up 2026-09-22 — the overlap is STRUCTURAL, not an edge case
Dylan asked for the fraction of chr21 affected. Measured over all 340,824 records,
taking each record's span as [POS, POS + len(REF) - 1]:

| measure | value |
|---|---|
| records whose POS lies inside an earlier record's span | 31,487 (**9.2%**) — matches the LV>=1 nesting fraction exactly |
| ...inside a record of REF length >= 10 kb | 19,211 (5.6%) |
| ...inside a record of REF length >= 100 kb | 4,608 (1.4%) |
| blocks containing >=1 overlapped record | 15,540 of 100,757 (**15.4%**) |
| variants living in an affected block | 78,937 (**23.2%**) |
| HMM-path blocks (>=2 anchors) affected | 1,229 of 10,502 (**11.7%**) |

**Verdict: not documentable as an edge case.** Roughly a tenth of records and a
quarter of variants-by-block are involved, so block disjointness — which both the
knapsack's additive accounting and our FFBS chain assume — fails broadly, not locally.

**The duplication's structure, measured.** Allele 58 of the parent bubble walks
10,008 node steps over 5,014 distinct nodes; 4,994 nodes appear exactly twice and 20
appear once. At step ~5004 the walk reaches `>102277684` (beside the closing
boundary) and jumps back to `>102270112` (beside the opening boundary). So it is a
**whole-region tandem repeat**, and the two copies are NOT identical — the 20
singleton nodes are ~10 interior positions where the copies diverge. That is exactly
why the interior variants cannot be given one allele per haplotype for this sample,
and why `CONFLICT=HG01891` appears at POS 14577263.

**Why VCF cannot express this compactly.** The graph already stores it compactly —
the traversal simply visits the nodes twice, duplicating no sequence. The ~28 MB is
purely a VCF export artifact: a VCF allele is a flat string, so a cyclic traversal
must be linearized. And a duplication CANNOT be re-expressed as a plain insertion
without loss, because the inserted copy carries its own variation and **VCF positions
are reference coordinates** — there is nowhere to place a variant that lives inside
inserted sequence. This is the limitation pangenome graphs exist to remove.

### Clarification 2026-09-22 — CONFLICT is mostly NOT caused by revisiting
Dylan asked whether the two-alleles-at-a-child-variant situation arises whenever a
repeat revisits a bubble. Measured on chr21:

| | count |
|---|---|
| records whose own AT contains a revisiting traversal | 37 (0.011%) |
| records carrying a CONFLICT tag | 1,855 (0.54%) |
| both | 11 |
| revisit but NO conflict | 26 |
| **conflict but NO revisit** | **1,844** |

CONFLICT names all 45 samples, roughly 200-470 records each — far too widespread to
be explained by tandem duplications, of which only 37 records chr21-wide show any
sign. The dominant cause is therefore something else; the most plausible candidate
is assembly fragmentation (each HPRC haplotype is many contigs, and two contigs of
one sample traversing the same region give that sample multiple paths), but **this
cannot be determined from the VCF alone and has not been verified.**

⚠ **Limitation of that test.** A record's `AT` lists the traversals of ITS OWN
bubble, so it cannot reveal that a sample passed through that bubble twice — that
is only visible at the enclosing parent. So the cross-tabulation UNDERSTATES the
link for child records inside a duplicated region. The mechanism is real; its
genome-wide share is not what the 11/37 overlap suggests.

Note the asymmetry it exposes: the 257 kb parent HAS a revisiting allele but is NOT
flagged CONFLICT, because at the parent level HG01891 has exactly one allele (58)
whose walk happens to loop. The conflict only materialises at CHILD records, where
that single parent allele corresponds to two passes that may disagree.

**Design implication for `T`.** A snarl-tree segmentation would make this one region
ONE block instead of 370 overlapping ones, and blocks disjoint by construction. The
usual objection — that snarls are not LD blocks — is weak for OUR design, because in
a Li-Stephens copying model the linkage is carried by the donor panel, not by the
block boundaries; and switching donors at snarl boundaries is graph-valid by
construction, which is what `output_support` needs anyway.


## 2026-09-18 — A-EVL-readmapping: reproduced the published read-mapping BASELINE, close but not exact

**Goal.** The sharpest available test of our installation. The paper's read-mapping
experiment is chr21-only and uses five EXTERNAL donors, so scope cannot excuse a
miss, and the BASELINE (unobfuscated graph) needs no obfuscation sweep.

**Method.** Extracted chr21 reads for the paper's five donors (HG00138 EUR,
HG00635 EAS, HG01112 AMR, HG02698 SAS, NA18853 AFR) from the 1000G 30x CRAMs;
built the graph per `vg_prep.py` and aligned single-end per `quick_align.py` with
vg 1.68.0; scored with `vg stats -a` using the metric definitions in
`external/PanMixer/plot_util/plot_helpers.py` (each of perfect / gapless / MAPQ60 over total aligned).

**Result — the right regime, but a systematic offset.**

| metric | published baseline | ours | difference |
|---|---|---|---|
| perfect | 77.83% | 79.05% | +1.22 |
| gapless | 95.62% | 95.96% | +0.34 |
| MAPQ60 | 77.01% | 79.23% | +2.22 |

Supporting signals that the pipeline itself is sound: ~10.9-11.4M reads aligned per
donor, matching the scale of the authors' own sample output preserved in
`external/PanMixer/tools/downstream/utility_out/giraffe_parser.py` comments (total
aligned 11,521,845; that sample works out to 77.68% perfect / 95.82% gapless, close
to the published figures, so those comments come from this pipeline). That sample
also shows `Total paired: 0`, confirming the authors mapped SINGLE-END as we did.
The per-donor ordering is biologically sensible: NA18853 (AFR) is lowest on every
metric, as expected for an African genome against a reference-biased graph.

**Provenance.** chr21 graph VCF 287,147 SNP records · vg 1.68.0 · 5 donors ·
`logs/pm_vg.log`, `logs/pm_vg2.log` · read slices in
`external/PanMixer/read_fastqs/`.

### Correction 6 — my explanation for the offset was wrong, and the test refuted it
[INCORRECT] - the offset is likely caused by the reference build: I used chr21 from the 1000G GRCh38 analysis set, which hard-masks 2.05 Mb more of chr21 (8,671,409 N vs 6,621,364 N) than the UCSC hg38 chr21 that clean_fa.py implies; masking the false-duplication regions removes mapping ambiguity and should inflate perfect alignment and especially MAPQ60
[CORRECTION - 2026-09-18]: Tested by rebuilding the entire graph and re-aligning all five donors against UCSC hg38 chr21. The result is essentially IDENTICAL: 79.05 / 95.96 / 79.23 (UCSC) versus 79.03 / 95.96 / 79.26 (analysis set) — a shift of at most 0.03 points, against an offset of 1.2-2.2. The hypothesis is refuted. The reason is clear in hindsight and should have been predicted: the reads were extracted from CRAMs that were themselves aligned to the analysis set, so no reads originate in the extra-masked regions under EITHER graph. The masking difference was real but could not matter. **The measured N-count difference stands; only the causal claim was wrong.**

### Hypotheses tested, and what remains unexplained
1. **Reference build masking** — REFUTED by direct rebuild (above).
2. **`convert_2_vcf` genotype rewriting** — MEASURED, at most a minor contributor.
   That converter rewrites any sample field lacking `|` to `.|.`
   (`external/PanMixer/tools/common/convert_2_vcf.py`:71-73), which would shrink allele support and
   drop variants at the `-c 1` step, lowering mapping rates. On 200,000 chr21
   records, 2.38% of sample fields lack `|`, but **zero records lose all support**,
   so few if any variants would be dropped. Direction is right; magnitude is not.
3. **`vg stats` version** — UNTESTABLE, and the attempt produced a finding of its
   own: the published Methods say alignment used Giraffe v1.68.0 and the GAMs were
   analysed with `vg stats` **v1.36.0**. Those two are incompatible. Running the
   real v1.36.0 binary on a v1.68.0 GAM aborts with "obsolete, invalid, or corrupt
   protobuf input". The stated combination cannot have been used as written.

[INCORRECT] - **VERDICT.** The installation runs the pipeline correctly and lands in the right regime with the right internal structure, but this is **NOT an exact reproduction** of the published baseline, and the residual +1.2 to +2.2 point offset is **unexplained**. Do not describe the read-mapping baseline as reproduced.
[CORRECTION - 2026-09-26]: Verdict unchanged; the range was misstated. The table directly above gives +1.22 (perfect), +0.34 (gapless) and +2.22 (MAPQ60), so the offset spans **+0.34 to +2.22** points. "+1.2 to +2.2" dropped the gapless column. Caught by the 2026-09-26 primer audit.

**What this does NOT establish.** Nothing about whether the paper's numbers are
wrong. An offset of this size is consistent with an undocumented difference in read
extraction or filtering, in the exact graph VCF the baseline was built from, or in
the tool versions actually used — all invisible from the paper and the released code.


## 2026-09-18 — A-EVL-published-validation: checked our run against the PUBLISHED paper

**Goal.** Dylan supplied the PUBLISHED PanMixer paper (Nature Communications,
DOI 10.1038/s41467-026-77591-0, "Article in Press", 12 pages). All prior analysis
in this project was made against the bioRxiv preprint. Validate our installation
against what the authors actually report, and re-check our findings against the
version of record.

**Method.** Extracted the published text, diffed it against the preprint, and ran
a 9-agent workflow (4 readers + adversarial verifiers + synthesis) over both
versions and the code. Then hand-verified every load-bearing claim.

**VERDICT: the installation is running correctly**, in the precise sense that it
faithfully executes the released code at `c182c38` and produces internally
coherent results. Four independent signals, none guaranteed in advance:

| check | published | ours | verdict |
|---|---|---|---|
| cohort size | 44 individuals | 44 / 88 haplotypes | match |
| capacity -> normalized utility loss | budget is a fraction of max loss | 928.92/9288.85 = 10.000% at capacity 0.1 | exact |
| per-subject cost, longest chromosome | 5 min / 14 GB (Table 1) | 68 s / 1.7 GB on chr21; x5.47 to chr2 = 6.2 min / 9.3 GB | consistent |
| determinism / invariants | not reported | byte-identical per seed, 0 allele-bound violations | pass |

**Provenance.** `paper/s41467-026-77591-0_reference.pdf` (sha256 `2bbba736...`) ·
`external/PanMixer` @ `c182c38` · chr21 on the 30x GRCh38 panel · HG00438 ·
workflow `wf_c1f52889-672`, 9 agents, 0 errors.

**What the verdict does NOT say.** We have reproduced NO reported result. Every
headline number — `eps_private` = 0.002, average utility loss 0.28, AF divergence
0.006/0.006/0.005, LD 0.003, read mapping 77.82/95.61/77.01 — is unrun. "Correct
installation" means "we run their code as shipped", NOT "the code implements the
published method"; we have four confirmed places where it does not.

### Correction 4 — the "47 individuals" discrepancy was preprint-only
[INCORRECT] - The cohort is **44 individuals / 88 haplotypes**, not the paper's "47" (the HPRC VCF carries 45 samples, one being `chm13`, and the pipeline drops it and chrX).
[CORRECTION - 2026-09-18]: The PREPRINT says 47; the PUBLISHED paper says 44 (page 1 line 20, and "n = 44" in the Fig. 2 caption and Table 1). There is no discrepancy against the version of record — the published cohort size matches our measurement exactly. Recorded because the original line reads as though the authors' paper disagreed with their own data, which is true only of the superseded preprint.

### Correction 5 — the HMM's anchor set is the PanGenie overlap, not the 1000G overlap
[INCORRECT] - Anchors are pangenome records that exactly match a PanGenie "bi_all" record ... the HMM-vs-AF routing is therefore exposed to the GRCh37/GRCh38 build mismatch along with everything else derived from the 1000G match set
[CORRECTION - 2026-09-18]: The first half is right and the implication is wrong. `hmm.py:32-33` loads `pangenome_to_thousand_g_alignments.pickle`, and in `get_mappings.py:28` the variable `thousand_g_alignments` is bound to **`PG.vcf.gz`** — the PanGenie callset — NOT to 1000G. The name is a trap. Because PanGenie is GRCh38, the per-block HMM-vs-AF ROUTING is NOT degraded by the build mismatch. What the mismatch does corrupt is the LD block BOUNDARIES (plink runs on the 1000G panel), the attack database, and the `af_loss`/`ld_loss` site sets. Two different anchor notions; do not conflate them. The paper meanwhile defines anchors as V_1000G (page 8 lines 944-945), which the code does not implement.

### The build inconsistency is in the PUBLISHED paper, not only the released code
The published Methods pin V_1000G to "the Phase 3 phased 1000 Genomes dataset" as
a SUBSET of V_0, the GRCh38 backbone (page 8 lines 941-945), and the NEW Data
availability section pins the exact GRCh37 EBI URL, release 20130502 (page 11).
The backbone is stated as GRCh38 three times. The full text contains zero
mentions of GRCh37, liftover, 30x, high-coverage or NYGC. Our 1,025 vs 258,610
chr21 match measurement shows the stated containment cannot hold as written. The
preprint had no Data availability section, so the published version sharpens the
inconsistency rather than resolving it. Note also that the official GRCh38
liftover of Phase 3 was WITHDRAWN in 2021, with 1000 Genomes redirecting users to
the 30x GRCh38 collection — the panel we substituted, and the one the repo's own
unwired `get_blocks_grch38.sh` downloads.

### Numbers that CHANGED between preprint and published — use the published ones
- `eps_private`: 0.001 -> **0.002**. Its DEFINITION also changed, from "the
  smallest privacy risk value at which all targets could no longer be linked" to
  "the **largest**" — the preprint's wording was backwards.
- Wasserstein AF divergence at `eps_private`: 0.004 / 0.004 / 0.002 ->
  **0.006 / 0.006 / 0.005** (all alleles / SNPs / rare SNPs).
- The published version adds a membership-inference analysis absent from the
  preprint (0 mentions in the preprint, 4 in the published version).
Any reproduction target must come from the published version.


## 2026-09-18 — A-DAT-grch38-repin: chr21 rebuilt on GRCh38; PanMixer executed and measured

**Goal.** Make PanMixer runnable here, on build-correct data, and replace the
inferred claims from the ingest with measurements on real chr21.

**Method.** Wrote a blocking local `sbatch` shim (`scripts/sbatch`) so the pinned
checkout runs UNMODIFIED, and re-pointed the hardcoded `1000g_phased.vcf.gz` slot
at the 30x GRCh38 panel. Ran the shipped preprocessing for chr21, then
`obfuscate.py` for HG00438 at two capacities and two seeds, then diagnostics.

**Result — the repin works.** `get_mappings`, PanMixer's own matching code, reports
**258,610 strict / 268,851 relaxed** matches, against 1,025 / 2,176 on the wired
GRCh37 panel. That is exactly the number an independent `bcftools` comparison
predicted, so two independent routes agree. The attack database is now
(3202, 258610, 2) instead of one built on 1,025 spurious matches.

**Result — PanMixer runs and is well-behaved.** 68 s per subject-chromosome, peak
RSS 1.7 GB. Same seed gives byte-identical `new_haplotypes.npy` and `xsol.npy`;
different seeds differ. Output shape (340824, 2), **zero** allele-bound violations,
VCF correctly bgzipped and indexed. Capacity is respected: 0.1 -> 10.0% utility
loss. At capacity 0.5 it spends only 33.1%, i.e. the knapsack saturates — all
positive-value moves fit — so capacity is not a tight knob at the top end.

**Result — block structure.** 14,137 plink LD blocks but 100,757 entries in
`blocks_dict.json`, **88.4% singletons**. Only **10.4% of blocks** take the HMM
path, though those hold **71.4% of variants**. [INCORRECT] - For us: `T_blocks` ~ 100,757 on
chr21 with `K_states` = 88, so FFBS is ~9M operations.
[CORRECTION - 2026-09-28]: The block structure is correct; reading it as OUR chain length is
not. Model C was adopted 2026-09-28 — one chain position per VARIANT — so T is 340,824 raw
and 306,480 after the `chain_span` cut, roughly 3x the figure above, and FFBS is ~30M
[CORRECTION - 2026-09-29]: 306,480 was the ONE-SIDED cut (p arm only). The telomeric end
was closed on 2026-09-29 for the same reason — 592 variants above the map end formed a
593-position zero-distance run — so the figure is **305,887** (340,824 - 34,345 - 592),
89.75%. The FFBS cost is unchanged in magnitude.
operations. `K_states` = 88 is right and stands. LD blocks play no role in our mechanism at
all under C.

**Provenance.** `external/PanMixer` @ `c182c38` · chr21 rebuilt on the 30x GRCh38
panel (sha256 `a925c112...`) · HG00438 · capacity 0.1 and 0.5 · seeds 123, 456 ·
`logs/pm_chr21.log`, `logs/pm_obf.log`.

### Confirmed on real data — the sampler is inert
300 of 300 randomly chosen HMM-path blocks returned a block **identical to a single
donor haplotype**. The synthetic-panel result from the ingest holds on real chr21.

### Correction 3 — I overstated the eps_pmi defect, and the measurement refutes it
[INCORRECT] - the forward panel's self-inclusion (obfuscate.py:130,132) has flattened the privacy score, so eps_j is pinned near log(2N) = 4.4773 nats for essentially every HMM block regardless of true rarity
[CORRECTION - 2026-09-18]: Measured on 150 real HMM-path blocks for HG00438: eps_j runs 0.071 to 4.508 with median 1.242, and only 9/150 sit within 0.05 nats of log(88) = 4.4773. It is NOT pinned. (The value 4.4773 itself remains correct — it is log(2N) and is still the ceiling; what was retracted is the claim that eps_j SITS there.) What self-inclusion actually does is impose a CEILING at ~log(2N): because the target always matches itself, p(h_bj) >= 1/2N, so eps_j <= ~4.48 however rare the block truly is. Re-scoring the same blocks with the target REMOVED from the emission panel gives median 1.299 but max 78.1. For the 9 ceiling-bound blocks the leave-one-out value has median 12.06 and max 78.12 — an understatement of up to 17.4x. So the defect is real but its shape is the opposite of "flattened": typical blocks are almost unaffected (median ratio 1.03x) while the RAREST blocks — the most identifying, exactly the ones a privacy mechanism must prioritise — are truncated hardest. The knapsack therefore systematically under-values the highest-risk blocks.

**What this does NOT establish.** Still nothing about whether the PAPER's numbers
are affected: eps_j enters the LP only through a ranking, and we have not measured
how much the selection changes when the ceiling is removed. Do not claim their
frontier is wrong; claim only that the shipped scorer truncates rare blocks.


## 2026-09-18 — A-DAT-panmixer-ingest: PanMixer pinned, its data acquired, and its mechanism mapped

**Goal.** Bring the upstream dependency in, pin it, run it, and establish how much
of it RanPanMixer can build on — the proposal names PanMixer as the source of its
target-independent prior, its blocks, and its utility information.

**Method.** Cloned G2Lab/PanMixer, pinned `c182c38d5bc8bb6f00f4b0b101207c4a009ca045`
(the only release, 2026-06-06, MIT). Relocated the checkout to `/data` because
PanMixer writes inside its own tree. Built its `environment.yaml` env. Downloaded
its three input datasets (16 GB) and its two external tools. Fetched the paper
(bioRxiv DOI 10.64898/2026.02.16.706152, CC-BY-NC-ND 4.0). Then ran a 22-agent
workflow: 8 readers over the paper and six code subsystems, each followed by an
independent adversarial verifier, plus 4 checks of what OUR proposal asserts about
PanMixer, a synthesis and a completeness critic. Finally re-verified every
load-bearing claim by hand — reading the cited lines, executing the shipped HMM
class on a synthetic panel, and measuring cross-build variant overlap directly.

[INCORRECT] - **Result.** The cohort is **44 individuals / 88 haplotypes**, not the paper's "47" (the HPRC VCF carries 45 samples, one being `chm13`, and the pipeline drops it and chrX).
[CORRECTION - 2026-09-26]: 44 / 88 is right and matches the PUBLISHED paper exactly; "47" is the superseded preprint only (Correction 4). Marked in place on 2026-09-26; the copy under Correction 4 was made without marking this original.
Two findings dominate, both independently reproduced:

1. **The released preprocessing pipeline matches GRCh38 against GRCh37.** Measured
   on chr21: 1,025 of 340,824 pangenome records match the wired 1000G Phase 3 panel
   exactly on (POS,REF,ALT) — 0.30%. The repo ships an unwired GRCh38 path.
2. **The Li-Stephens sampler never recombines.** Four compounding defects, present
   identically in the standalone tool AND in the pipeline copy the paper's
   precomputation imports. Executing the shipped class: total switch mass 5.8e-12
   at a 1-index gap; the sampled block equals exactly one donor haplotype. With the
   paper's own Fig. 6 constants it would be 0.13-0.99.

Both are in `docs/bugs.md` with evidence. The consequence for us is the headline:
**"reuse PanMixer's cohort HMM as R_D" is not available.** What we inherit is its
DATA, its BLOCK structure, its forward recursion's algebra, and its EVALUATION
suite — roughly 55-65% of ~8,660 lines by volume, but under 10% of the mechanism.

**Provenance.** `external/PanMixer` @ `c182c38` · HPRC v1.0 PGGB GRCh38 VCF
(sha256 `ead25541...`) · 1000G Phase 3 chr21 (sha256 `1942e070...`) · PanGenie
Zenodo 7669083 (sha256 `e8c0c48d...`) · paper sha256 `99ac1bd2...` · workflow run
`wf_a30688a6-3d7`, 22 agents, 0 errors.

**What this does NOT establish.** Nothing about RanPanMixer's own mechanism — no
line of it exists yet. It also does not establish that the PanMixer PAPER's
published numbers are wrong: the defects are in the RELEASED code, and the unwired
GRCh38 script plus the code branches preferring 3,202-sample artifacts suggest the
authors' own runs used a path the release does not wire in. We have not reproduced
their figures, so we cannot say which code produced them. **Do not state or imply
that their results are invalid.** What we CAN say is narrower and sufficient: the
released pipeline as documented does not reproduce the paper, and its sampler as
shipped is not the model the paper describes.

### Control — the build mismatch is measured, not inferred

Downloaded the GRCh38 30x panel (the one the repo's unwired `get_blocks_grch38.sh`
fetches) and repeated the identical comparison against the same pangenome chr21
records. Same pangenome, same method, only the panel's build differs.

| metric | Phase 3 (wired, GRCh37) | 30x (unwired, GRCh38) |
|---|---|---|
| exact (POS,REF,ALT) matches | 1,025 | 258,610 |
| relaxed (POS,REF) matches | 2,176 | 268,851 |
| match rate of 340,824 pangenome records | 0.30% | 75.88% |

252x. The deficit is the reference build, not the data. The residual ~24% is
expected — the pangenome carries graph-only and structural variants absent from a
SNV/INDEL/SV panel — but has not been characterised and should not be assumed benign.

**Still not established:** which code produced the paper's numbers. The unwired
GRCh38 script and the branches preferring 3,202-sample artifacts are consistent
with the authors having run a corrected path, but that is an inference, not a
measurement. Do not write that PanMixer's results are invalid.

### Correction 1 — the proposal's characterization of PanMixer is wrong in four places
These were checked against paper and code by dedicated agents and by hand.

[INCORRECT] - The HMM implementation by PanMixer uses a transition model with one probability for remaining in the same donor state and another probability for switching donors, and the per-target running time is O(TK).
[CORRECTION - 2026-09-18]: The stay/switch FORM is real and its forward recursion is genuinely O(K) per step (obfuscate.py:140-152). But PanMixer's SAMPLER is O(A*K^2) in time and memory — `get_transition_probabilities_without_subject` materialises every dense K'xK' matrix for a block up front (:117-124). More importantly the chain runs WITHIN a block over anchor SNPs, with a fresh uniform start state per block and NO transition model between blocks; there is no T-blocks x K-states chain, no backward pass and no FFBS anywhere in the repo. Our across-block chain and backward sampler are new work, not inherited.

[INCORRECT] - PanMixer protects individuals contributing to the pangenome by selecting and applying block-level haplotype obfuscations. The proposed mechanism protects an external target individual's mapped path without modifying the graph.
[CORRECTION - 2026-09-18]: "Without modifying the graph" is not the contrast. PanMixer does not modify graph topology either — its release is the same cohort VCF with one sample column rewritten. The real contrast is WHO is protected (a cohort member vs an external individual) and WHETHER the release is written back into the reference.

[INCORRECT] - PanMixer-compatible LD blocks or graph subpaths and the cohort HMM can serve as the target-independent prior R_D.
[CORRECTION - 2026-09-18]: Only the block BOUNDARIES (plink on an external panel) and the transition FORM are target-independent as shipped. The allele-frequency table includes the target at PanGenie-matched and novel sites (get_af.py:41-53,68-73); the support counts include the target (obfuscate.py:208); the forward/PMI panel includes the target (:130,132); and the site-to-block assignment ranges over sites private to the target. Every one must be recomputed leave-one-out before R_D satisfies the theorem's hypotheses.
[CORRECTION - 2026-09-28]: The catalogue above is accurate about PANMIXER and the final sentence is wrong about US. It reads as a precondition on our R_D; it is not. Our target is EXTERNAL to G, D, every panel and the attack database, so none of those artifacts can contain it and there is nothing to recompute: K = 88, and constraint 1 holds by construction. Leave-one-out is a requirement of the HEAD-TO-HEAD ARM ALONE, where a cohort member stands in as target and externality must be simulated. This entry is directly contradicted by the 2026-09-25 entry in this same file ("In v1 the target is a PARAMETER: it never enters the matrix, `support_D`, or the state set"). See CLAUDE.md constraint 1 and `loo_cohort` in docs/terms.md, both corrected in commit 1f8ea3f.

[INCORRECT] - PanMixer samples a candidate constrained to differ from the target's own block.
[CORRECTION - 2026-09-18]: The PAPER says this (p.10 lines 653-654); the released CODE implements no such constraint — `sample_block_prior` is an unconditional prior draw with no comparison to the original and no rejection loop. So its candidate generator is already target-independent given the panel. PanMixer's target-dependence comes instead from the knapsack SELECTION and from releasing unselected blocks verbatim. Note the consequence: because a no-op resample is possible while eps_j stays positive, the LP can buy privacy credit at zero utility cost from moves that change nothing.

### Correction 2 — an error this session made and caught
[INCORRECT] - the shipped transition code gives a total switch mass of 1.4e-12 at a 1-index gap, quoted as the sampler's behaviour
[CORRECTION - 2026-09-18]: That figure belongs to `get_transitions` (the FORWARD pass, no factor 4). The SAMPLER uses `get_transition_matrix_without_subject`, which includes `* 4 *` and gives 5.8e-12. Measured directly. The degeneracy conclusion is unchanged, but the number was attributed to the wrong function. Relatedly, the first draft claimed the forward pass "omits" a factor 4 that the sampler "includes", implying 4 is correct; paper Fig. 6 gives `d = Delta_x * Ne * r` with NO factor 4, so the SAMPLER is the deviation.

### Decision — the integration seam is `new_haplotypes.npy`
Every PanMixer attack and utility evaluator consumes one artifact: a
`(n_sites, 2)` int16 array, `-1` = missing, row-aligned to the chromosome's VCF.
Our sampler will emit exactly that file, making it effectively a third stacker
strategy, so their entire benchmark runs on our output unchanged. This is the
single highest-leverage design decision from this session and it removes most
implementation-difference confounds from the head-to-head.

### Decision — the head-to-head is leave-one-out, aligned on measured attack success
`eps_pmi` and `tau` are formally incomparable, and so are `capacity` and `tau`.
The comparison aligns the arms on MEASURED empirical attack success and reports two
frontiers. A hard limit to state in any write-up: the PGGB graph topology was built
from all HPRC assemblies, so even a full leave-one-out only APPROXIMATES our
external-target threat model on their data.

### Finding — their cohort-level utility metrics cannot see target fidelity
A pure-prior draw (our tau = 0) is a real cohort haplotype: it would score
near-perfectly on U_AF, U_LD and external read mapping while carrying zero
information about the target. Any comparison must carry a `target_fidelity` axis or
our maximal-privacy point looks like a free lunch. The critic established that such
a measurement already exists and is mechanism-agnostic
(`MIA_privacy.py:80-83,140`), correcting the synthesis's claim that it had to be built.

