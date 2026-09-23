# memory.md — the linear ledger

**DO NOT read this file whole.** `grep` it.

```
grep -n "^## 20" docs/memory.md          # list entries
grep -n "INCORRECT\|CORRECTION" docs/memory.md   # what we got wrong
```

`docs/plan.md` exists so a new session does not have to come here.

## Rules

1. **Append only.** Never delete, never overwrite. Newest at the bottom.
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
2. **The most identifying sites score ZERO.** Every allele here has population
   frequency 0 in any external callset, so `-log(0)` -> infinity -> zeroed by the
   `#FIX ME` in `get_support_and_pmi.py`. The sites where each haplotype is unique
   — a perfect fingerprint — contribute nothing to PanMixer's privacy score.
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

**VERDICT.** The installation runs the pipeline correctly and lands in the right
regime with the right internal structure, but this is **NOT an exact reproduction**
of the published baseline, and the residual +1.2 to +2.2 point offset is
**unexplained**. Do not describe the read-mapping baseline as reproduced.

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
path, though those hold **71.4% of variants**. For us: `T_blocks` ~ 100,757 on
chr21 with `K_states` = 88, so FFBS is ~9M operations.

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

**Result.** The cohort is **44 individuals / 88 haplotypes**, not the paper's "47"
(the HPRC VCF carries 45 samples, one being `chm13`, and the pipeline drops it and
chrX). Two findings dominate, both independently reproduced:

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

