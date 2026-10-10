# bugs.md — symptom -> proven fix

**DO NOT read on startup.** `grep` by symptom.

Different retrieval key from `docs/memory.md`: results are looked up by DATE, bugs
by SYMPTOM. Folding them together scatters one symptom across years.

⚠ The highest-severity class is **SILENT DEGRADATION** — a resource missing, a
gate defaulting to its passing value, a flag accepted and ignored, an extractor
returning empty instead of raising. These produce NUMBERS, not errors. When you
fix one, record **what the wrong output looked like**, because that is what the
next person will recognise.

Every entry below has been fixed and verified at least once. Do not re-derive.

---

## Tooling / reading source documents

### `pypdf` returns empty strings for every page of a PDF
- **[2026-08-29] Image-only PDF extracts as blank, with no error**
  **[Symptom]** `pypdf.PdfReader(f).pages[i].extract_text()` returns `""` or a
  couple of zero-width spaces for all 7 pages. No exception. The file is 5.3 MB, so
  "empty" is obviously wrong, but nothing in the return value says so.
  **[Cause]** `archive_docs/The_Six_File_Lab_Record.pdf` is a rendered HTML page —
  the content is a raster image with no text layer. `pdftotext`/`poppler-utils` and
  `pdftoppm` are not installed on this machine, so the Read tool's PDF path also
  fails, with a clearer message.
  **[Proven fix]** Render to PNG with PyMuPDF (installed: 1.27.2) and read the
  images: `fitz.open(f)[i].get_pixmap(dpi=150).save(out)`. 150 dpi is legible for
  body text; use `dpi=300` with a `clip` rect for equations.
  **[Severity]** Silent. An agent that trusted the empty extraction would conclude
  the document was blank and proceed to invent its contents.

### A PDF's text layer disagrees with what the page shows
- **[2026-08-29] Superscripts and inequality operators lost in extraction**
  **[Symptom]** Extracted text for the proposal's Theorem 1 reads
  `e−ετ ≤ Qτp(y)/Qτq(y) ετ` — a malformed inequality. It is not obvious whether
  the document or the extractor is at fault.
  **[Cause]** Both are possible, and they need different responses. Here the
  rendered page shows the same defect, so it is genuinely in the document
  (recorded in `docs/memory.md`); but elsewhere in the same file `e−ητ ≤ Zq/Zp ≤ eητ`
  extracted correctly, proving the extractor handles the construct in general.
  **[Proven fix]** Never conclude a formula is wrong from extracted text. Render
  the region at 300 dpi and read the image before recording anything about it.
  **[Severity]** Would have produced a false claim about someone's manuscript.

## The docs verifier itself

### A contract check reports SKIPPED when it should report FAILED
- **[2026-08-29] Retired-term check skipped because of an em-dash**
  **[Symptom]** `test_no_retired_term_is_a_live_endpoint` printed
  `SKIPPED [1] ... nothing retired yet`, while `docs/terms.md` had a populated
  `## Retired` section. The suite was green. Nothing was red, and nothing was
  checked.
  **[Cause]** The extraction regex matched `-` and `--` but not the en/em dashes
  the prose actually uses (`\u2013`, `\u2014`). Zero terms extracted, so the
  guard clause took the "nothing retired yet" branch.
  **[Proven fix]** `[-\u2013\u2014]+` in the pattern. Verified by injecting a
  retired term into `docs/plan.md` and confirming the check fails.
  **[Severity]** SILENT. A skip is not a pass, but a green summary line reads
  like one. **Rule:** any check with a `pytest.skip` guard must be armed by
  deliberately breaking what it protects, or the guard is untested and the check
  may never have run at all.

### The path-existence check flags prose as a missing path
- **[2026-08-29] `n/a` reported as a nonexistent path**
  **[Symptom]** `test_every_referenced_path_exists` failed on
  ``CLAUDE.md: `n/a` `` — the reporting-format instruction telling authors to
  write `n/a` in an empty table cell.
  **[Cause]** The `<dir>/<rest>` alternative in `_PATHY` accepted a one-character
  first segment, so `n/a` parsed as a directory reference.
  **[Proven fix]** Require a two-character-or-longer first segment in the
  directory alternative of the pattern. (Written out rather than quoted as a
  regex: a backticked brace quantifier trips the brace-shorthand check, which
  cannot tell a glob from a repetition count.)
  **[Severity]** Noise, not corruption — but a verifier that cries wolf gets
  suppressed, and a suppressed verifier is the framework's whole failure mode.

## Our own run discipline

### A stale status sentinel reports a job that is still running as FAILED
- **[2026-09-18] Re-running a job without clearing its sentinel**
  **[Symptom]** A background `until [ -f run.status ]; do sleep; done` returned
  immediately with `exit=1`, so a still-running `get_af` was reported as failed. I
  then started a SECOND copy of the same job to "debug" it, doubling a 13 GB-RSS
  process and the I/O behind it.
  **[Cause]** The first attempt (which really did fail, on the dateutil import)
  had already written that run's .status sentinel under logs/. The re-run reused
  the same path,
  so the sentinel existed before the job finished. `[ -f status ]` cannot tell
  "finished" from "left over".
  **[Proven fix]** `rm -f <sentinel>` at launch, in the same command that starts
  the job. When diagnosing, check the sentinel's mtime against the job's start
  time before believing it — that is what exposed this (status 13:40:31 predated
  the run's own log at 13:42:01).
  **[Severity]** SILENT, and the same class as the upstream defects in this file:
  the wrong answer arrived as a plausible value, not an error. Existence is not a
  completion signal; a sentinel must be newer than the launch that it reports on.

### Our GRCh38 repin put the target inside its own allele frequencies
- **[2026-09-24] The 30x panel contains 39 of the 44 HPRC donors**
  **[Symptom]** None visible — the frequencies look normal.
  **[Cause]** `get_af.py` branch 1 draws counts from the 1000G panel alone, and its
  comment asserts "None of the pangenome subjects appear in 1000g_phased". True for
  Phase 3 (0 of 44 overlap) but **false for the 30x panel we substituted: 39 of 44
  HPRC donors are present, including HG00438**. So the target's own alleles are
  counted in the frequencies used to score it.
  **[Proven fix]** A `loo_cohort` must drop the target from the **1000G panel** as
  well as from the pangenome. Recorded as a requirement of the comparison design.
  **[Side effect worth knowing]** This also suppresses the `#FIX ME` zeroing: with
  the target in the panel its allele is always counted, so f_v is never 0 and that
  branch never executes. Under the shipped Phase 3 panel it would.
  [INCORRECT] - **[Severity]** HIGH for us — it is a constraint-1 violation introduced by our own fix.
  **[CORRECTION - 2026-09-28]** **[Severity]** HIGH for the COMPARISON ARM; **not applicable
  to our mechanism.** The measurement is right — 39 of 44 HPRC donors including HG00438 are
  in the 30x panel — and it matters wherever a cohort member stands in as target. It is not
  a constraint-1 violation of our design, whose target is external to every panel by
  construction. Doubly moot since 2026-09-28: under model C the mechanism uses no external
  panel and computes no allele-frequency table at all. ⚠ The overlap still must be removed
  from the ATTACK DATABASE, and the obvious fix is incomplete: all 39 overlaps sit in the
  698 samples the 30x release ADDED beyond Phase 3 to complete trios, so dropping the 39
  named samples leaves their first-degree relatives in the database.

### Killing a duplicate process
- **[2026-09-18] Never `pkill -f` a pattern that matches your own command line**
  **[Proven fix]** Resolve the exact PID first (`ps -o pid=,etimes=,cmd= -C python`,
  pick by elapsed time) and `kill` that PID. The younger of two identical jobs is
  the accidental one when the older is the pipeline's.
  **[Severity]** LOW here, but `pkill -f get_af` would have matched the shell
  running the check as well.

### grep silently returns NOTHING on a file containing a NUL byte
- **[2026-09-18] Every grep over the extracted preprint text failed silently**
  **[Symptom]** `grep -c '47' panmixer_paper.txt` printed nothing and exited 1,
  while `python: open().read().count('47')` returned 5. Earlier greps over the same
  file for section headings had returned empty and I wrongly attributed that to
  formatting.
  **[Cause]** PyMuPDF extraction of the preprint emitted 2 NUL bytes from a math
  glyph in the PMI formula (`PMI\x00hbj, ...`). GNU grep classifies a file with a
  NUL as BINARY and, rather than erroring visibly, prints nothing for `-o`/`-c`
  and exits 1 — indistinguishable from "no matches".
  **[Proven fix]** Strip NULs at extraction time
  (`data.replace(b'\x00', b' ')`), or pass `grep -a`. **Rule: never conclude a
  string is ABSENT from extracted PDF text on the strength of a silent grep.**
  Confirm absence with a python `.count()` before recording it.
  **[Severity]** SILENT and high-impact: it produces "not present" conclusions
  about a source document, which is exactly how a false claim about someone else's
  paper gets recorded as fact.
  **[Related, different cause]** A grep pattern of the form `.\{80\}TERM.\{80\}`
  also returns nothing on PDF-extracted text, because the extraction wraps lines
  short and the required context does not fit on one line. Same empty output, a
  completely different reason. Check the pattern before blaming the file.

### An audit reports a true number as fabricated
- **[2026-09-26] "not found in `docs/`" was read as "invented"**
  **[Symptom]** A checker greps `docs/` for a number quoted in teaching material,
  finds nothing, and recommends deleting it — or worse, writes a published retraction
  of it. Four primer numbers were flagged this way (274,458; 177,865; 3,077.60; 8,661).
  One was already "retracted" in a draft, with an arithmetic argument that was itself
  wrong (it treated a subset count as if it had to equal its superset).
  **[Cause]** The numbers were measured during a session and quoted directly, never
  logged. Absence from `docs/` is a provenance gap, not evidence of falsehood.
  **[Proven fix]** Re-measure from the artifact BEFORE deleting anything (the mask
  file, the run log, the checkout). All four were exact. Log them, then quote them.
  And the upstream fix: anything measured in a session and quoted anywhere gets a
  `docs/memory.md` line the same day. See `docs/memory.md` 2026-09-26.

## PanMixer (upstream, commit c182c38)

Defects in the DEPENDENCY, not in our code. Recorded here because we will hit each
one again, and because several produce plausible NUMBERS rather than errors — the
silent-degradation class. "Proven fix" = what we do about it. Every entry was
verified by reading the cited line or by executing the shipped code, on 2026-09-18.

### The privacy score ignores most of what is in the block
- **[2026-09-25] eps_j scores only anchor variants, so rare non-anchor alleles are free**
  **[Symptom]** A block containing a genuinely rare allele scores as unremarkable.
  Worked case: block 1083 holds two common anchors (AF 0.934) and one non-anchor at
  AF 0.00385. The allele-frequency branch would give 5.6956 nats; the code returns
  **0.0349** — a 163x difference.
  **[Cause]** `forward_algorithm` builds `observed` and the emission panel from
  `anchor_indices` only, so non-anchor variants are not positions in the chain and
  contribute nothing to the likelihood.
  **[Measured]** Over 400 random HMM-branch blocks: median eps_j 1.079 while the
  all-variant allele-frequency sum is 2.769; **in 32% of blocks the ignored
  non-anchor rarity exceeds eps_j itself**.
  **[Proven fix]** Do not treat eps_j as "the self-information of the block" — it is
  the self-information of the block's ANCHOR alleles. For v1 we drop anchors entirely
  and make every variant a chain position (see `A-IMP-v1-sampler` in docs/plan.md).
  **[Severity]** HIGH, and silent: the number is plausible and monotone, just blind
  to the rarest material.

### Scoring and sampling disagree about which blocks are HMM blocks
- **[2026-09-25] The two branch thresholds differ by one anchor**
  **[Symptom]** A block gets an HMM-derived privacy score for a replacement that was
  drawn from allele frequencies.
  **[Cause]** Scoring takes the HMM when anchors >= 1
  (`get_support_and_pmi.py`:87); sampling takes it when anchors >= 2 (:159).
  **[Measured]** chr21: 10,954 blocks scored by HMM, 10,502 sampled by HMM, **452
  blocks scored one way and sampled the other**.
  **[Proven fix]** Use one threshold. v1 removes the question by routing every block
  through one sampler, with single-variant blocks as a chain of length 1.
  **[Severity]** MEDIUM.

### The optimizer books privacy for moves that change nothing
- **[2026-09-25] eta_j is 0 for a no-op, so the LP takes it unconditionally**
  **[Symptom]** Selected moves that leave the released alleles identical to the
  target's own. Worked case: block 1077, eps_j = 0.0711, replacement identical to the
  original, so eta_j = 0.00000; `xsol` selects it on both strands and nothing changes.
  **[Cause]** `eta_j` is `sum(utility_loss * (new != orig))`, which is zero when
  nothing differs, while `eps_j` is always positive and is computed from the target's
  ORIGINAL block regardless of what was sampled. Infinite value per unit cost.
  **[Measured]** capacity 0.1, HG00438, chr21: **14,526 of 47,070 selected moves
  (30.9%) changed nothing**, carrying **7.0%** of the total eps. They concentrate on
  LOW-eps blocks (median eps 0.39 vs 2.19 for real moves) because a randomly drawn
  donor is likely to match on common patterns — block 1077 has 82 of 88 donors
  carrying the target's exact pattern.
  **[Proven fix]** Report actual allele change, never `pmi_gain`, as the privacy
  numerator. The utility axis is NOT inflated — eta_j is correctly 0 — so it is the
  privacy-per-utility ratio that is overstated.
  **[Severity]** MEDIUM, and it compounds the sampler degeneracy: the more often the
  sampler returns the target's own alleles, the more free privacy the LP books.

### Utility accounting cannot see the size of what it changes
- **[2026-09-25] 3.13 Mb rewritten, reported as "10.5% of alleles"**
  **[Symptom]** A run reports a modest allele-change percentage while rewriting
  megabases of sequence.
  **[Cause]** `1/support(v)` is per RECORD. A 257 kb allele and a 1 bp SNP cost the
  same. Nothing in the metric suite is sensitive to base pairs.
  **[Measured]** capacity 0.1, chr21: **3,129,810 bp rewritten — 6.7% of the
  chromosome** — across 222 changed alleles over 1 kb, largest 257,485 bp. The same
  run reports "68,335 of 650,541 alleles changed (10.5%)". Non-SNP records are
  changed at a HIGHER rate than SNPs (23.3% vs 14.9%).
  **[Proven fix]** If we adopt this weighting as `phi_t`, we inherit the blind spot.
  Length-weighting is tracked in open question 8.
  **[Severity]** HIGH for interpreting any published utility number.

### The released pipeline matches GRCh38 against GRCh37
- **[2026-09-18] LD blocks and variant mappings built across two reference builds**
  **[Symptom]** Nothing errors. `pangenome_mask.npy` is tiny, `ld_loss` returns
  `(0.0, 0)` (ld_loss.py:123), and every `af_loss` MAF stratum is 0/0 — all
  reported as legitimate values.
  **[Cause]** `get_1000g_phased.sh:3` downloads the 1000G Phase 3 20130502 panel,
  whose header is `assembly=b37`, `hs37d5.fa`, contig `21` length 48,129,895 —
  **GRCh37**. The pangenome is **GRCh38** (`grch38#chr21`, max POS 46,699,788).
  `get_blocks.sbatch:13-16` runs `plink --blocks` on the GRCh37 panel and
  `get_blocks_simple.py:29-34` bins GRCh38 positions into those intervals;
  `get_mappings.py:45-49` matches (POS, REF, ALT) across the two builds.
  **[Measured]** chr21: **1,025** of 340,824 pangenome records match Phase 3
  exactly on (POS,REF,ALT) (0.30%); 2,176 on (POS,REF) (0.64%).
  **[Proven fix]** Use the GRCh38 30x panel via the repo's own **unwired**
  `external/PanMixer/starting_data/scripts/src/get_blocks_grch38.sh`, which downloads
  `1kGP_high_coverage_Illumina.chr{N}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz`
  and renames `grch38#chrN` -> `chrN` before intersecting. The existence of that
  unwired script, plus code branches preferring 30x artifacts, suggests the
  authors' own runs used the GRCh38 path that the released pipeline does not wire in.
  **[Severity]** BLOCKER for any use of their data. **Canary:** print
  `len(pangenome_to_thousand_g_phased)` after `get_mappings`; near-zero confirms it.

### The Li-Stephens sampler never recombines
- **[2026-09-18] "Cohort HMM" is in effect a uniform draw over donors, copied verbatim**
  **[Symptom]** Sampling succeeds and returns plausible haplotype blocks. Nothing
  indicates the mosaic never mixes.
  **[Cause]** Four independent defects in the same path, present identically in
  BOTH copies (`external/PanMixer/tools/panmixer/obfuscate.py` and the pipeline's
  `external/PanMixer/starting_data/scripts/src/hmm.py`, which `get_support_and_pmi.py` imports):
  (a) `EFFECTIVE_N = 1.0/10_000.0` (obfuscate.py:41) is the RECIPROCAL of the
  paper's Ne = 10,000; (b) `get_anchor_snps` returns `(snp, block[1][j])` with
  `snp == block[1][j]` (:82-90), so the "distance" fed to the transition model is a
  difference of VCF ROW INDICES, where the paper specifies centiMorgans;
  (c) the forward pass (:93) and the sampler (:100) use d values a factor 4 apart —
  the paper's Fig. 6 has no factor 4, so the sampler is the deviation;
  (d) the non-anchor fill keys `anchor_pos_to_state` by row indices but iterates bp
  positions (:182-190), so the "is this an anchor" guard never fires and
  `argmin` always selects the last anchor.
  **[Measured]** Executing the shipped class on a 44x200x2 synthetic panel: total
  switch mass 5.8e-12 at a 1-index gap, 5.8e-6 at 1e6; the sampled block equalled
  exactly one donor row of the without-target panel. With the paper's constants the
  switch mass would be 0.13-0.99. That is 10-11 orders of magnitude.
  **[Proven fix]** Do not inherit the sampler. Re-derive the transition constants
  from paper Fig. 6 (Ne = 10,000, r = 1.26, Delta_x in cM, d = Delta_x*Ne*r) and
  supply a real genetic map. Keep ONLY the log-space forward recursion
  (obfuscate.py:140-152), which is correct and O(K) per step.
  **[Severity]** BLOCKER for reusing it as our `baseline_model`, and the reason
  "reuse PanMixer's cohort HMM" is not a plan.

### PanMixer's privacy score is truncated for the rarest blocks
- **[2026-09-18] eps_pmi is capped at ~log(2N), hitting the most identifying blocks hardest**
  **[Symptom]** Plausible per-block privacy values, with no indication that the
  rare tail has been clipped.
  **[Cause]** `forward_algorithm` scores the target against a panel that INCLUDES
  the target (`self.pangenome_haplotypes`, `self.num_states`,
  `external/PanMixer/tools/panmixer/obfuscate.py`:130,132) while normalizing
  transitions by the without-target count. The target always matches itself, so
  p(h_bj) >= 1/2N and therefore eps_j <= ~log(2N).
  **[Measured]** 150 real HMM-path blocks, HG00438, chr21: shipped eps_j spans
  0.071-4.508 (median 1.242); log(88) = 4.4773. Re-scored with the target removed
  from the panel: median 1.299 but max 78.1. For the 9 blocks at the ceiling, the
  leave-one-out value has median 12.06 and max 78.12 — understated up to 17.4x.
  **[NOTE — an earlier version of this entry said eps_j was "pinned near log(2N)
  for every block". That was wrong; see Correction 3 in `docs/memory.md`.]** The
  effect is a ceiling, not a flattening: typical blocks are barely affected
  (median ratio 1.03x), the rarest are truncated hardest.
  **[Related — SUPERSEDED, see the NOTE directly below]** Alleles with AF exactly 0 — the novel, most identifying ones — get
  eps_j = 0, because `-log(0)` is +inf and then zeroed by
  `pmi_subject[pmi_subject > 1e200] = 0` (:232-240). The pipeline copy carries a
  literal `#FIX ME` at get_support_and_pmi.py:65-74.
  **[NOTE — that zeroing branch never executes on our data: 0 of 319,092 chr21
  sites have f_v = 0, because the target's own allele is counted in every frequency
  branch. See Correction 7 in `docs/memory.md`. Hypervariable sites DO score zero,
  but by being excluded as non-anchors, not by this code path.]**
  **[Proven fix]** Never use eps_pmi inside our mechanism (it is a function of the
  private input). If we plot on their axis, recompute it leave-one-out and disclose both.
  **[Severity]** HIGH for interpreting their published frontier.

### The gap-score attack cannot run from the documented pipeline
- **[2026-09-18] apply_mask is never called; the 30x attack DB has no producer**
  **[Symptom]** Following the README end-to-end, `gap_score` fails to find its
  attack database; or, with some files present, it announces a fallback and
  produces numbers against a different, weaker database than the paper's.
  **[Cause]** `apply_mask.sbatch` writes `1000g_phased_masked.npy` and
  `..._posref.npy`, but is NOT referenced by `run_get_starting_data_pipeline.sh`.
  Separately, `diploid_gap_score.py:167-178` prefers a 30x (3,202-sample) attack DB
  — `pangenome_mask_posref_30x.npy`, `1000g_30x_phased_masked_posref.npy` — and no
  script in the repo produces either. `constants.py:38` hardcodes
  `NUM_1000G_SUBJECTS = 3202` while the pipeline downloads 2,504.
  **[Proven fix]** Run `apply_mask` explicitly, and PIN which attack DB is used in
  every reported result. A denser panel is a strictly stronger attacker, so the
  choice moves the privacy axis.
  **[Severity]** HIGH — silent in the sense that it still yields numbers.

### VCF-to-numpy conversion destroys every reason a genotype is missing
- **[2026-09-25] Five different conditions all arrive as -1**
  **[Symptom]** Downstream code cannot tell "this haplotype has no such DNA" from
  "we failed to assemble it", and silently treats them alike.
  **[Cause]** `external/PanMixer/tools/common/VCFtoNP.py` reads only `fields[1]`
  (POS) and `fields[9+i].split(':')[0]` (the GT subfield). The entire INFO column —
  `LV`, `PS`, `CONFLICT`, `AT` — is discarded, as are REF, ALT and the record ID.
  It also maps a HAPLOID genotype to `-1` on strand 1 (3.135% of GT fields), so a
  conversion decision becomes indistinguishable from a data property.
  **[Consequence]** `1/support(v)` counts all five alike, so a variant looks
  "low-support" whether the DNA is absent by construction or merely unassembled; and
  `pmi[haplotype == -1] = 0` zeroes the privacy score for all five.
  **[Proven fix]** For our pipeline, emit an auxiliary REASON array beside the allele
  matrix during conversion. The causes are all recoverable at that moment
  (LV/PS for structure, run-length for gaps, the CONFLICT tag, GT arity) and none of
  them is recoverable afterwards.
  **[Severity]** HIGH for our design — it is the input to `missing_policy`, which
  cannot be made cause-aware if the cause has already been thrown away.

### A variable named `thousand_g_alignments` is the PanGenie callset
- **[2026-09-18] The HMM's anchor set is not what its name says**
  **[Symptom]** Reading the code top-down, `pangenome_to_thousand_g_alignments.pickle`
  reads as "the pangenome-to-1000G map", and the HMM loading it looks like it
  implements the paper's V_1000G anchors.
  **[Cause]** `get_mappings.py:28` binds `thousand_g_alignments` to **`PG.vcf.gz`**,
  the PanGenie callset. The actual 1000G map is a different file,
  `pangenome_to_thousand_g_phased.pickle`, which `hmm.py` never opens.
  **[Consequence]** Anchors come from a GRCh38 source, so the per-block HMM-vs-AF
  routing is NOT degraded by the GRCh37/GRCh38 mismatch — only the LD block
  BOUNDARIES, the attack DB and the af_loss/ld_loss site sets are. Do not conflate
  the two anchor notions. Note the code also does not implement the published
  definition of anchors (V_1000G).
  **[Severity]** MEDIUM, and purely a comprehension trap — it cost this project one
  wrong inference before it was caught (see Correction 5 in `docs/memory.md`).

### The gap-score attack does not exclude the target from the database
- **[2026-09-18] Published equation (7) says max over i != target; the code takes a plain max**
  **[Symptom]** Gap scores systematically depressed; a target can be its own rank-1 hit.
  **[Cause]** `external/PanMixer/tools/downstream/privacy/gap_score.py`:241-242 takes
  `np.max` over the ENTIRE attack database with no self-exclusion. HPRC individuals
  largely appear in 1000 Genomes, so the target's own entry can win.
  **[Proven fix]** Exclude the target row before the max when reproducing
  `eps_private`. `external/PanMixer/tools/common/utils.py`:357-361 does delete the self row, but only
  when the name is present in the panel's sample list — check which panel is loaded.
  **[Severity]** HIGH for any attempt to reproduce the published 0.002 threshold.

### The paper's stated vg versions are mutually incompatible
- **[2026-09-18] A v1.68.0 GAM cannot be read by `vg stats` v1.36.0**
  **[Symptom]** Attempting the published tool combination aborts:
  `terminate called after throwing an instance of 'std::runtime_error' what():
  obsolete, invalid, or corrupt protobuf input`, then `Signal 6 ... VG has crashed`.
  **[Cause]** The published Methods state alignment with Giraffe **v1.68.0** and
  analysis of the resulting GAM files with `vg stats` **v1.36.0**. The GAM
  serialization changed between those releases. Verified by downloading the real
  v1.36.0 binary and running it on a v1.68.0 GAM.
  **[Proven fix]** None available for exact reproduction — the stated combination
  cannot have produced the published numbers as written. Use one version for both
  and STATE which, since the metric definitions may differ across versions.
  **[Severity]** HIGH for reproduction; it removes a route to explaining the
  residual offset in our baseline (see `docs/memory.md`, 2026-09-18).

### vg_prep normalizes against a reference whose contigs cannot match
- **[2026-09-18] `bcftools norm -f` runs before the contig rename**
  **[Symptom]** The step fails to find the contig, because the VCF still uses
  PanSN naming (`grch38#chr21`) while every reference FASTA in the repo names it
  `chr21`.
  **[Cause]** `external/PanMixer/tools/downstream/utility_out/vg_prep.py` orders the
  steps view -c1 -> norm -f -> trim-alt -> **rename-chrs**, so the rename happens
  AFTER the step that needs it. `external/PanMixer/starting_data/references/clean_fa.py` only strips
  contigs whose names contain an underscore; it does not rename to PanSN form.
  **[Proven fix]** Rename first, then normalize. Equivalent, and it runs.
  **[Severity]** MEDIUM — blocks the read-mapping arm until worked around.

### bcftools -r needs an index and silently cannot stream
- **[2026-09-18] `bcftools view -r` in a pipe**
  **[Symptom]** `Could not retrieve index file for '-'` / `Input is not detected as
  bcf or vcf format`, mid-pipeline.
  **[Cause]** `-r/--regions` seeks via the index, so it cannot read stdin.
  **[Proven fix]** Use `-t/--targets`, which streams. Also set `bcftools sort -T`
  to a directory with room; it defaults to /tmp.
  **[Severity]** LOW, and it fails loudly — recorded because it will recur.

### PanMixer writes inside its own checkout
- **[2026-09-18] A pinned dependency mutates itself**
  **[Symptom]** `starting_data/chr*/`, `experiments/`, `slurm/`, `plots/` appear
  inside the pinned clone; a "read-only" dependency is no longer byte-identical.
  **[Cause]** `constants.py:4` sets `BASE_PATH = dirname(abspath(constants.py))`
  with no environment override; `starting_data/scripts/src/paths.py:4` does the same.
  **[Proven fix]** Keep the checkout on `/data` and reach it through the
  `external/PanMixer` symlink; never put it on a full filesystem. Verify the pin
  with `git -C external/PanMixer rev-parse HEAD` and treat tracked-file
  modifications as an error, ignoring the gitignored data dirs.
  **[Severity]** MEDIUM — cost us a relocation before any data landed.

### Shell and exit-code traps in PanMixer
- **[2026-09-18] Failures that report success**
  **[Symptom]** `obfuscate.py` reports success while leaving an UNCOMPRESSED `.vcf`.
  **[Cause]** It resolves `bgzip`/`bcftools` with `shutil.which` and then calls
  `os.system` with the return code IGNORED (obfuscate.py:489-493).
  **[Proven fix]** Assert the `.vcf.gz` exists after the call. Do not trust the
  process exit status.
  **[Also]** `starting_data/scripts/get_genetic_maps.sh:1` is `!/bin/bash` — the
  `#` is missing, so the first line executes. `optimizer.py:293` references an
  undefined `pmi_gain` (the random baseline cannot run) and `stacker.py:140-141`
  raises `NotImplementedError` for its advertised `to_random`. Two registered
  subcommands (`beagle_stats`, `get_alignment_accuracy`) have no dispatch branch and
  silently do nothing; `gather_results` swallows every per-step exception, so a
  missing metric shows up as an absent column rather than an error.
  **[Severity]** MEDIUM, all silent.

## Mechanism / sampler

### ⚠ A GATE STATED AGAINST A BOUND THAT IS NEVER APPROACHED HAS NO POWER — third instance

- **[2026-10-10] `p_succ_bound` cannot detect a miscalibrated mechanism**
  **[Symptom]** A positive control for the optimal attacker required that a mechanism
  miscalibrated with `epsilon_tau` in place of `eta_tau` push the attacker's accuracy
  past `p_succ_bound`. It FAILED at every `tau`. Measured on the toy: at `tau` = 0.9
  the doubled exponent reaches **0.8412** against a bound of **0.9500**; at `tau` = 0.5,
  0.6476 against 0.7500.
  **[Cause]** `p_succ_bound = (1 + tau)/2` is derived from `TV <= tau`, and the worst
  exact TV on this toy reaches only **27.4% to 43.1% of `tau`**. Doubling a quantity
  that sits at 40% of its ceiling leaves it inside the ceiling. The bound is not
  wrong; it is simply never approached, so it has no discriminating power.
  **[Proven fix]** Make the control a PAIRED COMPARISON rather than a threshold test:
  the break must move the attacker's own advantage relative to its correct-mechanism
  value. That fires at every `tau` (+0.0069 at 0.05 rising to +0.1473 at 0.9). The
  test also asserts the threshold version would still fail, so if TV ever does
  approach `tau` the stale reasoning is caught rather than silently inherited.
  **[⚠ THE PATTERN, and this is the third instance in three days]** Clause (b) of the
  haploid toy gate (`tv_empirical <= tau`) passes a doubled exponent at every `tau`.
  The genome-wide `log_z_p` bracket catches one only at `tau` = 0.9. And now
  `p_succ_bound`. **Before writing any gate against a theoretical bound, measure how
  close the real quantity gets to it.** If the answer is "a third of the way", the
  gate cannot detect a doubling. What works in every case here is a saturating
  construction (`ratio_spread`, which the sharp pair attains exactly) or a paired
  comparison against the unbroken mechanism.

### Correlated RNG streams across chromosomes keep every marginal right and break the joint

- **[2026-10-10] Found by the demonstration script's own assertion**
  **[Symptom]** `scripts/toy_genome_demo.py` compared the sampler's empirical joint
  distribution over two chromosomes against the exact enumeration. It passed at
  `tau` <= 0.5 and FAILED at 0.9 and 0.99 — worst deviation 0.03066 against a 4-sigma
  band of 0.02342 at `tau` = 0.99.
  **[Cause]** The loop passed the SAME `seed` to `ffbs.release` for both chromosomes.
  `release` derives its streams from that seed alone, so chrA and chrB drew from the
  same stream and their outputs were correlated. **Every marginal stayed exactly
  correct and only the JOINT was wrong**, which is why it passed at low `tau`: with
  near-uniform distributions there is little structure for a correlation to distort.
  The chromosomes of one release are independent chains and must be seeded as such.
  **[Proven fix]** A distinct seed per chromosome. Confirmed to be a bias rather than
  noise by the scaling: at 10x the draws pure noise would have improved the worst
  deviation only to ~0.0097, still outside the band, and the measured value after the
  fix is **0.00137**. All seven `tau` values now pass.
  **[How to recognise it]** This is the same shape as the 2026-09-30 diploid defect,
  where the two haplotype streams were `seed` and `seed + 1` so one release's
  haplotype 1 was bit-identical to the next release's haplotype 0. **A correlated
  stream is invisible in every marginal and in every per-chain assertion.** Only a
  joint comparison catches it, and only where the distributions carry structure.

### A test that checks the same bound an assertion checks cannot detect the assertion being loosened

- **[2026-10-10] 20 clauses passing, and 3 of 5 mutations of the new guards survived**
  **[Symptom]** `tests/test_genomewide_budget.py` asserted that `log_z_p`, the
  per-chromosome bracket and the genome-wide share total all landed in range, at
  four `tau` values, and it included a clause named for the composition error. All
  20 passed. Then: reverting the diploid bracket to the full `eta_tau`, dropping
  the `budget` argument from `_draw_one`, and disabling `GenomeRelease`'s share
  total — **all three passed every clause.**
  **[Cause]** Two distinct failures, and the first is general.
  1. **Loosening an assertion changes no value.** Every clause checked that the
     computed quantities sit inside their bounds. A mutation that only widens a
     bracket leaves every quantity exactly where it was, so the clause still
     passes while the guard is silent. Verifying a guard needs a test that puts
     the code in a state where the guard SHOULD fire.
  2. **One guard was hidden behind another.** The natural input for "shares must
     sum to 1" — two chromosomes each at share 1 — also doubles the total
     `log Z`, so the genome-wide `log_z` bracket fires on the same input for a
     different reason. "`GenomeRelease` raised" was being read as evidence about
     the share check.
  **[Proven fix]** Tests that assert the guard fires, each ISOLATED so exactly one
  can: the share check tested with shares of 0.3 + 0.3 (every `log_z` comfortably
  inside its bracket, only the total wrong); the `log_z` check tested with correct
  shares and a mis-stamped `eta_tau`; the diploid bracket tested with a `Utility`
  whose declared `share` understates its own weights. All five mutations now fire.
  **[⚠ And one guard could only be covered by injection]** `_draw_one`'s bracket is
  `log Z <= eta_tau * sum(beta)`, which on correct inputs is a THEOREM — no target
  and no `tau` can push a correct forward pass past it, so no black-box test can
  tell the tight bracket from the loose one. Its actual job is to catch a broken
  forward pass, so the test breaks the forward pass: `log_norm` inflated by an
  amount strictly between this chain's share and the whole budget. **A guard whose
  precondition is a theorem has no input that exercises it; it must be injected.**

### A falsification gate that checks the reference against itself catches nothing

- **[2026-10-10] 57 clauses passing, zero power against the two errors it was built for**
  **[Symptom]** `tests/test_diploid_gate.py` in its first form enumerated the exact
  diploid and multi-segment release distributions and checked TV, `ratio_spread`,
  the `log_z` bracket and the budget invariant across a `tau` grid. All 57 clauses
  passed. Injecting the full `eta_tau` per haplotype instead of half into
  `src/mechanism/ffbs.py` — the defect the gate was commissioned for — **still
  passed all 57.** Same for `epsilon_tau` in place of `eta_tau`.
  **[Cause]** Every clause compared `exact.py` to `exact.py`. Both mutations live in
  `ffbs.py`, and nothing in the gate imported `ffbs.py`. A reference implementation
  agreeing with its own gate is two unverified things agreeing — the same defect this
  file records upstream, where PanMixer's shape assertion compares two expressions
  that are equal by construction and can never fire.
  **[Proven fix]** Two clauses that tie the reference to the mechanism, and the
  stronger needs no sampling: **`log_z_p` is deterministic**, so comparing
  `ffbs.release(..., ploidy=2).log_z_p` against `exact.composite_release` pins the
  exponent, the per-haplotype factorisation and the baseline subtraction to 1e-12
  with zero Monte-Carlo error. Verified: with those clauses in place, both mutations
  fire at 5 of 6 `tau`. The sampled-marginal clause is the complement — it checks the
  DRAW follows the factorisation, which is what a shared RNG stream would break.
  **[How to recognise it elsewhere]** Ask of any new gate: *which module does it
  import, and is the thing I am afraid of in that module?* If the answer is only the
  reference, the gate is decoration whatever its clause count.

### A fixture invariant written against donor 0 instead of any pair

- **[2026-10-10] "some donor is another donor minus a cell" checked only donor 0**
  **[Symptom]** `toy.toy_segments` rejected a valid fixture: segment B's third donor
  is the SECOND donor minus a cell, not the first. The check as written would also
  do the reverse and worse — silently PASS a fixture whose only collapsing pair sat
  between two non-zero donors, leaving the output-path collapse untested.
  **[Cause]** The invariant is a property of the donor SET (two state paths emit one
  output path), but it was expressed as a property of donor 0. The collapse is what
  makes TV a statement about OUTPUT paths rather than state paths, so losing it would
  make the headline clause answer a question nobody asked.
  **[Proven fix]** Quantify over all ordered pairs. ⚠ **`toy.toy_cohort` carries the
  same narrow check** and is currently satisfied by donor 0, so widening it is a
  no-op today — but it is the identical latent weakness and should be widened when
  that fixture is next touched.

### ⚠ LATENT AT GENOME SCALE ONLY: a running-float `log_c` accumulator breaches the `log_z_p` tolerance

- **[2026-10-09] Found by arithmetic, before the code was written. NOT yet a bug in
  the repo — recorded so it is not introduced.**
  **[Symptom it would produce]** `log_z_p` fails its `[0, eta_tau]` bracket, or worse
  passes it with a wrong value, **only on a genome-wide run.** Every chr21 test
  passes. There is no error and no symptom at the scale the suite covers.
  **[Cause]** `log_z_p = tilted.log_norm - baseline_log_norm` is a difference of two
  sums whose magnitude grows with T, while the result stays bounded by `eta_tau` —
  catastrophic cancellation that worsens with chromosome length. At chr21 each sum is
  about -584.81; genome-wide it is of order 1e4. Today `log_c` is materialised as a
  `(T,)` array and summed by `np.sum`'s PAIRWISE algorithm, which holds the error to
  ~1.5e-10 at genome scale. The obvious memory optimisation — stop materialising the
  array and accumulate `log_z += math.log(c)` as a running Python float — gives
  sequential summation, whose worst case reaches **~9.9e-05 against the 1e-6
  tolerance**: a 100x breach.
  **[Proven fix]** Keep `log_c` as an array, or use `math.fsum`/Kahan. Both are free.
  Add a regression test that the naive sequential form FAILS at T ~ 1.5e7, so the
  optimisation cannot be reintroduced by someone who only runs chr21.
  **[Why it is in this file]** It is the exact silent-degradation shape this file's
  banner warns about, and it is invisible to the entire current suite.

### ⚠ float64 stops being a preference and becomes a correctness requirement genome-wide

- **[2026-10-09] Arithmetic, not a run. The guard already exists; what was missing is
  the reason it is load-bearing.**
  **[Symptom it would produce]** The mechanism silently becomes an UNTILTED draw from
  `baseline_model`. Every release is a valid path, every gate passes, `tau` appears to
  be honoured, and the target's influence is simply absent. No error.
  **[Cause]** The tilt enters as `psi = exp(eta_tau * beta_v * phi_v)`, per position,
  so the exponent is at most `mean(beta_v) = 1/T`. At chr21 (T = 305,886, tau = 0.5)
  that is 1.80e-06 against float32 eps of 1.19e-07 — only **15.1x** the noise floor.
  Under a genome-wide `beta_v` (standing constraint 4) the exponent falls **BELOW**
  float32 eps, so every potential rounds to exactly 1.0.
  **[Proven fix]** The dtype refusal in `ffbs.forward` already raises on anything but
  float64 and must never be relaxed. Record the arithmetic beside it — done
  2026-10-09 — so the next person reading "float64 is load-bearing" knows why.
  float64 itself is comfortable: the boundary is `exp(x) == 1.0` iff `x <= eps/2` =
  1.110e-16, which a genome-wide exponent clears by ~1.6e+08.

### A dense-transition memory figure was wrong by 1000x, in two files

- **[2026-10-09] `21 TB` where the arithmetic gives `21.1 GB`**
  **[Symptom]** `src/mechanism/ffbs.py` and `src/mechanism/hmm.py` both justified the
  rank-one decomposition by claiming a dense `A_t` per position would be "about
  21 TB on chr21". `hmm.py` printed its own arithmetic one line later —
  `340,849 * 88 * 88 * 8 B` — which is **21,116,277,248 B = 21.1 GB**.
  **[Cause]** A unit slip, copied between the two modules, that survived because the
  conclusion it supports (never materialise `A`) is correct either way, so nothing
  downstream disagreed with it.
  **[Proven fix]** Corrected in both files with the arithmetic kept inline. ⚠ The
  lesson is the one this project keeps relearning: **a number that supports a correct
  conclusion is never checked.** When a docstring shows its own arithmetic, recompute
  it — this one was off by three orders of magnitude for six weeks.

### A diploid assertion message printed the wrong field of the wrong chain

- **[2026-10-09] The bracket fires correctly and then misreports why**
  **[Symptom]** The diploid `log_z_p` bracket in `ffbs.py` raised with
  `f"(haplotype values {draws[0][2]!r} and {draws[1][3]!r})"`. Index 2 is `log_z`;
  index 3 is `u_achieved`. So the message printed chain 0's `log_z` beside chain 1's
  **utility**, as though both were log-normalisers.
  **[Cause]** An index typo on the message path only. The assertion's own arithmetic
  uses `lz`, which is correct, so the guard worked and only its explanation lied.
  **[Proven fix]** `draws[1][3]` -> `draws[1][2]`. ⚠ **No test covers this** — the
  message is only produced on the failure path, and the suite asserts that the bracket
  raises, not what it says. That is the general hole: an error message is code that
  only ever runs when something is already wrong, so it is the least-exercised and
  most-trusted text in the repo.

The classes most worth watching were known from the proposal before any code
existed and are written as constraints in `CLAUDE.md`: a target-dependent
`output_support` or `baseline_model`, an unnormalized `u_path`, `beta_t` not
summing to 1, a `diploid_utility` that sums instead of averaging, and argmax
substituted for `ffbs`. Every one of them produces a plausible path and no error.
Each now has a test. The entries below are what the first build actually hit,
which was mostly a DIFFERENT set.

### Six defects an 88-test suite passed over, found by an independent audit
- **[2026-09-30] The implementer wrote the code, the docstrings AND the tests**
  **[Symptom]** 88/88 tests pass, the toy gate is sharp and saturated, and the
  headline finding is confidently wrong. Every check the implementer thought to
  write, passed.
  **[Cause]** Self-audit. Six independent auditors were given the docs as the
  authority and told explicitly NOT to trust the docstrings, because the same
  agent wrote the claims and the tests that check them. 30 gaps survived
  three-lens adversarial verification. The ones that were real bugs:
  1. **`epsilon_tau` returned `2*arctanh(tau)/delta_u`.** `delta_u` CANCELS in the
     budget. At `delta_u` = 2 it returned exactly `arctanh(tau)` — the expression
     `docs/terms.md` forbids in capitals. Unreachable and untested at any
     non-unit range, so the one scenario the parameter exists for was the one that
     corrupted the headline privacy number.
  2. **The diploid `log_z_p` was never bracketed.** Each haplotype was checked
     against half of `eta_tau`; their SUM was returned unchecked. Verified by mutation:
     restoring the full `eta_tau` per haplotype was ACCEPTED at `tau` = 0.25, 0.5
     and 0.9, returning a record whose published `log_z_p` (0.763 at `tau` = 0.5)
     violated its own stated invariant of [0, 0.549].
  3. **The diploid seeds were `seed` and `seed + 1`.** So `release(seed=s)`
     haplotype 1 was bit-identical to `release(seed=s+1)` haplotype 0 — 0 of
     305,886 positions differed. "3 independent draws" were 4 correlated haplotype
     draws; `sd(u)` was inflated 2.2x and landed coincidentally on the HAPLOID
     value when the diploid prediction is that over sqrt(2).
  4. **`release()` took the allele matrix as a free argument.** `CohortHMM` carries
     only the availability mask, so a caller could hand the sampler a
     target-pruned matrix of the same shape — the standing-constraint-1 violation
     the cohort/target directory split exists to prevent and the one shape it
     cannot catch. The post-draw `-1` check caught only the occupied subset.
  5. **`phi_v` was never called by the mechanism.** `tilt_potentials` and
     `u_of_output` both re-inline the indicator for speed, so the function
     carrying the documented semantics — and the only one with a test — was on no
     live path. Two implementations of one definition with no edge between them.
     The same was true of `exact.utility_vector` against `utility.u_of_output`,
     which differ in their null guards and which the gate and the mechanism use
     respectively.
  6. **`sum_v beta_v^2` was reported as `mean(beta_v)`.** Two statistics that
     coincide for uniform weights. The real value is 4.556723e-06, not 3.269192e-06
     = 1/T. It is the number the "lever" claim rested on.
  **[Proven fix]** All six fixed; 1-4 verified by mutation before the fix was
  accepted. Seven regression tests added, including one that PARSES the
  calibration table out of `docs/terms.md` — the previous three "checks" of it
  were three independent transcriptions with no edge to the document.
  **[Severity]** ⚠ The class is what matters, not the six. A suite written by the
  same agent as the code tests the agent's model of the code, not the spec. The
  cheapest counter is an auditor holding the spec who is told the docstrings are
  the claim under audit.

### Two no-op-rate measurements of the same run disagree by 2.2x
- **[2026-10-08] 30.9% recorded on 2026-09-24; 66.5% measured now, same capacity and seed**
  **[Symptom]** Both figures describe "selected moves that changed nothing" for
  PanMixer at `capacity` 0.1, HG00438, seed 123. The 2026-09-24 entry says 14,526 of
  47,070 (30.9%); recomputing from `xsol.npy` and `new_haplotypes.npy` gives 31,315
  (66.5%).
  **[Cause]** UNRESOLVED. Candidates: a different notion of "changed" (the earlier
  pass may have looked only at a block's ANCHOR variants rather than all its sites);
  a different block mapping than `simple_blocks_idx.npy`; or different treatment of
  blocks where the target's genotype is missing throughout, which the current
  definition counts as no-ops.
  **[Proven fix]** None yet. Both are recorded with their definitions and neither is
  quotable as settled. The cells-changed counts do NOT depend on the no-op definition
  and are unaffected.
  **[Severity]** ⚠ A 2.2x disagreement on a number used to characterise the other
  arm. It is in the comparison arm, not our mechanism, so nothing of ours rests on
  it -- but it must be closed before either figure appears in a writeup.

### Model B silently created two coordinate systems and could not release
- **[2026-10-08] Dropped block units left `unit_spans` and `n_positions` disagreeing**
  **[Symptom]** `release()` refused with a shape error under model B on chr21.
  Marginals and `expected_mismatches` worked, so the comparison table looked fine.
  **[Cause]** `from_arrays` dropped the 4 units no donor can cover but left
  `unit_spans` indexing the full 305,886-position chain, while `n_positions`
  summed only the surviving spans (304,664). `expand_to_positions` then wrote past
  the end of its output array. ⚠ **And the fidelity shares were wrong**: the
  numerator covered 304,664 positions while the denominator was the full chain's
  506,788-plus-1,222 known cells.
  **[Proven fix]** Refuse instead of dropping. `from_arrays` raises and names the
  remedy; `hmm.dead_block_positions()` returns the mask to compose into `in_chain`,
  so the exclusion lives on the axis every artifact shares -- the same treatment
  `support == 0` already gets from `genmap.chain_mask`. Re-measured on the common
  axis, the floors are C 14.8073% and B 14.5959%, so B's advantage was real and not
  a denominator artifact.
  **[Severity]** ⚠ The guard caught it rather than producing a plausible wrong
  path, which is the guard working. But the quantity that DID run -- the marginals
  -- was reported against a denominator from the other coordinate system, and
  nothing flagged that. **A derived quantity that works while the main path is
  broken is a warning, not a reassurance.**

### A docstring that described the model and code that implemented a different one
- **[2026-10-07] Model B's transition used start-to-start, the docstring said the gap**
  **[Symptom]** Model B ran, every test passed, and `expected_switches` came out
  **5,270.9** against the **3,532** an independent calculation had predicted on
  2026-09-28. A 49% discrepancy in a quantity with a published prediction.
  **[Cause]** `from_arrays` set each unit's genetic position to its FIRST position
  and let `switch_mass_from_cm` diff those, which is start-to-start spacing. That
  silently charges a unit's own internal genetic length to the transition OUT of
  it. The docstring two paragraphs above said "the transition uses the GAP between
  units, from the last position of one to the first of the next, so the within-unit
  genetic length is discarded" -- which is also what `docs/plan.md` specifies and
  what makes model B "discard the within-block recombination". The prose was right
  and the code was wrong.
  **[Proven fix]** Compute the inter-unit gaps explicitly and build the transition
  array from their cumulative sum, keeping each unit's real first-position cM in
  `genetic_pos` for provenance. `expected_switches` then lands at **3,394.1**
  against the predicted 3,532. A test now asserts the gap form AND asserts that the
  start-to-start form would differ on the same fixture, so the two cannot be
  confused again.
  **[Severity]** ⚠ The whole model-B result set was wrong before this: the floor
  moved from 14.78% to 14.57% of cells and the tilt's slope from -1.3355 to -1.7738
  cells per unit `eta_tau`. **Caught only because a prediction existed to compare
  against** -- no test failed, and the numbers looked entirely plausible. Where a
  quantity has an independently predicted value, check it.

### A measurement script that did not measure the mechanism
- **[2026-10-03] fid2.py tilted EACH haplotype at the full `eta_tau`, not half**
  **[Symptom]** A plausible differ-count table at every `tau`, with sensible
  magnitudes. Nothing in the output said the mechanism being measured was not the
  mechanism in `src/`.
  **[Cause]** The script built potentials and called `F.forward` DIRECTLY instead of
  going through `ffbs.release()`. `release()`'s diploid path runs each chain at
  half of `eta_tau` so the two `log_z_p` values sum into the bracket; the script used
  the full `eta_tau` on both. **That is the "SUM instead of MEAN" aggregator error
  the project has a standing constraint about** — and calling `forward` directly
  also bypassed the diploid `log_z_p` bracket added on 2026-09-30, which exists to
  catch exactly this and would have fired.
  **[Measured]** Every effect in that run is exactly 2x the mechanism's: slope
  -2.4761 cells per unit `eta_tau` against the true -1.2381, confirmed by running
  both conventions side by side in `scripts/exact_fidelity.py`.
  **[Proven fix]** Measure through the public entry point. A diagnostic that
  reconstructs the mechanism's internals by hand is a second implementation with no
  test, and here it silently became a different mechanism. Where a diagnostic must
  touch internals, assert the invariant the entry point would have checked.
  **[Severity]** ⚠ Silent, and it defeated a guard built specifically for this error
  two days earlier. Found by independent audit, not by the suite.

### 50 draws spent on a quantity that is exactly computable
- **[2026-10-03] E[differ] is a linear functional of the posterior marginals**
  **[Symptom]** Two sampling runs (n = 10 then n = 40), a sign flip between them,
  and an argument about what sample size would settle it.
  **[Cause]** `differ` is a sum over positions of a function of the state at that
  position, so `E[differ] = sum_t sum_k gamma_t(k) 1[hap[k,t] != target_t]` with
  `gamma ∝ alpha * beta` — exact, zero Monte-Carlo error, O(TK). **The project had
  already written this lesson down** for expected utility, in the
  common-random-numbers entry below: "report the exact quantity instead ... three
  forward passes and no sampling". It was not applied to the metric that replaced
  it.
  **[Proven fix]** `ffbs.state_marginals` and `ffbs.expected_mismatches`, validated
  against the repo's own O(K^T) enumerator on the toy to **4.8e-15**. chr21 runs in
  ~30 s per `tau`. The exact answer: the tilt reduces expected mismatches by
  **1.2381 cells per unit `eta_tau`**, perfectly linearly — 0.680 cells at
  `tau` = 0.5 and 1.823 at `tau` = 0.9, of 507,692.
  **[Severity]** Wasted effort, and worse, it put a sample-size argument where a
  closed form belonged. ⚠ The general rule now: before estimating anything about
  this mechanism by sampling, ask whether it decomposes over positions. If it does,
  it is exact.

### A pattern in n=10 noise that reversed at n=40
- **[2026-10-03] Right-signed and monotonic at n=10; opposite sign at n=40**
  **[Symptom]** Paired `differ(tau) - differ(0)` came out -0.1, -7.0, -21.2 at
  `tau` = 0.1/0.5/0.9: correctly signed (the tilt should IMPROVE fidelity), monotonic in
  `tau`, and with 10 of 10 seeds moving. A tidy, publishable story. At n = 40 the same
  measurement gives **+11.8 and +9.8** — the wrong sign. No n reaches significance.
  **[Cause]** The true effect is **1.80 cells of 507,692**, derivable exactly from
  `d log Z / d eta` without sampling, against a paired sd of 51-61 cells. Neither n could
  resolve it. The monotone-looking sequence was three draws from noise, and the ordering
  was read as structure because theory PREDICTED that ordering — so the prediction supplied
  the pattern and the data merely failed to contradict it.
  **[Proven fix]** Derive the expected effect size BEFORE interpreting the sample, from the
  exact route where one exists. Here one did, and it says no feasible n matters. Also:
  increase n before believing a sign, especially when the sign is the one you expected.
  **[Severity]** ⚠ Would have been a published finding. It is the same failure as the
  common-random-numbers entry below but one level up: not a bad measurement, a bad
  inference from a defensible measurement. The tell was available — three points, a
  predicted direction, and |t| < 2 throughout.

### A headline finding measured under common random numbers
- **[2026-09-30] "0 of 611,772 positions differ" was true of one seed and predicted in advance**
  **[Symptom]** The released path was reported bit-identical at `tau` = 0, 0.1, 0.5
  and 0.9. It looked like a strong, clean result.
  **[Cause]** Two compounding errors. `_draw_one` opens a fresh generator from the
  same seed at every `tau`, so the paths were compared under COMMON RANDOM NUMBERS
  — with a per-position tilt of 1.8e-06, agreement is the expected outcome, not a
  finding. And n = 3 was too small: at ten disjoint-stream draws, **8 of 10 move**,
  by a mean of 417 state positions per haplotype. The one seed checked was one of
  the two that happen not to. ⚠ Worse, the docs PREDICTED the result they then
  reported as a discovery — "the expected number of positions the tilt moves is of
  order `eta_tau`" — so the measurement could not have surprised anyone.
  **[Proven fix]** Report the exact quantity instead: `E_Q[u] = d log Z_p/d eta`
  and `Var_Q(u) = d2 log Z_p/d eta2`, three forward passes and no sampling. The
  finding's substance survives unchanged — expected utility moves 2.94e-06 across
  the whole grid — and is now stated as "the tilt reshuffles the path and the
  reshuffling is utility-neutral", which is both true and a stronger claim.
  **[Severity]** ⚠ Silent, and in the most dangerous direction: a wrong headline
  with a correct conclusion underneath it. Also a lesson about the project's own
  reporting rule — a row whose columns are all identical deserves suspicion before
  it deserves a synthesis paragraph. The rounded MEAN concealed movement that the
  per-draw records showed: `donor_switches` 14,147/14,211/14,147 at `tau` = 0.9
  against 14,148/14,212/14,148 below it, both printing as 14,169.

### `log_z_p` lands ~1000x outside `[0, eta_tau]` and the tempting fix is fatal
- **[2026-09-30] The forward normaliser is not the tilt normaliser when a mask is applied**
  **[Symptom]** `log_z_p` from `sum_t log c_t` is **-584.81** on chr21 against a
  bracket of `[0, 0.5493]` at `tau` = 0.5. On the toy it is -0.575 against
  `[0, 0.5493]`. The sampler is otherwise correct and every other check passes, so
  the obvious reading is that the assertion is too strict.
  **[Cause]** `missing_policy` = RENORMALISE enters as a 0/1 potential, so `R_D` is
  the Li-Stephens chain CONDITIONED on never occupying a donor with no called
  allele — and that conditioning has its own normaliser, which `sum_t log c_t`
  contains along with the tilt.
  **[Proven fix]** `log_z_p` = the tilted pass's `sum_t log c_t` MINUS
  `baseline_log_norm(hmm)`, an untilted pass whose only potentials are the mask.
  It is target-independent by construction, so compute it once per cohort and
  cache it. Verified against brute-force `log sum_y R_D(y) exp(eta u)` on the toy:
  agree to < 1e-12 across the whole `tau` grid and 417 inputs.
  **[Severity]** ⚠ Silent, and the WORKAROUND is the danger, not the bug. Widening
  the bracket to accommodate the wrong quantity would have retired the cheapest
  correctness assertion in the project while leaving it in place as decoration.
  `tests/test_mechanism.py::test_log_z_p_needs_the_baseline_subtraction_and_says_so_when_it_is_missing`
  asserts the failure fires.

### RENORMALISE as a renormalised transition row divides by zero
- **[2026-09-30] A masked transition row sums to exactly 0 across a zero-cM step**
  **[Symptom]** Reading `missing_policy` as "mask row `j` to the allowed set and
  renormalise" works everywhere except at a zero-distance step, where it produces
  `0/0`. On a chromosome where zero-cM steps are rare it would surface as a handful
  of NaNs late in a long run.
  **[Cause]** Across a step with `Delta_x = 0` the transition IS the identity, so a
  donor available at `t-1` and unavailable at `t` has nowhere to go: the masked row
  is all zeros and there is no chain to renormalise. Measured: **1,058
  zero-distance steps inside chr21 `chain_span`** against 4.26% mean missingness,
  so it is reachable on the real data rather than a contrived case.
  **[Proven fix]** Apply the mask as a 0/1 POTENTIAL, not as a modified transition
  — which is also the form `docs/plan.md` had already pre-registered. A genuinely
  dead-ended chain then raises, naming the position and the cause, instead of
  returning NaN. See the `log_z_p` entry above for the consequence.
  **[Severity]** Would have been silent on a chromosome with no zero-cM steps and
  wrong on every chromosome that has them.

### A test that reads prose rather than behaviour keeps lying
- **[2026-09-30] The "no argmax" gate fired on the docstrings explaining why there is no argmax**
  **[Symptom]** `test_no_map_decoding_anywhere_in_the_mechanism` reported five
  offenders in `ffbs.py`, all of them lines of prose warning against MAP decoding.
  **[Cause]** The check searched the source TEXT. This is the SECOND time this
  exact mistake has been made here — `tests/test_preprocess.py` records it for the
  shared-parser check and states the lesson in the docstring — which suggests the
  text-grep formulation is the natural thing to reach for and needs a standing
  answer, not another one-off fix.
  **[Proven fix]** Walk the AST and match `ast.Name.id` / `ast.Attribute.attr`.
  Docstrings are string constants and comments are absent from the tree, so prose
  cannot trip it while real uses cannot hide from it.
  **[Severity]** Noisy rather than silent — but a gate that cries wolf gets
  weakened, and weakening this one would have removed the only structural defence
  against a decoder that carries no guarantee.

### A Monte-Carlo gate reports 19.9 sigma on a provably correct sampler
- **[2026-09-30] Largest deviation divided by the SMALLEST standard error**
  **[Symptom]** `test_ffbs_matches_the_exact_release_distribution` failed at
  `tau` = 0 — the one setting where the sampler is just a draw from `R_D` — while
  the state-path JOINT test on the same sampler passed comfortably.
  **[Cause]** `np.abs(got - want).max() / se.min()`. The max deviation is at a
  high-probability path and the min standard error at a low-probability one, so the
  ratio compares two different paths.
  **[Proven fix]** Divide elementwise, then take the max:
  `(np.abs(got - want) / se).max()`. Reports 1.96-2.67 over 40,000 draws.
  **[Severity]** A false FAILURE, which is the benign direction — but it pointed at
  the sampler, and the natural response to a failing statistical gate is to loosen
  its threshold, which would have blinded it for real.

### A utility figure quoted against a ceiling of 1 that is actually 0.815
- **[2026-09-30] `utility_retained` has no meaning without `u_ceiling`**
  **[Symptom]** `utility_retained` = 0.6940 on HG00096 looks like 69.4% of
  achievable. It is 85.1%.
  **[Cause]** `docs/terms.md` defined the ceiling as "the target's own path,
  u = 1 by construction". True for a target with a call at every weighted position;
  false for a real external target. Under `beta` scope `chain`, 16.4% of the chr21
  weight sits at positions where the target has no call, so the ceiling is
  **0.815455** over 253,846 of 305,886 chain positions — deferred decision R1
  measured rather than predicted.
  **[Proven fix]** `u_ceiling` is now a PRIMARY term, computed per target and per
  `beta` scope, stamped on every release record, and emitted as its own row by the
  scoring script. Report both or neither.
  **[Severity]** Silent and bidirectional — against 1 it understates the mechanism,
  against 0.8155 unlabelled it overstates it.
