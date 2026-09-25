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
  **[Severity]** HIGH for us — it is a constraint-1 violation introduced by our own fix.

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
  **[Related]** Alleles with AF exactly 0 — the novel, most identifying ones — get
  eps_j = 0, because `-log(0)` is +inf and then zeroed by
  `pmi_subject[pmi_subject > 1e200] = 0` (:232-240). The pipeline copy carries a
  literal `#FIX ME` at get_support_and_pmi.py:65-74.
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

No entries yet — no code exists. When the first one lands, the classes most worth
watching are already known from the proposal and are written as constraints in
`CLAUDE.md`: a target-dependent `output_support` or `baseline_model`, an
unnormalized `u_path`, `beta_t` not summing to 1, a `diploid_utility` that sums
instead of averaging, and argmax substituted for `ffbs`. Every one of them
produces a plausible path and no error.
