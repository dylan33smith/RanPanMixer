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

### A PanMixer evaluator dies on a missing `xsol.npy`, and the fix is a comparability decision

**[Symptom]** `scripts/beagle_attack_one.py --experiment 6` raises
`FileNotFoundError: .../exp_6/data/chr21/0/xsol.npy` when run on PanMixer's own releases.

**[Why it is not a missing file]** `accuracy_stats.py:162` loads `xsol.npy` and builds
`xsol_expanded`, then ANDs every variant subset with it, so **every metric is computed only
over cells inside selected blocks**. Our releases were staged with an ALL-ONES mask — the
declared substitution in `scripts/beagle_attack_setup.py`, because a path release has no
selection step. PanMixer's releases carry a REAL mask from their optimiser (14,342 blocks
at `capacity` 0.01, 177,865 at saturation). Dropping their mask in would score their
releases over a SELECTED SUBSET while ours were scored over EVERY CELL — different
denominators on two sides of one table.

**[Proven fix]** Stage BOTH and report both: `exp_6` with an all-ones mask (matched to
ours, the comparable column) and `exp_7` with their own `xsol` (their native protocol,
reported separately and labelled not-comparable). `scripts/arm1_pm_beagle.py --stage` does
this. ⚠ The crossover where Beagle beats the baseline MOVES between the two masks —
`capacity` 0.01 and 0.05 under all-ones, 0.01 only under theirs — so the choice is not
cosmetic.

**[Found by]** Validating ONE row before launching eight, 2026-10-15.

### A constant-answer attacker reports balanced accuracy 0.5 and reads as a clean null

**[Symptom]** A3's `D1_extremal` row is exactly 0.5000 at all seven `tau` — the tidiest
null in the table.

**[Why it is wrong]** Measured over all 1,750 trials: **TPR = 1.0000, FPR = 1.0000**. The
attacker answers "p" every single time regardless of which hypothesis produced the release.
Balanced accuracy of an attacker that always gives the same answer is 0.5 BY CONSTRUCTION,
so the number carries no information about the mechanism. The cause is structural:
PanMixer's score credits agreement with non-negative weights and never penalises
disagreement, and `D1_extremal` differs from the target at 100% of called positions.

**[Proven fix]** Check TPR and TNR separately before reading any balanced accuracy as a
null. A genuine null has BOTH near 0.5 (A1 `D1_extremal`: 0.4947 / 0.4767; A3 `D3_random`:
0.4566 / 0.5326). Report A3's informative conditions (`D2_nearest`, `D3_random`) and label
`D1_extremal` degenerate wherever it appears.

**[Found by]** Dylan asking for the `D1_extremal` row specifically, 2026-10-15. ⚠ Same
failure family as A3's positive control, which was also algebraically unable to fail.

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

### A finished 3.4-hour fan-out is destroyed by its own reporting step
- **[2026-10-25] one knob, two names and two spellings, in `exact_kl_stable.py`**
  **[Symptom]** All 22 autosomes print ✅ over 3.4 hours, then
  `KeyError: '0.0'` at `kl = math.fsum(per[c]["kl"][t][str(tau)]["kl"] ...)`, and
  NOTHING is saved — no JSON, no partial results.
  **[Cause]** Two independent mismatches at one call. (1) The worker keyed its results
  `f"{v:g}"` (`eta` 0.0 -> `'0'`) and the aggregator asked for `str(tau)` (`'0.0'`).
  (2) The aggregator read `a.taus` while the run was driven by `--etas`, so it walked a
  grid the workers had never evaluated. The results were lost rather than merely
  unreported because they were held in MEMORY and `a.json.write_text(...)` sat
  DOWNSTREAM of the crash — a formatter could destroy a finished computation.
  **[What the wrong output looked like]** 22 green ticks and a traceback. Every
  expensive thing had SUCCEEDED; only the cheap assembly failed, which is the most
  expensive place for a bug to be.
  **[⚠ The trap when fixing it]** Repairing the reported `KeyError` alone buys another
  wasted 3.4 hours. An audit found the run would then die at `p_succ_bound(250.0)`
  (`ValueError`, it requires `tau` in [0,1)), and after that at `_report`'s
  `max()` over an empty generator — both still upstream of the JSON write. Three
  crash sites, one cause, and the first two are invisible until hour three.
  **[Proven fix]** `knob_key()` — ONE spelling (`repr(float)`), round-trip asserted,
  called by producer and consumer so they cannot drift; `:g` was also lossy at 6 s.f.
  and would have silently merged two near grid points into one cell. One `grid` object
  used at every former `a.taus` site. `p_succ_bound` fed `tanh(eta)` with the saturated
  case reported vacuous. **JSON write moved ABOVE the report**, carrying `per_chrom`.
  Per-chromosome checkpoint JSONL with `--resume`, reusing a chromosome only when arm,
  weighting, grid and target set all match.
  **[Gate]** `tests/test_exact_kl_eta_mode.py`, 13 tests, 4 mutation-tested (each
  injected defect fired the matching gate and was restored).
  **[Generalisation]** Occurrence FIVE of one shape. Four CRN defects addressed a
  knob's value by its POSITION; this addressed it by two SPELLINGS. A value written in
  one place and read in another needs exactly one function that names it, and that
  function needs an assertion. Separately: any fan-out over ~20 minutes must checkpoint
  — the cost of a late bug is not the bug, it is the hours.


### A keyword default silently selected a different utility function for a whole attack suite
- **[2026-10-25] `utility.beta_v(sup, in_chain=ic)` with no `weights=`**
  **[Symptom]** NONE. Nothing raised, nothing warned, the controls passed and the
  tables printed plausible numbers. `logs/B_eta_sweep/a1_eta.md` and `a34_eta.md`
  looked like finished results for the released tool and were not.
  **[Cause]** `utility.beta_v` has `weights: str = "inv_support"`. `arm2_attack.py:99`
  called it positionally and never passed the keyword, so the whole Arm 2 attack suite
  ran at `inv_support` while every metric it was to share a table with had been re-run
  at `uniform` on 2026-10-24. `arm2_attack_a3a4.py` imports `load_state` from that
  module, so one omission took out BOTH attacks.
  **[What the wrong output looked like]** Identical in FORM to a correct run. The only
  visible tell was arithmetic: `ETA_GRID`'s top point was 855, and 1/sqrt(1.3741e-06)
  = 853 is `inv_support`'s chr21 `eta*` while `uniform`'s is 756. A grid built from the
  wrong weighting's `eta*` is the fingerprint; the accuracies themselves gave nothing
  away.
  **[Why it matters]** chr21 `var_q_u` is 1.3741e-06 against 1.7500e-06, so one `eta`
  sits at two different `acc_ceiling`. Constraint 5 forbids two utility functions in a
  comparison and constraint 7 requires matched `acc_ceiling`, so the attack rows and
  the privacy rows could not have gone in one table.
  **[Proven fix]** `arm2_attack.WEIGHTS`, set from `--weights` on both entry points
  (default `uniform`), threaded into the single `beta_v` call; `load_state` now RAISES
  if its cached state was built under a different weighting rather than serving it.
  A3/A4 sets the global through the module object (`A1.WEIGHTS = a.weights`) because
  `load_state` reads A1's global, not its own.
  **[Gate]** `tests/test_beta_weighting_explicit.py` — an AST walk asserting every
  PRODUCTION `beta_v` call passes `weights=`, plus a frozen `LEGACY` census so a NEW
  omitting file fails. Mutation-tested twice (restored the omission; flipped the CLI
  default) and both named the right test. ⚠ The walk found **54 omitting sites in 23
  files**; 5 production fixed, 23 legacy frozen.
  **[Generalisation]** A keyword default is a silent wrong answer whenever two values
  are both valid and name DIFFERENT quantities. Grep for callers, not for errors.


### An index array read as a boolean mask, so a results table never ran at all
- **[2026-10-22] `pangenome_mask.npy` is site INDICES, not a boolean**
  **[Symptom]** `scripts/arm1_attacks_results.py` raised `ValueError: operands could
  not be broadcast together with shapes (340824,2) (258610,2)` in `load_changed`.
  **[Cause]** The file holds **258,610 site INDICES** into the 340,824-site chr21
  axis, and the script did `.astype(bool)` on it. That is wrong twice: it is True at
  every nonzero index (258,609 of 258,610 entries) and it is the wrong LENGTH, which
  is what raised. PanMixer's own `tools/beagle/accuracy_stats.py:43` builds the
  boolean correctly (`snp_mask = zeros(n, bool); snp_mask[snp_indices] = True`) and
  `gap_score.py:101` uses it as indices.
  **[Proven fix]** Build the boolean the way their code does. ⚠ Nothing published was
  affected — it raised before returning, and `snp_only` is the default subset — but
  the consequence is worth recording: **the 308-row pinned-grid attack matrix had
  never been tabled**, so 396 evaluator runs plus 308 Beagle runs existed only as
  unread result files on disk. A script that cannot run is not a script that produced
  a conservative answer. See `docs/memory.md` 2026-10-23.

### A directory namer refuses the parameter it is asked to name
- **[2026-10-22] `releases.tau_dir` was pinned to two decimals, by design**
  **[Symptom]** Every release at `tau` = 0.999, 0.9999 or 0.99999 was refused with
  "does not survive the two-decimal directory name".
  **[Cause]** It formatted `%.2f` and raised unless the name read back exactly —
  deliberate, and right: without it, 0.995 and 0.99 would both render `tau_0p99` and
  the second release would land in the first one's directory, a silent collision
  between two PRIVACY LEVELS. The guard was correct; its two-decimal assumption was
  not.
  **[Proven fix]** Name losslessly — two decimals where that is exact (so every
  pinned directory keeps its name), otherwise the shortest exact decimal — and keep
  raising on any name that does not round-trip. Four scripts that built the same path
  by hand with `f"tau_{tau:.2f}".replace(".", "p")` now call the one function, so the
  rule cannot be bypassed. `tests/test_releases.py` asserts the names are distinct and
  each reads back exactly.

### A `tau` grid re-run redraws with ANOTHER `tau`'s seed, and the releases come out bit-identical
- **[2026-10-22] Seeding a cell from the `tau`'s POSITION in `--taus`, not its VALUE**
  **[Symptom]** `scripts/arm1_release.py --taus 0.999 0.9999 0.99999` wrote 132
  plausible releases, and then PanMixer's `gap_score` on them returned
  `genotype_g_to_gstar` identical to the PINNED grid's first three columns to 13
  significant figures (-2,935.9 / -2,998.2 / -3,554.3). `obfuscated_accuracy` matched
  to six decimals. Nothing errored; every file, sidecar and completeness assertion
  passed.
  **[Cause]** `genome_seed(target, tau_index)` took `tau_index` from
  `enumerate(taus)`, i.e. the caller's list position. A subset run therefore passed
  0, 1, 2 and reproduced the seeds of `tau` = 0, 0.05 and 0.10 — 4238719590,
  2209799234, 1378656232. Because the tilt is far too weak to move the draw on one
  chromosome, the chr21 slices were **BIT-IDENTICAL** to those releases (0 cells
  differing of 681,648), so every DRAWN-path quantity was the pinned grid's number
  wearing a new `tau` label. ⚠ This is occurrence FOUR: the entry below records it
  for a headline finding (2026-09-30), `docs/memory.md` records seed reuse across
  chromosomes, and `scripts/arm2_attack.py` carries fourteen lines of comment about
  it. **A comment in another file is not a gate.**
  **[What was NOT affected, and it is measured]** The exact endpoint.
  `ffbs.state_marginals` uses no RNG and `scripts/arm1_exact_endpoint.py` has no seed
  at all, so redrawing all 132 releases under different seeds reproduced **792/792
  exact values bit-identical** (`e_differ`, `e_q_u`, `target_fidelity`, `u_ceiling`,
  `e_q_u_over_ceiling`, `called_cells`; worst relative difference 0.000e+00).
  **[Proven fix]** Key the seed on the value's index in a FIXED
  `SEED_GRID = PINNED_TAUS + EXPLORATORY_TAUS`, and REFUSE a `tau` outside it rather
  than invent an index. Append, never prepend: the seven pinned points keep the
  indices their 308 releases were drawn under, verified against the seeds recorded on
  disk (42 of 42 pinned directories agree; the 18 exploratory ones disagree because
  they were the bad draws). Gated by `tests/test_arm1_release_seeding.py`, five tests,
  mutation-tested — the positional seed fails 3, prepending fails the reproducibility
  test, dropping the target fails 3. See `docs/memory.md` 2026-10-23.
  **[How to recognise it again]** Two `tau` columns of a DRAWN quantity agreeing to
  more digits than a median over 44 targets could. Check the seed in the release
  filename: it embeds `_seed<n>`, so a collision is visible with `ls`.

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

### The forward pass dead-ends at a zero-distance step once the posterior concentrates
- **[2026-10-23] An `eta`-DEPENDENT zero-mass failure, distinct from decision K**
  **[Symptom]** `scripts/eta_sweep.py --genome` at `eta` = 1.997e8 raised, on chr7
  and chr10 only, from `ffbs.forward`: *"chain position 1154456 (pos 142603528) has
  zero total mass: every state reachable from position 1154455 is excluded there."*
  chr7 at chain position 1,154,456; chr10 at 676,991 (pos 79,896,712). The other 20
  autosomes completed, and the same two chromosomes succeed at EVERY `eta` up to
  6.6566e7 — so it is the top point alone.
  **[Cause]** NOT the structural untraversable-block case decision K excluded. At a
  step with zero genetic distance the transition is the identity, so no mass mixes
  between donors. At low `eta` the forward message holds mass on many donors and
  some are callable at the next position. At extreme `eta` the posterior
  concentrates so hard that the donors carrying non-negligible mass are ALL
  excluded at that one position, and the total goes to exactly zero. "Reachable" in
  the message means "still carrying representable mass", which is `eta`-dependent.
  **[Proven fix / disposition]** NOT fixed, deliberately, and the curve is reported
  only up to `eta` = 6.6566e7 — the highest `eta` all 22 autosomes reach. Removing
  the limit needs a LOG-SPACE forward pass, and `ffbs.forward` is load-bearing for
  every release and every gate including the signed pre-registration's toy
  enumeration, so it is not changed on the side of an exploratory sweep. ⚠ The
  genome driver's completeness assertion REFUSED to print a 21-chromosome curve as
  the genome-wide result, which is that assertion working as intended.
  **[How to recognise it again]** A "zero total mass" raise that appears only at
  high `eta`, only on some chromosomes, and goes away when `eta` is lowered. If it
  fires at an `eta` you need, the answer is a log-space forward pass, not a wider
  tolerance. ⚠ Any sweep over a parameter must NAME THAT PARAMETER in its failure:
  the first traceback did not, and whether the curve died at its second point or
  its eighth is the only thing that decides what is reportable. See
  `docs/memory.md` 2026-10-23 (later still).


### A variance route reports a constant floor instead of a collapsing variance, in the flattering direction
- **[2026-10-23] `haplotype_variance` cancels below ~1e-33**
  **[Symptom]** `scripts/eta_sweep.py --validate` printed `var_q_u` =
  **7.703720e-34** at toy `eta` = 2000 AND at `eta` = 4400, while enumeration gave
  3.522099e-127 and 4.424217e-276. Identical to seven digits at two very different
  `eta`, which is the tell: a real variance does not stop moving.
  **[Cause]** The centred-moment recursion in `scripts/exact_variance.py`
  subtracts quantities of order 1 to obtain a variance of order 1e-130, so the
  result is pure cancellation noise once the posterior concentrates onto a single
  path. Nothing errors; the number is finite and plausible.
  **[Why it matters more than its size]** `Var_Q(u)` IS the slope of the response
  curve (`d E_Q[u] / d eta = Var_Q(u)`), so an overstated floor **overstates how
  much the next unit of `eta` would buy** — the flattering direction. At `eta` =
  4400 it overstates by 10^242.
  **[Proven fix]** Do not quote `var_q_u` below ~1e-28. The gate now compares the
  mean at every `eta` (agrees to 4.76e-16) and the variance only above that probe,
  prints the floored rows explicitly with the overstatement in powers of ten, and
  asserts the floor runs ABOVE the truth so the direction of the error is known
  rather than assumed. ⚠ `e_q_u` is NOT affected: it is a linear functional of the
  marginals with no centring and no cancellation. See `docs/memory.md` 2026-10-23
  (later).

### A toy gate compares nan to nan and reports agreement
- **[2026-10-23] The enumerator overflows before the sampler does**
  **[Symptom]** The first `eta` sweep gate asked for toy `eta` = 1e4 and 1e6, and
  printed `nan | nan | nan` for the enumerated value, the route and the relative
  difference — then `worst relative difference: 1.311e-14`, a PASS, because `nan`
  never exceeded the tolerance.
  **[Cause]** TWO different walls, and the lower one was not the one being guarded.
  `exact.release_distribution` computes `prob * exp(eta * u)` and `u` reaches 1 by
  construction, so it overflows at `eta` ~ **709 whatever `beta` is**;
  `utility.tilt_potentials` computes `exp(eta * beta_v)` and survives to
  709/max(`beta_v`), which is 4,963 on the toy and 6.1e6 on chr21. The gate had
  been written against the sampler's wall and the REFERENCE died first.
  **[Proven fix]** Compute the reference with log-sum-exp, local to the gate, which
  removes its wall entirely and is the better reference anyway because it owes
  nothing to the code under test. ⚠ `src/mechanism/exact.py` is deliberately NOT
  changed: it is load-bearing for the signed pre-registration's toy enumeration
  gates. The gate now also asserts every compared value is finite before comparing,
  so a nan can never again read as agreement, and asserts that past the wall the
  potentials DO go non-finite so the sweep's own finiteness guard has a real defect
  to catch.
  **[How to recognise it again]** Any gate whose PASS line sits under a row of
  `nan`. A tolerance comparison against nan is false, so `max()` silently ignores
  it.


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

### Every chromosome of a genome release drew the same uniforms
- **[2026-10-12] the RNG stream was keyed by `seed` alone, not by the chromosome**
  **[Symptom]** None visible. The release records would look correct, each part well
  formed, `log_z_p` in range, shares summing to 1. ⚠ Caught by reading the seeding path
  while building arm 1, not by any check.
  **[Cause]** `GenomeRelease` requires all 22 parts to agree on `seed`, because one genome
  release is one draw of one person — correctly. But `release()` derived its streams as
  `SeedSequence(seed).spawn(2)`, which depends on `seed` ALONE, so every chromosome of one
  genome release got bit-identical streams (verified). `_pick` samples by inverse CDF on
  one uniform per position, so every chromosome would have chosen its donor at the SAME
  QUANTILE of its own posterior at every step — a systematic cross-chromosome coupling
  that nothing in the record would reveal.
  **[Proven fix]** `release()` takes `stream_key`, and `ffbs._streams` derives
  `SeedSequence([seed, sha256(stream_key)]).spawn(n)`. The arm 1 runner passes the
  chromosome. The RECORDED `seed` stays the genome's, so `GenomeRelease`'s one-seed
  invariant is untouched and a part replays from `(seed, chrom)` — both already fields of
  `Release`, so no record schema changed. `stream_key=""` reproduces the old derivation
  EXACTLY, pinned by a test, so every single-chromosome release on disk stays
  reproducible.
  **[Severity]** Would have been invisible and uncorrectable after the fact. ⚠ **Nothing
  on record is affected**: no genome-scope release existed when it was found — the two
  records on disk are single-chromosome, where no cross-chromosome sharing can occur.
  **[Severity note]** FOURTH common-random-numbers defect in this project, and the FIRST
  caught before it reached an artifact. The previous three were found in the results.
  **[Found by]** Asking what `seed` keys, before launching a 308-release sweep.

### A parallel run's wall clock extrapolated from an isolated measurement
- **[2026-10-12] 2.5 h predicted, 8.7 h actual, on a 24-worker sweep**
  **[Symptom]** I timed ONE genome-wide release alone on an idle machine at 12.3 min,
  multiplied by the 308 cells, divided by 24 workers, and told Dylan ~2.5 h. Measured at
  51.7 min the workers were all on chr3 — 19.8% of one target — implying **8.7 h**, wrong
  by ~3.5x.
  **[Cause]** The isolated measurement had 48 cores and the full memory bandwidth to
  itself. These are O(T*K) streaming passes over GB-sized arrays, so with 24 workers the
  bottleneck is MEMORY BANDWIDTH, not cores: each worker ran **4.1x slower** than in
  isolation, giving ~6x aggregate throughput and about **25% parallel efficiency**.
  Adding workers cannot fix a bandwidth cap.
  **[Proven fix]** No code fix — an estimation rule. **Time a parallel run from the
  parallel run**, not from an isolated sample: per-worker speed under contention is the
  only number that extrapolates. ⚠ And instrument long runs so progress is observable:
  this one wrote nothing for its first 4.4 h because each target saved all 7 `tau` at the
  end, so the only way to locate the workers was reading `/proc/<pid>/maps` for which
  chromosome each had mapped. `scripts/arm1_release.py` now prints per-chromosome.
  **[Severity]** No wrong number was published, but it drove a planning decision: at the
  true cost, computing the releases and the primary endpoint in SEPARATE passes was
  ~17.4 h against ~11.5 h sharing one forward pass. The correction is what surfaced the
  saving.
  **[Found by]** Asking where the workers actually were, instead of trusting that 0/44
  after 51 min matched the prediction.

### A normalisation gate whose tolerance is computed from the thing it checks
- **[2026-10-12] `sum beta_v` was unchecked on the whole release path**
  **[Symptom]** A hand-built `beta` summing to 2.0 runs the entire release path with no
  exception and no symptom. The release is stamped `tau` = 0.5 while really carrying
  TV = tanh(`eta_tau` * 2) = 0.80.
  **[Cause]** Two gaps composing. `beta_v` asserts `sum beta_v` == `share`, but nothing
  requires a caller to route through `beta_v`, and `scripts/beta_concentration.py`
  legitimately hand-builds allocations. The one downstream check that looks protective is
  not: `ffbs._draw_one` passes `log_z_p` a bracket of `eta * float(beta.sum())`, so the
  tolerance is DERIVED FROM THE QUANTITY UNDER TEST and widens in lockstep with the error.
  `tilt_potentials` validated only array lengths and `eta_tau` >= 0.
  **[Proven fix]** `tilt_potentials` rejects `sum beta_v` > 1 + 1e-9 — the one function
  every tilted release must call. `<=` and not `==`, because a per-chromosome release
  carries that chromosome's `beta_share` and the 22 shares sum to 1 (decision J).
  Mutation-verified: removing it lets `sum beta_v` = 2.0 through.
  **[Severity]** Silent and unbounded in the dangerous direction — it inflates the real
  `tau` above the stamped one, and `p_succ_bound` cannot catch it because measured TV
  reaches only 27-43% of `tau`.

### Two figures pinned to the superseded coordinate-cut chain
- **[2026-10-12] `u_ceiling` and the float64 headroom argument both quote T = 305,886**
  **[Symptom]** `u_ceiling` on HG00096 chr21 at scope `chain` measures **0.692509** on the
  live build, against **0.815455** on record in three places. Separately,
  `src/mechanism/ffbs.py`'s float64 justification quotes a per-position tilt of 1.80e-06
  and "15.1x the noise floor".
  **[Cause]** Both were measured on the pre-decision-C coordinate-cut chain, T = 305,886.
  The live decision-K chain is T = 332,521. The rebuild added 26,635 positions of which
  only 4,362 are called for HG00096, so the called fraction fell 0.829871 -> 0.776516 and
  the inverse-support-weighted ceiling fell further still. Restricting to the old span
  reproduces 305,886 and 0.815455 exactly, which pins the cause to the build, not the
  arithmetic.
  **[Proven fix]** Corrected in place at `src/mechanism/ffbs.py` (both the docstring and
  the runtime error message) with the live T, the `eta_tau`/2 halving, and the
  genome-wide figure. ⚠⚠ Correcting it STRENGTHENS the guard: at T = 22,134,562 the
  realised per-haplotype tilt at `tau` = 0.5 is 1.2408e-08 = **0.104x float32 epsilon**,
  an order of magnitude BELOW the representable step, so in float32 the tilt vanishes
  outright rather than losing precision.
  **[Severity]** Documentation-only for the headroom claim, and in the conservative
  direction. For `u_ceiling` it is a live denominator: anything divided by 0.815455 on
  this build overstates by 17.8%.
  **[Found by]** The first row of a new experiment failing to reproduce a recorded number,
  which is only diagnosable because the record carried its DENOMINATOR (253,846 of
  305,886) rather than the ratio alone.

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

- **[Symptom]** A whole-tree checker dies with `FileNotFoundError` on the first
  chromosome that lacks an artifact, so the other 21 are never checked at all.
  **[Cause]** `scripts/verify_cohorts.py` loaded every per-site array unconditionally.
  chr21 was built BEFORE decision C, so it has no `cm_source.npy` — the coordinate
  cut had nothing to record — and the loader raised before any check ran. The same
  artifact was then read a SECOND time in `main()` to build the results table, so
  fixing only the first reader reproduced an identical traceback.
  **[Proven fix]** A missing artifact is now a reported FINDING naming the
  chromosome, the file and the manifest's `chain_rule`, so a pre-decision-C tree
  says so instead of crashing. Both readers are guarded; `tests/test_verify_cohorts.py`
  asserts the report, not an exception.
  **[Severity]** Exactly the "found nothing vs ran nothing" confusion the fan-out
  rule exists for — one stale chromosome suppressed the verification of all 22.

- **[Symptom]** A checker reports a length mismatch and then raises
  `ValueError: operands could not be broadcast together` on the next check.
  **[Cause]** Every check after the length comparison is POSITIONAL — it compares
  arrays element by element — so once two of them have different lengths the rest
  are not merely uninformative, they raise. The length check appended its problem
  and fell through. Found by the mutation test, not by use: truncating
  `support.npy` by one row is a corruption nothing in the real tree had yet.
  **[Severity]** Would have turned the single most likely real corruption — a
  per-site array off by one row — into a crash rather than the finding it is.
  **[Proven fix]** `check()` returns immediately on any length disagreement, and
  `test_length_mismatch_returns_before_positional_checks` pins that the report is
  the ONLY problem returned, so a later check cannot creep back in front of it.

- **[Symptom]** `FileNotFoundError: data/temp/chrN/chrN.hprc-v1.0-pggb.gfa` raised
  from `store.sha256_file`, AFTER `vg deconstruct` had already run for one to two
  hours. Six chromosomes (chr10, chr11, chr12, chr13, chr14, chr17) lost in one
  batch; chr15, chr21 and chr22 were in flight and hit the same wall.
  **[Cause]** `src/cohort/vcf.py:deconstruct` unlinked the decompressed GFA and then
  passed the SAME PATH to `store.sha256_file` three lines later, inside
  `store.attach_prov`. Introduced by commit `8d56346`, whose purpose was unrelated
  (reclaiming 4-9 GB per chromosome); the cleanup was inserted above the provenance
  call rather than below it.
  **[Why the suite could not see it]** Every pre-existing test passed a PLAIN
  `.gfa`, so `decompressed_here` was False, the unlink never fired, and the ordering
  was unreachable. The failure needs the `.gz` path AND real provenance AND the
  cleanup, together. ⚠ And the traceback named `store.sha256_file`, so it read as a
  provenance bug when the cause was three lines above in cleanup code added for a
  different reason.
  **[Why it hit mid-run]** ⚠⚠ **THE COMMIT LANDED WHILE THE 22-CHROMOSOME BUILD WAS
  RUNNING.** Each chromosome runs in a fresh subprocess that imports `vcf.py` at
  launch, so chromosomes whose process STARTED before 15:37 ran the old code and
  succeeded — chr9's manifest is timestamped 16:09, 32 minutes AFTER the commit, and
  it is fine — while every chromosome started after it failed. The split is by
  process start time, not by wall clock, which is why 11 chromosomes succeeded and
  the next 9 did not.
  **[Proven fix]** `gfa_digest = store.sha256_file(gfa)` is taken BEFORE the unlink
  and the stored value is passed to `attach_prov`. `tests/test_deconstruct_prov.py`
  asserts the order by its observable consequence — a sidecar carrying the digest of
  the decompressed bytes — rather than by "no exception raised"; reinstating the
  original ordering fails 3 of its 6 tests. A control test asserts the cleanup still
  happens, so the ordering test cannot pass trivially on a function that never
  deletes.
  **[Recovery, and why it is cheap]** vg writes the VCF BEFORE the crash, so all
  nine VCFs survived complete (543 MB - 2.15 GB each). Re-running vg to recover a
  small JSON sidecar would cost about 6 min of wall time and 2 CPU-hours per
  chromosome, so
  `vcf.attach_deconstruct_prov` attaches the sidecar to an existing VCF without
  invoking vg, recomputing the GFA digest with `store.sha256_gz_plain` — which
  streams the `.gz`'s decompressed bytes and writes nothing, since gzip
  decompression is deterministic. It REFUSES to overwrite an existing sidecar (real
  provenance must not be silently replaced by a reconstruction) and stamps
  `prov_recovered_after_crash: true` so a reconstructed sidecar says so.
  **[Severity]** Loud, not silent — `run_jobs` raises and the driver exits non-zero,
  which is the one thing that went right. But it cost 9 chromosomes' worth of
  scheduling and would have cost about 12 CPU-HOURS of vg (roughly an hour of wall
  time at this parallelism) without the recovery path.
  ⚠ **[CORRECTION - 2026-10-10]** An earlier version of this entry said "1-2 hours
  per chromosome", which was WALL time read off `ps -o time=` -- a CUMULATIVE CPU
  figure summed over 32 threads, not elapsed time. The measured elapsed cost of a
  whole chromosome from scratch is **437.9 s (chr18) and 449.9 s (chr20)**, running
  3-wide. `ps -o time=` is CPU time; `ps -o etime=` is elapsed. Same instrument
  error class as the four false findings recorded earlier today.

- **[Symptom]** The new cohort verifier reported `chr9: positions are not strictly
  increasing (1 violations)` on a chromosome that is correct.
  **[Cause]** MY CHECK, not the data. It demanded a STRICT increase. chr9 has exactly
  one tied position, 37,752,316, where two genuinely distinct records share a start:
  a parent bubble (`lv` = 1, 302 bp REF, 4 ALTs) and a 1 bp deletion NESTED inside it
  (`lv` = 2, GG -> G). Nesting is routine rather than anomalous — 10-18% of records on
  every autosome carry `lv` > 0 — so a shared start coordinate is legal
  `vg deconstruct` output. Every other built chromosome happens to have zero ties,
  which is why only chr9 tripped it.
  **[Severity]** A FALSE POSITIVE on real data, and the expensive kind: the natural
  next action is to rebuild or discard a chromosome that was fine. Same class as the
  four false findings recorded on 2026-10-10 — a check that fails for the wrong
  reason is worse than one that cannot fail.
  **[Proven fix]** The check now forbids a DECREASE, which is the error that actually
  matters (it makes a genetic distance negative and puts the Li--Stephens switch
  probability outside [0,1], and `numpy.interp` would not complain). Ties are counted
  and REPORTED as a `shared pos` column rather than failed, so a jump in the count is
  still visible. `tests/test_verify_cohorts.py` has both directions: a decrease must
  be caught, and a tie must NOT be reported — the second fails if anyone tightens it
  back.

- **[Symptom]** Re-running a chromosome after a late-pipeline failure redid
  `vg deconstruct` from scratch -- about 6 min of wall time and 2 CPU-hours -- although
  a complete VCF was on disk.
  **[Cause]** `src/cli.py`'s resume map listed `arrays`, `support`, `genmap` and
  `cohort-vcf` but NOT `deconstruct` — the single most expensive step. So any crash
  anywhere later in the pipeline cost the whole chromosome rather than the step.
  **[Proven fix]** `deconstruct` is in the map, and the skip branch REPAIRS missing
  provenance instead of carrying the gap forward: a complete VCF with no sidecar gets
  one from `vcf.attach_deconstruct_prov`, which needs no vg. Verified on real data —
  chr10 recovered in 9 s against a 1-2 h rebuild, and `store.sha256_gz_plain`
  reproduced the digest the ORIGINAL code recorded, byte for byte, on all three
  chromosomes where both exist (chr9, chr16, chr19). ⚠ The map now also requires a
  NON-EMPTY output, because a zero-byte file is what a killed step leaves behind and
  treating it as complete feeds a truncated artifact to the next step.

- **[Symptom]** The Arm 2 attack harness would have reported a privacy VIOLATION with
  Theorem 1 fully intact, and its own comment said such a measurement "falsifies the
  implementation" — so the next action would have been hunting a bug that does not exist.
  **[Cause]** `p_succ_bound` = (1 + `tau`)/2 bounds the EQUAL-PRIOR accuracy
  `0.5 * P_p(Lambda > 0) + 0.5 * P_q(Lambda < 0)`, which needs draws from BOTH input
  paths. `scripts/arm2_attack.py` drew every release from `Q_p` and compared the
  ONE-SIDED rate `P_p(Lambda > 0)` to that bound. From
  `balanced = 0.5*TPR + 0.5*(1 - FPR) <= (1 + TV)/2` we get `TPR <= TV + FPR`, so an
  attacker fully respecting the theorem with FPR 0.5 reaches `TPR = tau + 0.5`, which
  EXCEEDS (1 + `tau`)/2 at every `tau` > 0; and a degenerate attacker that always
  answers "p" scores **1.0** one-sided while sitting at exactly chance balanced.
  **[Root cause, and the generalisable part]** `src/attacks.py` ALREADY had a correct
  balanced `empirical_accuracy`. The script reimplemented the accuracy inline, wrongly.
  ⚠ A statistic defined twice will drift; the second definition is where the error
  lives, and it is invisible because the first one is right and tested.
  **[Why 265 tests could not see it]** Every test exercised the correct helper. Nothing
  tested the SCRIPT's arithmetic, and nothing asserted the relationship BETWEEN the
  one-sided and balanced quantities.
  **[Proven fix]** The script calls `src/attacks.py:empirical_accuracy` and asserts it
  agrees with `0.5*TPR + 0.5*TNR`. The trial draws one release from `Q_q` per condition
  with an independent seed. `--one-sided` survives as a quarter-cost diagnostic that
  REFUSES to compare anything to `p_succ_bound`, printing why instead of a number. Two
  new tests encode the derivation and both fail if `empirical_accuracy` is made
  one-sided.
  **[Severity]** The worst class available: a false POSITIVE on a privacy claim. It
  would not have been caught by inspection, because the number would have looked like a
  finding.

- **[Symptom]** Four stale claims in `docs/terms.md` — including "there is no pinned
  attacker set yet for either arm" when four are pinned, and `src/attacks.py` marked
  `PLANNED` three days after it was written — survived a passing 21-test docs contract.
  **[Cause]** The contract's path check skips any line carrying an absence marker, on
  the reasoning that such a line is honestly declaring a gap. It had NO COUNTERPART:
  nothing ever verified that a declared gap is still a gap. So a line could claim a file
  was unbuilt forever and the suite would agree by looking away. The `<date>`-style
  placeholder was doubly invisible, because the check also skips anything containing
  `<` or `>`.
  ⚠⚠ **AND THIS WAS A KNOWN HOLE.** `docs/terms.md` recorded it on 2026-10-08 — "the
  row survived only because the contract test's `_ABSENT_MARKERS` skips any line
  containing PLANNED" — and it was left open. A known blind spot left open is worse than
  an unknown one, because the record creates the impression it is handled.
  **[Proven fix]** Two converse checks: a path claimed absent must really be absent, and
  a `<placeholder>` attached to an absence claim must not resolve to a real file. Both
  are mutation-tested and both fail on their own injected claim with a clean control.
  ⚠ **The escape for corrected text had to be made POSITIONAL.** A `docs/data.md` row
  runs to thousands of characters and can carry a `[CORRECTION]` about one claim while
  making a stale claim about another; a line-wide escape swallowed exactly that and the
  gate stayed SILENT on a verified mutation — decoration, by this project's own
  standard. An escape now excuses only what follows it, which is also how the in-place
  correction rule reads: the marker introduces the superseded text it quotes. ⚠ The
  escape is keyed to the ABSENCE MARKER's position, not the path's — the path is merely
  what the claim is about and usually appears earlier on the row.

- **[Symptom]** ⚠⚠ **BLOCKER, FOUND 2026-10-10 BY ACTUALLY RUNNING THE RELEASE.** The
  forward pass raises on the REBUILT chr21: "chain position 29210 (pos 10,543,014) has
  zero total mass: every state reachable from position 29209 is excluded there."
  Nothing can be released genome-wide until it is settled.
  **[Cause]** NOT a code defect. Where two chain positions have ZERO genetic distance
  the transition is the IDENTITY, so a donor that becomes unavailable at the next
  position has nowhere to go; when all surviving mass sat on the donors that disappear,
  the total is zero and the path genuinely dead-ends. The model is internally
  consistent — `docs/terms.md` already states P(switch) = 0 on a zero-cM plateau is
  CORRECT. It is the DATA that makes the state space disconnected there.
  **[Why it is new]** ⚠ These are the **interior zero-cM plateaus, which decision C
  never addressed.** Decision C's scope was map-SILENT regions, and it is clean there:
  all 9,384 of chr21's zero-cM steps carry `cm_source` = 0, i.e. an INTERPOLATED cM
  from the map, not a declared-rate one. But retiring the coordinate cut also admitted
  the plateau positions BELOW the old span's lower bound, which nothing ever decided to
  include — chr21 has 397 zero-cM steps with a shrinking donor set and **395 of them
  are outside the old cut, all below its lower bound**. The old build had only 2 such
  steps inside its span, which is why every chr21 release to date worked.
  **[Scale — measured, all 22 autosomes]** **357,918 zero-cM steps genome-wide, 9,971
  of them with a shrinking donor set**, worst single drop **87 of 88 donors**. It is not
  a chr21 quirk and not an acrocentric quirk: chr9 is the worst at 3,366 and chr20 has
  962. ⚠ A step only RAISES when the surviving mass is entirely on the vanishing
  donors, so 9,971 is an upper bound on the failures, not a count of them — but one is
  enough to stop a chromosome.
  **[Status]** ✅ **RESOLVED the same day as DECISION K (Dylan, 2026-10-10):** the whole
  block leaves the chain. 20 blocks, 172,142 positions, 0.772% of the chain; the release
  then ran. `src/cohort/genmap.py:untraversable_blocks` implements it,
  `tests/test_untraversable.py` has 15 tests and three mutations are caught, and
  `scripts/verify_cohorts.py` recomputes the blocks INDEPENDENTLY so a chain that was
  never rebuilt fails there rather than at release time.
  ⚠ Two holes in my own checks surfaced doing it, both found by the opposite-direction
  test: the verifier asserted `in_chain` was exactly `support > 0` (falsified by K), and
  its replacement only recomputed the blocks when something was ALREADY excluded — so a
  chain never rebuilt skipped the branch and passed, which is the one case it exists for.
  Any fix must stay target-independent (standing constraint 1); all three candidates
  are, since they depend only on the map and on support(v).
  **[Severity]** Loud, not silent — it raises rather than returning a path, which is
  the one thing that went right. ⚠ Found only by RUNNING the release end-to-end; 290
  tests pass and none of them touches a real genome-scale chain.

### Common random numbers, THE THIRD OCCURRENCE — and I diagnosed it backwards first
- **[2026-10-11] The Arm 2 attack reused 25 seeds across all 70 (target, `tau`) cells**
  **[Symptom]** `E_p[Lambda]` came out at **-1.46e-04, z = -3.81**, where theory
  requires `E_p[Lambda] = KL(Q_p||Q_q) >= 0`. All ten targets negative. I reported it
  as an implementation defect and stopped the chunk.
  **[Cause]** NOT an implementation defect. `scripts/arm2_attack.py` built its job list
  as `(t, tau, 50_000 + i, ...)`, seeding from the DRAW INDEX ALONE, so all 10 targets
  and all 7 `tau` reused seeds 50000..50024. With a per-position tilt of
  `eta * mean(beta)` = 0.2554 / 332,521 = **7.7e-07** the tilt cannot move the draw:
  **the released path is BIT-IDENTICAL at `tau` = 0, 0.05 and 0.25**, `u_p` agreeing to
  nine digits. The seven `tau` columns were ONE experiment reported seven times, the
  ten targets were near-replicates (cross-target correlation **+0.305**), and the three
  decoy conditions share each `Q_p` release (pairwise **+0.57 to +0.63**).
  Under clustered standard errors — analytic, bootstrap and a t-test on seed means, all
  agreeing — z goes from **-3.81 to -1.99** (p = 0.058). There was no anomaly.
  **[What the implementation actually does]** ✅ The exact, sampling-free route gives
  `E_Q[Lambda]` = **+1.784e-07** for HG00096 and POSITIVE for all ten, equal to KL as
  theory requires. The sampler is unbiased: 3,320 fresh draws sit at **z = -0.16**
  against that exact value, and the shift needed to produce the reported mean is
  **15.0 SE** away. `log_z_p` is bit-equal to the recorded value. Instrument validated
  four ways, including against the repo's own enumerator to 1e-16.
  **[⚠⚠ THIS IS THE THIRD OCCURRENCE OF A TWICE-DOCUMENTED FAILURE]** The entry "A
  headline finding measured under common random numbers" (2026-09-30) describes the
  same mechanism, and `docs/memory.md` records seed reuse correlating two chromosomes
  of one release. I wrote the script that reproduced it AFTER both were logged. The
  lesson that did not transfer is not "watch for seed reuse" — it is that **a tilt of
  1e-06 per position means ANY two runs sharing a seed are the same draw**, so seed
  structure is the FIRST thing to check whenever this mechanism's outputs look alike or
  move together.
  **[⚠ AND I DIAGNOSED IT BACKWARDS]** Three claims of mine were wrong, and the middle
  one actively misdirected the search:
    1. "These numbers can't be right" — the (1+`tau`)/2 floor bounds the TRUE optimal
       accuracy, not a 250-sample estimate of it. Estimator confused with estimand,
       for the second time in one session.
    2. "The bias is EXACTLY first order in `eta`, which rules out variance terms" — an
       ARTEFACT presented as a fingerprint. The sampled `(u_p - u_q)` is literally one
       number at every `tau` and `(log Z_p - log Z_q)/eta` is constant to 8 figures, so
       their ratio is `tau`-independent BY CONSTRUCTION. It ruled out nothing, and I
       used it to narrow the hypothesis space.
    3. "Ten of ten negative is 1-in-1024" — void. Their z-values have sd **0.643**, not
       the 1.0 of ten independent tests; per seed the negative count runs 4/10 to 10/10
       and the mean flips sign.
  **[Proven fix]** A `SeedSequence` over `(base, target index, tau index, draw index)`,
  giving 1,750 distinct streams over the 1,750 cells, verified distinct. ⚠ This
  RESTORES the pre-registered design rather than changing it: `delta` = 0.0556 assumed
  250 independent trials per cell and clustering had cut the effective n to about 67
  (`delta` ~ 0.107). Four tests; both defects caught on reintroduction.
  **[Severity]** The null was the right answer for the wrong reason, and the error
  flattered us — the bias pushed the attacker BELOW chance. ⚠ And no sample size would
  have fixed the diagnostic I chose: `sd(Lambda)` per draw is 6.15e-04 against a true
  value of 1.78e-07, so resolving `E[Lambda]` by sampling needs n ~ 1.2e7. The exact
  marginals route gives it for free, and A2 already pre-registers an exact route.

- **[2026-10-11] Decoy selection was not reproducible: Python's `hash()` is salted**
  **[Symptom]** None visible. The recorded run looked internally consistent.
  **[Cause]** `rng = np.random.default_rng([DECOY_SEED, hash(target_id) % (1 << 31)])`.
  `PYTHONHASHSEED` randomises `str.__hash__` per interpreter. MEASURED:
  `hash('HG00096') % 2**31` returned **1730677557, 783898805, 1043054338** on three
  invocations of the same binary. The run was self-consistent only because `fork`
  inherits the parent's hash seed, so every worker agreed with its parent — and a
  re-run picks a different D3 decoy.
  **[Severity]** ⚠ A pre-registered experiment whose decoy set cannot be reproduced
  from its recorded seed is not pre-registered in the sense standing constraint 6
  means. The cell was irreproducible while looking fully determined.
  **[Proven fix]** `hashlib.sha256(target_id.encode())`, stable across interpreters,
  machines and versions. ⚠ The test SPAWNS THREE FRESH INTERPRETERS and compares,
  rather than asserting stability in-process — an in-process check cannot detect a
  per-process salt, which is exactly why this survived.

### A checksum computed from the download verifies STORAGE, not TRANSFER
- **[2026-10-11] chr10's panel arrived corrupt, and every check we had passed it**
  **[Symptom]** Ten of 220 target-preprocessing jobs failed, all on chr10, all with
  `bcftools view exited 9`. The real message was `[E::bgzf_uncompress] CRC32 checksum
  mismatch`, in the block at offset **393,732,438** of a 1,493,904,721-byte file.
  **[Cause]** A corrupt TRANSFER. `scripts/get_panel.sh` did
  `curl ... -o "$f.part" && mv "$f.part" "$f" && sha256sum "$f" > "$f.sha256"`, so the
  sidecar was computed FROM the downloaded bytes. ⚠⚠ **It therefore matches whatever
  arrived.** VERIFIED on the corrupt file: recorded and actual sha256 are
  `3b54d73d6a2b...` and `3b54d73d6a2b...` — identical. A self-computed checksum detects
  later bit-rot ON DISK; it cannot detect a bad download.
  **[Why every other check passed]** The corruption is MID-FILE, and
  `scripts/verify_panel.py` deliberately avoided reading whole files: the index was
  intact, the record count right (3,495,473), the sample set right (3,202), the
  reference-build comparison against our own pangenome passed, and the truncation test
  SEEKS to the last 3 Mb — which never touches the middle. It reported chr10 sound.
  ⚠ And I had told Dylan `--deep` "would have caught this". It would NOT: `--deep`
  compared the file to that same self-referential sidecar.
  **[Proven fix, in two places]** Every bgzf block carries its own CRC, so
  decompressing the whole file is a check that can FAIL ON ITS OWN TERMS rather than
  against a reference derived from the same bytes.
    * `scripts/get_panel.sh` runs `bgzip -t` on the `.part` BEFORE renaming it into
      place and DISCARDS a failed one. ⚠ Discarding matters: `curl -C -` resumes from
      the existing size, so a complete-but-corrupt file would be skipped forever.
    * `scripts/verify_panel.py --deep` runs `bgzip -t`; the sha256 check is kept but
      relabelled for what it can actually prove — that the bytes have not changed
      SINCE we recorded them.
  **[Swept]** All 22 panels tested with `bgzip -t`: **only chr10 corrupt, 21 pass.**
  **[Severity]** ⚠ It failed LOUDLY three steps downstream instead of at fetch time,
  which cost 10 of 220 jobs rather than a silent wrong answer — `bcftools` refuses to
  read a bad block rather than returning partial data. But the same corruption in a
  file read by something more forgiving would have produced a short result with no
  error at all. ⚠ The general lesson is the one in the title, and it generalises past
  downloads: a verification derived from the artifact it verifies proves only internal
  consistency. ⚠ It also sharpens the cohort-build cleanup rule recorded earlier --
  `store.sha256_file` sidecars across this project share the property, and are
  similarly storage checks rather than provenance checks about upstream content.
  **[Tests]** 12 in `tests/test_verify_panel.py`, including one that corrupts a real
  bgzf file by flipping a single mid-file byte and then asserts that its freshly
  computed sha256 MATCHES — the finding written as an assertion rather than prose.

### [Symptom] A "resolved" item and the stale text it resolves BOTH present in a document

Seen 2026-10-13 in `paper/draft_2026-10-13.md`: limitation 11 appeared twice — once as
`✅ RESOLVED` and once as the original `NOT COMPUTABLE BY THIS ROUTE / ~3.4e-11` text that the
resolution disproves by a factor of 13,000 — and limitation 12 appeared twice verbatim. A reader
following the list in order hits the CORRECTION first and the STALE CLAIM last, which is worse
than never having corrected it. Detected only because a later edit asserted `count(old) == 1` and
got 2.

[Proven fix] An insert-before-anchor edit CANNOT resolve anything. A resolution must REPLACE the
old text, and the edit must then assert the old distinctive phrase is ABSENT from the file, not
merely that the new one is present. Grep the distinctive NUMBER or PHRASE afterwards — here
`3.4e-11` and `NOT COMPUTABLE` — across `docs/`, `paper/`, `src/`, `tests/` and `CLAUDE.md`,
because the same stale conclusion was simultaneously live in a `docs/data.md` row and in two
places in `docs/plan.md`. Corrections do not propagate themselves.

### [Symptom] A pre-registered positive control passes at 100% on the first try

Seen 2026-10-13 on attack A3. The pre-registration requires the attack to "link the
UNOBFUSCATED target at near-100%", and it did -- 30 of 30. ⚠ The control was
ALGEBRAICALLY UNABLE TO FAIL. PanMixer's linkage score sums `-log P_HWE` only over sites
where the genotypes AGREE, mismatches contribute nothing, and the weights are
non-negative, so `score(y, q)` over the agreeing sites is a SUBSET SUM of `score(y, y)`
and `Lambda > 0` is forced the moment a single site differs. A 100% pass rate was the
only possible outcome.

[Proven fix] Before trusting a control, ask what INPUT would make it fail and construct
it; if none exists, the control is decoration and must be replaced, not reported. Here
the replacement is a graded POWER CURVE -- an untilted release with a declared fraction
`f` of positions overwritten by the target verbatim -- which answers the question the
control was meant to answer: A3 needs f = 0.01 (about 3,300 of 332,521 positions) to
reach 0.75 accuracy and is at chance at f = 0.001. ⚠ A perfect score on a first run is
evidence about the TEST at least as often as about the thing tested; the same reasoning
caught the `log_z_p` gate that could not fail in 2026-10.

### [Symptom] A variance comes back `nan` for a whole chromosome while the toy gate passes

Seen 2026-10-22 on `scripts/exact_variance.py`'s first run: `var_q_u`, `var_q_differ`
and the centred-mean gate all printed `nan` for every one of chr21's 332,521 positions,
while `e_q_u` in the same run reproduced its logged value to five significant figures —
so the setup was right and only the new recursion was poisoned. The enumeration gate on
the project's toy cohort had passed at 7.7e-15 immediately before.

The cause is a 0/0 in the backward kernel's normaliser,
`D_t(i) = s_t/K + (1 - s_t) * alpha_{t-1}(i)`. It needs TWO conditions at once: a step
with `switch_mass == 0` — zero genetic distance, where `unif = s/K` vanishes and the
transition is the identity — AND a donor with no call at the SOURCE position, so
`alpha_{t-1}(i) == 0` as well. chr21 has **1,058** zero-distance steps and **7** of them
sit beside an unavailable donor. One NaN then propagates through every later scalar
product, because `0 * nan = nan`, so a defect at 7 positions destroys all 332,521.

⚠ **THE TOY COULD NOT CATCH IT, AND THAT IS THE REAL LESSON.** `toy.toy_cohort()` has a
zero-distance step (index 2, `switch_mass` 0) and it has uncalled cells (2 of them), but
never at the same place — so `D` never reaches 0 there. A gate built from the richest
fixture in the project still had a hole exactly where the real chain does not.

**[Proven fix]** Guard the division, and ASSERT the reason it is safe rather than
assuming it: where `D == 0`, check that `alpha_t(i) == 0` too — with the identity
transition such a state is reachable only from itself at `t-1`, where it has no mass, so
it is genuinely unreachable and its conditional moments are weighted by zero everywhere
they are used — then raise if any such state does carry forward mass, and otherwise set
its moments to 0 so the dot products stay finite. Then build a fixture where the two
conditions COINCIDE: `validate_dead_kernel()` uses a 4x5 cohort with donor 3 uncalled at
position 1 and zero genetic distance from 1 to 2, agrees with enumeration to 5.2e-15,
and returns `(nan, nan)` without the guard — measured, so the gate is known to fire.
⚠ When a new gate passes on the existing fixture, ask which real-data property the
fixture lacks before trusting it; "the toy has feature A and feature B" is not "the toy
has A and B together".

## [A rounded display made four different numbers look identical, and the claim was written down as an identity] -> [State what the test tests; compare the objects, not their formatted form]

**Symptom.** `docs/data.md` and `docs/memory.md` both asserted that PanMixer's capacity
plateau is FOUR byte-identical settings (0.3/0.5/0.7/0.9). It is THREE. Measured 2026-10-28
by md5 over all 44 subjects' `new_haplotypes.npy`: 0.5/0.7/0.9 are byte-identical for 43 of
44 and 0.3 is DISTINCT (median 24 cells away); only HG03492 includes 0.3. The same holds on
PanMixer's OWN build, where the phase-A arrays survive: `cap0p3` hashes `fa2f721c377d`
against `5f3321ed2227` shared by the three above it.

**Why it survived.** `logs/pmsweep.log` printed
`changed=106534/649242 alleles (16.41%)` at capacity 0.3 and
`changed=106549/649242 alleles (16.41%)` at 0.5/0.7/0.9. Four lines, one percentage, two
different counts — **the percentage is rounded to two decimals and the counts are fifteen
alleles apart.** Reading the formatted column instead of the counts gave "four settings
identical". `scripts/fidelity_frontier.py:99` then made it durable: it computes the
plateau's start as every capacity within **0.01 PERCENTAGE POINTS** of the maximum and
printed "four settings return the identical count" — a TOLERANCE reported as an identity,
with the word "four" hand-written in the prose rather than computed. It never compared
bytes or counts. The claim was then copied forward into two doc files and a published page.

**Proven fix.** `scripts/fidelity_frontier.py` now computes `_plat` and prints how many
settings land within 0.01 points, says explicitly that this is NOT byte-identity and that
byte-identity is a separate test, and names its outcome. The corrections are in place in
`docs/memory.md` (2026-10-27 entry), `docs/data.md` (two rows) and `docs/plan.md`.

⚠ **The contradiction was in the same document.** The median table printed directly above
the false sentence in `docs/memory.md` reads 0.176759 at capacity 0.3 against 0.176790 at
0.5/0.7/0.9. The prose contradicted its own table, one line away. When writing a ✅ beside a
reproduction claim, re-read the table above it as the first check.
⚠ No reported figure moved: `arm1_headtohead.py` selects their column with `max(tf)` and
`pm_draw_noise.py` reads `draws[(s, 0.9)]`, both inside the real plateau. A documentation
defect, not a measurement defect — but it was quoted as a validation against their paper.

## [`local a=$1 b="${a}"` aborts under `set -u` with "a: unbound variable"] -> [One `local` per variable when a later one reads an earlier one]

**Symptom.** A four-job timing harness died immediately with
`line 8: name: unbound variable` four times, from
`local name=$1 exp=$2 row=$3 cwd="$SCR/timing_${name}_exp${exp}_row${row}"`.
Nothing ran, and with `set -u` absent it would instead have built paths from EMPTY strings
and silently written every job's output to the same file.

**Proven fix.** Split the declaration — `local name=$1` then `local exp=$2` and so on, with
the composed value last. Bash word-expands **every argument to `local` before the builtin
runs**, so `${name}` on the same line is still unset at expansion time; the left-to-right
assignment order people expect applies to the builtin's execution, not to the expansion.
⚠ This is the failure mode `set -u` exists for: loud and immediate instead of four jobs
overwriting one path. Keep `set -u` in every fan-out harness.

## ["Applied identically to both sides" is not "does not bias the comparison"] -> [Ask what property the identical rule keys on, and whether the two sides share it]

**Symptom.** Both of PanMixer's in-graph utility evaluators were run unmodified on our
releases and on theirs, against one build, with the same 44 people — and the resulting
headline figures were recorded, and committed, as a settled split: `af_loss` "they are
lower for 0 of 44, p = 1.14e-13" and `ld_loss` "we are lower for 40 of 44, ~7.5x". The
entry asserting them also asserted that the evaluators' shipped quirks, "applied to both
sides, so they do not bias the comparison". **That reasoning is wrong, and both headlines
move when it is tested** (`scripts/arm1_afld_sensitivity.py`,
`logs/B_arm1_afld_ourbuild/sensitivity.md`).

An identical rule has an **asymmetric** effect whenever it keys on a property the two
sides do not share. Both quirks do:

**`af_loss`** treats an uncalled cell as an ordinary allele, so a cell whose
called/uncalled status flips scores exactly like an allele swap. What that penalises is
divergence from the row being replaced, and the two mechanisms differ there **by
construction**: our release is a fixed-support path (uncalled 16,656..16,666, only 5
distinct values over 308 releases, correlation with the master row's own uncalled count
-0.32) while PanMixer edits a copy of the target's row (44 distinct values, +0.9991).
Genuine allele swaps are within 1% (111,615 against 112,922). Across three defensible
policies the verdict is 0/44, 1/44, then **31/44 the other way** (p = 0.00956) once the
missingness skeleton is held fixed. Nothing about either mechanism changes between them.

**`ld_loss`** computes `p_A` and `p_Bs` as `nanmean`s over the non-missing cells but
`p_AB = nanmean((xs == 1) & (ys == 1))` over a BOOLEAN array, which holds no NaN because
`np.nan == 1` is False — so it always divides by 88 and D mixes two denominators. Our
advantage is **1.30x, not 7.53x**; and the shipped figure's largest **6 pairs of ~7,500
carry 94%** of their total, so it is not a mean over 7,462 pairs in any useful sense.

**Proven fix.** `scripts/arm1_afld_sensitivity.py` measures each evaluator under its
shipped rule and under one counterfactual, paired over the same 44, and prints the
support-divergence table that explains the `af_loss` result. The shipped figures stay the
ones to publish — they are what a reader of PanMixer's paper reproduces — but every quote
carries its policy dependence.

⚠ **The test to apply, not the number to remember.** "Both sides get the same treatment"
licenses nothing on its own. Name the property the rule keys on, then measure whether the
two sides share it. Here they provably do not, and it was knowable before the fan-out ran.
⚠ **How it was caught, and how it was nearly missed twice.** A five-lens adversarial audit
of the staging. Four lenses found no defect; the fifth found both of these. It was then
nearly discarded: its evidence included a 4-subject toy demonstrating that -1 scores like
a swap, and that toy was mistaken for the source of its 31/44 figure, which actually came
from all 660 staged arrays. Its `ld_loss` prediction reproduced to every digit. **When a
verifier disagrees with your own counterfactual, first establish which question each of
you answered** — here both answers were right, and the distance between them was the
finding.

## [A gate that clears its own precondition can never fire, so it is decoration] -> [Key the shared state instead of clearing it, and make the test try to poison it]

**Symptom.** `arm2_attack_a3a4._A3` caches the bi-allelic site axis in an UNKEYED module
global. Arm 1's axis is per target and measured `t_chain` spans only 332,516..332,521, so
many targets share a T -- and where they do, a stale cache is served silently and every
A3 number for that target is computed on another target's site axis. `arm1_attack._a3_for`
was written to close that: clear `_A3`, call `a3_sites`, then assert the returned mask
length against this target's chain.

**The assertion could not fail.** Nothing can be stale if the cache is always empty, so
the length check only ever compared `a3_sites`'s output against the panel it had just
been derived from. `CLAUDE.md`'s rule names this exactly: a gate that cannot fail is
decoration. It was caught by writing the test FIRST as an injection -- poison `_A3` from
another panel, call `_a3_for`, expect a raise -- and watching it report DID NOT RAISE.

**Proven fix.** Key the cache rather than clearing it, on a sha256 of the panel bytes:
  * same key  -> serve it, which is the point of a cache (tested: the same panel twice
                 returns the identical object, so no target recomputes);
  * other key -> clear and recompute, the stale-serve path, caught BY KEY and not by
                 length (tested with two DIFFERENT panels of the SAME length, which is
                 the case a length check passes on);
  * no key    -> RAISE. A non-empty cache this driver did not stamp was populated by some
                 other caller, so whose chain it holds is unknowable (tested).

⚠ **The general test, which is the reusable part.** For any gate, ask what state would
have to be wrong for it to fire, then ask whether the code immediately before it makes
that state impossible. If it does, the gate is asserting its own setup. Write the
injection test before trusting the assertion; "the gate is there" and "the gate can fire"
are different claims and only the second one is worth anything.
⚠ A related trap in the same function: `_a3_for` returns the LIVE module global, not a
copy, so a caller holding it across a `clear()` finds an empty dict. That cost one failing
test too. The driver is safe (one target at a time, nothing clears mid-target) but the
aliasing is now stated at the call site.

## [`pytest ... | tail` returns tail's exit code, so `&&` proceeds and a failing gate is pushed] -> [Redirect to a file, echo `$?`, and never put a gate inside a pipeline]

**Symptom.** A commit went out, and was PUSHED, with `tests/test_docs_contract.py`
failing. The command was
`python -m pytest tests/test_docs_contract.py -q 2>&1 | tail -2 && git add ... && git commit ... && git push`.
A shell pipeline's exit status is the status of the LAST command, so `tail` returning 0
satisfied the `&&` and the whole chain ran on a red gate. The failure text was printed
right above the commit hash and the chain continued anyway.

**Why it kept happening.** This was the THIRD occurrence in one session. Earlier the same
pattern reported `PYTEST EXIT=0` for a run that had actually died on an unrecognised
`--timeout` argument, and again when a results script raised `SystemExit` with a message
that `| tail` swallowed into a 0. Each time the output LOOKED like evidence the step had
passed, because the error text and the success exit sat side by side.

**Proven fix.** Never put a gate inside a pipeline. Redirect and check explicitly:
```
python -m pytest tests/ -q > /tmp/out 2>&1; echo "EXIT=$?"; tail -3 /tmp/out
```
`set -o pipefail` also works inside a script, but the redirect form is what survives being
pasted into a one-off command, which is where this bug lives.

⚠ **The general shape.** A verification step whose exit code can be lost is not a
verification step. If a gate's result reaches you through a pipe, a `grep`, a `tail` or a
`|| true`, you are reading its OUTPUT and not its VERDICT -- and this repository's standing
lesson is that output which looks plausible is exactly what gets believed. The fix for the
plateau claim, for the broken one-liners and for this are all the same instruction: check
the instrument that reports the check.
