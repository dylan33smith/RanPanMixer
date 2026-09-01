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

## Mechanism / sampler

No entries yet — no code exists. When the first one lands, the classes most worth
watching are already known from the proposal and are written as constraints in
`CLAUDE.md`: a target-dependent `output_support` or `baseline_model`, an
unnormalized `u_path`, `beta_t` not summing to 1, a `diploid_utility` that sums
instead of averaging, and argmax substituted for `ffbs`. Every one of them
produces a plausible path and no error.
