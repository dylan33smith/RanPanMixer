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
