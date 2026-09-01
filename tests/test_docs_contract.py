"""The documentation contract.

A documentation system with no test rots. This runs at the end of every session
(see the Wrap-Up Protocol in CLAUDE.md) and fails loudly when the docs, the code
and the disk disagree.

Every failure message says what to do. The standing instruction is:
    the docs and the repo disagree. Fix the docs or fix the code --
    do not quote a number until this passes.

Runnable standalone:  python -m pytest tests/test_docs_contract.py
"""
from __future__ import annotations

import math
import re
from datetime import date
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
DOCS = REPO / "docs"

CLAUDE_MD = REPO / "CLAUDE.md"
PLAN, TERMS, DATA, MEMORY, BUGS = (
    DOCS / n for n in ("plan.md", "terms.md", "data.md", "memory.md", "bugs.md")
)
THE_SIX = [CLAUDE_MD, PLAN, TERMS, DATA, MEMORY, BUGS]

CLAUDE_MD_LINE_BUDGET = 150

EXCLUDED_DIR_NAMES = {".git", ".claude", "archive_docs", "__pycache__",
                      ".pytest_cache", ".ipynb_checkpoints", "node_modules", ".venv"}

FIX = ("the docs and the repo disagree. Fix the docs or fix the code -- "
       "do not quote a number until this passes.")


def _excluded(path: Path) -> bool:
    return any(part in EXCLUDED_DIR_NAMES for part in path.parts)


def _read(p: Path) -> str:
    return p.read_text(encoding="utf-8")


def _scan_files():
    """Every file a stale claim could be hiding in."""
    roots = [REPO / "src", REPO / "tests", REPO / "scripts", REPO / "paper"]
    out = list(THE_SIX)
    for root in roots:
        if not root.exists():
            continue
        out += [f for f in root.rglob("*")
                if f.is_file() and not _excluded(f)
                and f.suffix in {".py", ".md", ".yaml", ".yml", ".json", ".txt"}]
    return out


# ---------------------------------------------------------------------------
# 1. The board is not stale
# ---------------------------------------------------------------------------

def test_plan_last_updated_is_not_older_than_the_newest_memory_entry():
    """Catches "work was done and the board was never reset"."""
    m = re.search(r"\*\*Last updated:\*\*\s*(\d{4}-\d{2}-\d{2})", _read(PLAN))
    assert m, f"plan.md has no '**Last updated:** YYYY-MM-DD' stamp. {FIX}"
    plan_date = date.fromisoformat(m.group(1))

    entry_dates = [
        date.fromisoformat(d)
        for d in re.findall(r"^##\s+(\d{4}-\d{2}-\d{2})\s", _read(MEMORY), re.M)
    ]
    assert entry_dates, f"memory.md has no dated '## YYYY-MM-DD' entries. {FIX}"
    newest = max(entry_dates)

    assert plan_date >= newest, (
        f"plan.md Last updated ({plan_date}) predates the newest memory.md entry "
        f"({newest}). Work was archived but the board was never reset. Run step 6 "
        f"of the Wrap-Up Protocol."
    )


# ---------------------------------------------------------------------------
# 2-3. CLAUDE.md stays a contract, not a results file
# ---------------------------------------------------------------------------

def test_claude_md_is_within_its_line_budget():
    n = len(_read(CLAUDE_MD).splitlines())
    assert n <= CLAUDE_MD_LINE_BUDGET, (
        f"CLAUDE.md is {n} lines, over the {CLAUDE_MD_LINE_BUDGET}-line budget. It "
        f"costs context on every single turn. Cut a section -- do not raise the "
        f"limit. Detail belongs in docs/plan.md (state), docs/terms.md "
        f"(definitions) or docs/data.md (paths)."
    )


_RESULT_WORDS = (r"(utility_retained|attacker_accuracy|tv_empirical|p_succ|"
                 r"accuracy|retained|ceiling|floor|leakage rate)")
_NUMBER = r"\d*\.\d+"


def test_claude_md_contains_no_results():
    """The contract must stay findings-free or it becomes a stale results file."""
    offenders = []
    for i, line in enumerate(_read(CLAUDE_MD).splitlines(), 1):
        low = line.lower()
        for m in re.finditer(_NUMBER, low):
            window = low[max(0, m.start() - 60): m.end() + 60]
            if re.search(_RESULT_WORDS, window):
                offenders.append(f"  CLAUDE.md:{i}: {line.strip()}")
                break
    assert not offenders, (
        "CLAUDE.md contains result-like figures:\n" + "\n".join(offenders) +
        "\nThe contract holds ZERO findings, results and numbers. Move the value to "
        "docs/terms.md (as a definition) or docs/plan.md (as current state)."
    )


# ---------------------------------------------------------------------------
# 4. Every path the docs reference actually exists
# ---------------------------------------------------------------------------

# A row explicitly marked absent is not required to exist -- the doc is being
# honest about the gap, which is the whole point of the status column.
_ABSENT_MARKERS = ("MISSING", "PLANNED", "DEPRECATED", "Removed", "removed",
                   "no longer", "not yet", "none yet", "not built")

# The second alternative needs a >=2-character first segment: "n/a" is prose,
# not a path, and matching it produced this check's first false positive.
_PATHY = re.compile(r"`([A-Za-z0-9_][A-Za-z0-9_./<>*-]*\.[A-Za-z0-9_]+"
                    r"|[a-z_]{2,}/[A-Za-z0-9_./-]*)`")


def _candidate_paths(text: str):
    for line in text.splitlines():
        if any(mark in line for mark in _ABSENT_MARKERS):
            continue
        for raw in _PATHY.findall(line):
            if any(ch in raw for ch in "<>*") or raw.endswith("/"):
                continue
            if "..." in raw or raw.startswith("http"):
                continue
            # A bare basename in prose is a reference, not a path claim.
            if "/" not in raw:
                continue
            if ":" in raw:
                raw = raw.split(":", 1)[0]
            yield line, raw


def test_every_referenced_path_exists():
    missing = []
    for doc in THE_SIX:
        for line, raw in _candidate_paths(_read(doc)):
            p = REPO / raw
            if p.exists() or p.is_symlink():
                continue
            missing.append(f"  {doc.relative_to(REPO)}: `{raw}`  <- {line.strip()[:90]}")
    assert not missing, (
        "Documented paths do not exist:\n" + "\n".join(sorted(set(missing))) +
        f"\n{FIX} If the path is genuinely gone, mark its row MISSING in "
        "docs/data.md rather than deleting the row."
    )


def test_data_md_layout_rows_marked_ok_exist_on_disk():
    """The storage-layout table is a claim about the disk. Check it directly.

    The generic path check above only sees paths with a directory component, so a
    single-segment root would otherwise never be verified.
    """
    wrong = []
    for line in _read(DATA).splitlines():
        m = re.match(r"\|\s*`([A-Za-z0-9_./-]+)`\s*\|.*\|\s*`(OK|PLANNED|MISSING|STALE|DEPRECATED)`\s*\|\s*$", line)
        if not m:
            continue
        rel, state = m.group(1), m.group(2)
        exists = (REPO / rel).exists()
        if state == "OK" and not exists:
            wrong.append(f"  `{rel}` is marked OK but does not exist")
        if state == "PLANNED" and exists:
            wrong.append(f"  `{rel}` is marked PLANNED but DOES exist -- "
                         f"promote the row to OK and describe what is in it")
    assert not wrong, (
        "docs/data.md storage layout disagrees with the disk:\n" + "\n".join(wrong) +
        f"\n{FIX}"
    )


# ---------------------------------------------------------------------------
# 5. Every directory that exists is registered
# ---------------------------------------------------------------------------

def test_every_top_level_directory_is_registered_in_data_md():
    data_txt = _read(DATA)
    unregistered = [
        f"  {c.name}/" for c in sorted(REPO.iterdir())
        if c.is_dir() and c.name not in EXCLUDED_DIR_NAMES
        and not c.name.startswith(".") and c.name not in data_txt
    ]
    assert not unregistered, (
        "Directories exist but are not registered in docs/data.md:\n" +
        "\n".join(unregistered) +
        "\nRun step 5 of the Wrap-Up Protocol: every new output, dataset or "
        "directory gets a docs/data.md row with its state. An unregistered "
        "artifact is unattributable."
    )


# ---------------------------------------------------------------------------
# 6. Terms used in reports resolve to a glossary entry
# ---------------------------------------------------------------------------

# A glossary entry is a ### heading followed by at least one [tag]. That is what
# separates "### tau  [theory]" from a prose subheading like "### Calibration table".
_TERM_HEADING = re.compile(r"^###\s+([A-Za-z_][A-Za-z0-9_]*)\s+\[", re.M)


def _terms_defined() -> set[str]:
    return set(_TERM_HEADING.findall(_read(TERMS)))


def test_glossary_is_not_empty():
    assert _terms_defined(), (
        "docs/terms.md defines no terms. A '### name  [tag]' heading is what makes "
        "an entry visible to this contract."
    )


def test_report_row_and_column_labels_resolve_to_terms():
    """Row labels are the exact terms.md identifier. No prose synonyms."""
    defined = _terms_defined()
    used: set[str] = set()
    for line in _read(PLAN).splitlines():
        if not line.startswith("|"):
            continue
        for cell in line.strip("|").split("|"):
            used.update(re.findall(r"`([a-z][a-z0-9_]{2,})`", cell))
    unresolved = sorted(u for u in used if u not in defined)
    assert not unresolved, (
        "Report labels in docs/plan.md do not resolve to a docs/terms.md entry:\n" +
        "\n".join(f"  {u}" for u in unresolved) +
        "\nRun step 4 of the Wrap-Up Protocol: any new metric or quantity gets a "
        "full terms.md entry. Do not use a prose synonym in a table."
    )


def test_glossary_entries_carry_the_required_fields():
    """Is / Computed by / Status. Without them an entry is a name, not a definition."""
    text = _read(TERMS)
    blocks = re.split(r"^###\s+", text, flags=re.M)[1:]
    incomplete = []
    for block in blocks:
        head = block.split("\n", 1)[0]
        m = re.match(r"([A-Za-z_][A-Za-z0-9_]*)\s+\[", head)
        if not m:
            continue
        name = m.group(1)
        for field in ("Is:", "Computed by:", "Status:"):
            if field not in block:
                incomplete.append(f"  {name} is missing '{field}'")
        if "CHANGES MEANING WITH:" not in block and "Status:" in block \
                and "RETIRED" not in block:
            incomplete.append(
                f"  {name} is missing 'CHANGES MEANING WITH:' -- the field that "
                f"prevents the worst class of error")
    assert not incomplete, (
        "Incomplete docs/terms.md entries:\n" + "\n".join(incomplete))


def test_no_retired_term_is_a_live_endpoint():
    retired_section = _read(TERMS).split("## Retired", 1)
    if len(retired_section) < 2:
        pytest.skip("no Retired section in terms.md yet")
    # Match ASCII hyphens AND en/em dashes. The first version of this check
    # matched only "-" and silently skipped a genuinely retired term.
    retired = set(re.findall(r"`?([a-z][a-z0-9_]{2,})`?\s*[-\u2013\u2014]+\s*RETIRED",
                             retired_section[1]))
    if not retired:
        pytest.skip("nothing retired yet")
    plan = _read(PLAN)
    live = plan.split("## Dropped", 1)[0]
    offenders = [f"  {r}" for r in sorted(retired) if r in live]
    assert not offenders, (
        "Retired terms appear as live work in docs/plan.md:\n" + "\n".join(offenders) +
        "\nA retired name must only appear under Dropped, with its reason.")


# ---------------------------------------------------------------------------
# 7. The theory the project rests on is transcribed correctly
# ---------------------------------------------------------------------------
# This project's primary artifact is a proof, so its highest-severity docs error
# is a mis-transcribed constant. These recompute the published table from its own
# definitions. A silent factor-of-two here understates the leakage the mechanism
# actually permits.

def _calibration_rows():
    rows = []
    for line in _read(TERMS).splitlines():
        m = re.match(r"\|\s*([0-9.]+)\s*\|\s*([0-9.]+)\s*\|\s*([0-9.]+)\s*\|\s*([0-9.]+)\s*\|\s*$", line)
        if m:
            rows.append(tuple(float(g) for g in m.groups()))
    return rows


def test_the_tau_calibration_table_recomputes():
    rows = _calibration_rows()
    assert len(rows) >= 5, (
        "docs/terms.md has no usable tau calibration table (found "
        f"{len(rows)} numeric rows). It is the transcription check for the whole "
        "mechanism; do not delete it.")
    bad = []
    for tau, eta, eps, psucc in rows:
        assert 0.0 <= tau < 1.0, f"tau={tau} is outside [0,1)"
        want_eta = math.atanh(tau)
        want_eps = 2 * math.atanh(tau)
        want_ps = (1 + tau) / 2
        if abs(eta - want_eta) > 5e-7:
            bad.append(f"  tau={tau}: eta_tau documented {eta}, arctanh(tau) is {want_eta:.6f}")
        if abs(eps - want_eps) > 5e-7:
            bad.append(f"  tau={tau}: epsilon_tau documented {eps}, 2*arctanh(tau) is {want_eps:.6f}")
        if abs(psucc - want_ps) > 5e-5:
            bad.append(f"  tau={tau}: p_succ_bound documented {psucc}, (1+tau)/2 is {want_ps:.4f}")
    assert not bad, (
        "The tau calibration table in docs/terms.md does not recompute:\n" +
        "\n".join(bad) + f"\n{FIX}")


def test_epsilon_tau_is_never_defined_as_arctanh_alone():
    """epsilon_tau = 2*arctanh(tau). Dropping the 2 halves the stated budget.

    It is the exact error the source manuscript's own typesetting invites (see
    docs/memory.md, 2026-08-29), which is why it is checked rather than trusted.
    """
    offenders = []
    for f in _scan_files():
        try:
            lines = _read(f).splitlines()
        except (UnicodeDecodeError, OSError):
            continue
        for i, line in enumerate(lines, 1):
            if not re.search(r"eps(ilon)?_tau\s*=\s*(math\.)?arctanh|"
                             r"eps(ilon)?_tau\s*=\s*(math\.)?atanh", line):
                continue
            if re.search(r"\b2\s*\*", line) or re.search(r"NEVER|never|not\b|forbidden|wrong", line):
                continue
            offenders.append(f"  {f.relative_to(REPO)}:{i}: {line.strip()[:110]}")
    assert not offenders, (
        "epsilon_tau defined without its factor of 2:\n" + "\n".join(offenders) +
        "\nepsilon_tau = 2*arctanh(tau) = log((1+tau)/(1-tau)). Half of it is "
        "eta_tau, the tilt magnitude -- a different quantity. This understates "
        "the leakage the mechanism permits.")


# ---------------------------------------------------------------------------
# 8. Every tau-dependent number is stamped with its tau
# ---------------------------------------------------------------------------

# These quantities are meaningless without the privacy parameter they were
# measured at. An unstamped one is not imprecise -- it is uninterpretable.
_NEEDS_TAU = ("utility_retained", "attacker_accuracy", "tv_empirical", "log_z_p")


def test_tau_dependent_numbers_carry_their_tau():
    offenders = []
    for doc in (PLAN, MEMORY, DATA):
        lines = _read(doc).splitlines()
        for i, line in enumerate(lines, 1):
            if not any(t in line for t in _NEEDS_TAU):
                continue
            if not re.search(_NUMBER, line):
                continue
            ctx = " ".join(lines[max(0, i - 3): i + 2])
            if re.search(r"\btau\b|`tau`|eta_tau", ctx):
                continue
            offenders.append(f"  {doc.relative_to(REPO)}:{i}: {line.strip()[:110]}")
    assert not offenders, (
        "tau-dependent figures quoted without their tau:\n" + "\n".join(offenders) +
        "\nSee the reporting format in CLAUDE.md: every number is stamped with its "
        "tau, its denominator and its n.")


# ---------------------------------------------------------------------------
# 9. Corrections do not survive elsewhere (the known failure mode)
# ---------------------------------------------------------------------------

_ANNOTATION = re.compile(
    r"(INCORRECT|CORRECTION|superseded|SUPERSEDES|stale|STALE|RETIRED|DIAGNOSTIC|"
    r"DEMOTED|historical|no longer|was the|formerly|deprecated|DEPRECATED)", re.I)


def _retracted_probes() -> list[str]:
    """Distinctive tokens from retracted lines.

    Two rules inherited from the framework, both learned the hard way:
      * A number needs >=4 decimals to be distinctive.
      * A number the paired [CORRECTION] restates is NOT retracted -- what was
        wrong was the interpretation. Probing it generates pure noise.
    """
    lines = _read(MEMORY).splitlines()
    elsewhere = " ".join(l for l in lines if not l.startswith("[INCORRECT]"))
    probes: set[str] = set()
    for line in lines:
        if not line.startswith("[INCORRECT]"):
            continue
        for num in re.findall(r"\d+\.\d{4,}", line):
            if num not in elsewhere:
                probes.add(num)
        for ident in re.findall(r"[a-z][a-z0-9]*(?:_[a-z0-9]+){2,}", line):
            if ident not in elsewhere:
                probes.add(ident)
    return sorted(probes)


def test_corrections_do_not_survive_elsewhere():
    """Grep the refuted NUMBER, not the topic. Numbers are the high-signal probe."""
    if not any(l.startswith("[INCORRECT]") for l in _read(MEMORY).splitlines()):
        pytest.skip("nothing has been retracted yet -- this check arms itself on "
                    "the first in-place correction")
    probes = _retracted_probes()
    assert probes, (
        "memory.md has [INCORRECT] lines but none carries a distinctive probe (a "
        ">=4-decimal number or a multi-underscore identifier). Retractions that "
        "cannot be grepped cannot be propagated.")

    survivors = []
    for f in _scan_files():
        if f == MEMORY:
            continue
        try:
            lines = _read(f).splitlines()
        except (UnicodeDecodeError, OSError):
            continue
        for i, line in enumerate(lines, 1):
            for probe in probes:
                if probe not in line:
                    continue
                if _ANNOTATION.search(" ".join(lines[max(0, i - 4): i + 3])):
                    continue
                survivors.append(f"  {f.relative_to(REPO)}:{i}: {line.strip()[:110]}")
    assert not survivors, (
        "Retracted values are still asserted as fact outside memory.md:\n" +
        "\n".join(sorted(set(survivors))) +
        "\nA correction was recorded but never propagated. Either update the line, "
        "or annotate it as historical.")


def test_corrections_are_well_formed():
    lines = _read(MEMORY).splitlines()
    bad = []
    for i, line in enumerate(lines):
        if not line.startswith("[INCORRECT]"):
            continue
        nxt = " ".join(lines[i + 1: i + 3])
        if not re.search(r"\[CORRECTION - \d{4}-\d{2}-\d{2}\]", nxt):
            bad.append(f"  memory.md:{i + 1}: {line.strip()[:100]}")
    assert not bad, (
        "[INCORRECT] lines without a dated [CORRECTION - YYYY-MM-DD] directly "
        "below:\n" + "\n".join(bad) +
        "\nThe wrong version stays, but it must be followed by what replaced it.")


# ---------------------------------------------------------------------------
# 10-13. Naming and filesystem hygiene
# ---------------------------------------------------------------------------

def test_no_two_paths_differ_only_in_case():
    seen: dict[str, str] = {}
    collisions = []
    for p in REPO.rglob("*"):
        if _excluded(p):
            continue
        rel = str(p.relative_to(REPO))
        low = rel.lower()
        if low in seen and seen[low] != rel:
            collisions.append(f"  {seen[low]}  vs  {rel}")
        seen[low] = rel
    assert not collisions, (
        "Paths differing only in case:\n" + "\n".join(collisions) +
        "\nRename one. Case-only differences break on case-insensitive filesystems.")


def test_no_brace_shorthand_in_docs():
    """Shorthand is not greppable and a verifier cannot check it."""
    offenders = []
    for doc in THE_SIX:
        for i, line in enumerate(_read(doc).splitlines(), 1):
            for m in re.finditer(r"`[^`]*\{[^`}]*,[^`}]*\}[^`]*`", line):
                offenders.append(f"  {doc.relative_to(REPO)}:{i}: {m.group(0)}")
    assert not offenders, (
        "Brace/glob shorthand in docs:\n" + "\n".join(offenders) +
        "\nWrite the paths out in full (tau_0p10/, tau_0p25/) so they are "
        "greppable and checkable.")


def test_no_self_referential_symlinks():
    """A symlink pointing at its own ancestor makes every recursive walk infinite."""
    offenders = []
    for p in REPO.rglob("*"):
        if _excluded(p) or not p.is_symlink():
            continue
        try:
            target = p.resolve()
        except (OSError, RuntimeError):
            offenders.append(f"  {p.relative_to(REPO)} -> <unresolvable>")
            continue
        parent = p.parent.resolve()
        if target == parent or target in parent.parents:
            offenders.append(f"  {p.relative_to(REPO)} -> {target}")
    assert not offenders, (
        "Self-referential symlinks (infinite recursion in any tree walk):\n" +
        "\n".join(offenders) + "\nRemove the link. See docs/bugs.md.")


def test_the_six_files_exist_and_nothing_else_is_a_working_doc():
    for f in THE_SIX:
        assert f.exists(), f"Missing documentation file: {f.relative_to(REPO)}"
    strays = sorted(
        str(p.relative_to(REPO)) for p in DOCS.rglob("*")
        if p.is_file() and p.resolve() not in {f.resolve() for f in THE_SIX}
    )
    assert not strays, (
        f"Extra files in docs/: {strays}. The system is exactly six files. Content "
        "belongs in one of them, in paper/, or in archive_docs/.")


def test_nothing_was_lost_in_migration():
    """Source material moved out of docs/ must still exist somewhere."""
    for p in (REPO / "paper" / "Private_Genome_Path_Release.pdf",
              REPO / "archive_docs" / "The_Six_File_Lab_Record.pdf"):
        assert p.exists(), (
            f"{p.relative_to(REPO)} is gone. It was moved out of docs/ during the "
            f"2026-08-29 migration; see docs/memory.md. Restore it from git.")
