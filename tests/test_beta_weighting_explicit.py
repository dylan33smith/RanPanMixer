"""GATES on the 2026-10-25 defect: an attack script that silently took the DEFAULT
`beta_v` weighting while the table it fed was computed at another one.

WHAT WENT WRONG. `utility.beta_v` has a `weights=` keyword whose default is
`"inv_support"`. `scripts/arm2_attack.py` called it positionally and never passed
`weights=`, so A1 -- and A3/A4, which import `load_state` from it -- ran the whole Arm 2
attack suite at `inv_support` while `eta_sweep.py`, `exact_kl_stable.py` and
`arm1_release.py` all ran at `uniform`. Nothing crashed and nothing looked wrong: the
tables printed, the controls passed, and the numbers were plausible. Standing constraint
5 forbids mixing two utility functions in one comparison, and the mismatch is not
cosmetic -- chr21 `var_q_u` is 1.3741e-06 at `inv_support` against 1.7500e-06 at
`uniform`, so `eta*` = 1/sqrt(var) is 853 against 756 and one `eta` sits at two
different `acc_ceiling`.

WHY THE GATE IS STATIC. The bug was an OMISSION at a call site, so the gate that catches
it has to look at call sites. Checking the two scripts we already fixed would be
decoration -- it could never fail again. Walking the AST of every caller in the repo
fails for the NEXT script someone writes that forgets the keyword, which is the actual
defect class.
"""

from __future__ import annotations

import ast
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
# Where a `beta_v` call may legitimately live. `src/mechanism/utility.py` holds the
# definition and its own internal uses, so it is not a caller.
SEARCH = ("scripts", "src", "tests")


def _beta_v_calls() -> list[tuple[Path, int, bool]]:
    """Every `beta_v(...)` call in the repo: (file, line, passes `weights=`)."""
    out = []
    for sub in SEARCH:
        for f in sorted((REPO / sub).rglob("*.py")):
            if f.resolve() == Path(__file__).resolve():
                continue          # this file names the keyword only in prose
            try:
                tree = ast.parse(f.read_text(), filename=str(f))
            except SyntaxError as e:                       # pragma: no cover
                pytest.fail(f"{f.relative_to(REPO)} does not parse: {e}")
            for node in ast.walk(tree):
                if not isinstance(node, ast.Call):
                    continue
                fn = node.func
                name = (fn.attr if isinstance(fn, ast.Attribute)
                        else fn.id if isinstance(fn, ast.Name) else None)
                if name != "beta_v":
                    continue
                kw = {k.arg for k in node.keywords if k.arg}
                # `**kwargs` forwarding (arg is None) counts as explicit: the caller
                # above it is a separate call site this walk also visits.
                star = any(k.arg is None for k in node.keywords)
                out.append((f.relative_to(REPO), node.lineno,
                            "weights" in kw or star))
    return out


def test_the_walk_actually_finds_the_call_sites():
    """⚠ A gate whose finder returns nothing PASSES vacuously. Pin the floor."""
    calls = _beta_v_calls()
    assert len(calls) >= 4, (
        f"the AST walk found only {len(calls)} beta_v call(s); it used to find at "
        f"least 4 (eta_sweep, exact_kl_stable x2, arm1_release, arm2_attack). Either "
        f"the walk broke or the callers moved -- either way this gate is not looking "
        f"at the thing it claims to guard. Found: {calls}")


# The scripts that produce the PAPER'S numbers under the `eta_tau` knob. Every one of
# these must state its weighting; this is the surface the defect reached.
PRODUCTION = {
    "scripts/eta_sweep.py",
    "scripts/exact_kl_stable.py",
    "scripts/arm1_release.py",
    "scripts/arm2_attack.py",
    "scripts/arm2_attack_a3a4.py",
}

# ⚠ KNOWN-LEGACY call sites that still omit `weights=`, recorded 2026-10-25 so that a
# NEW one fails this gate. These are `tau`-era scripts and tests: the tests exercise the
# default DELIBERATELY, and the scripts are slated for removal in the paper-branch
# cleanup (docs/plan.md). This list may SHRINK, never grow. It is not an excuse -- it is
# the honest boundary of what has been fixed, so the gate can still fail for new code.
LEGACY = {
    "scripts/arm1_exact_endpoint.py", "scripts/arm1_headtohead.py",
    "scripts/arm1_symmetric_compare.py", "scripts/beta_concentration.py",
    "scripts/churn_across_tau.py", "scripts/compare_models_chr21.py",
    "scripts/exact_endpoint.py", "scripts/exact_fidelity.py", "scripts/exact_kl.py",
    "scripts/exact_kl_genome.py", "scripts/exact_utility.py",
    "scripts/exact_variance.py", "scripts/fidelity_frontier.py",
    "scripts/mismatch_map_chr21.py", "scripts/results_table.py",
    "scripts/run_mechanism_chr21.py", "scripts/toy_genome_demo.py",
    "tests/test_diploid_gate.py", "tests/test_exact_kl.py", "tests/test_exact_out.py",
    "tests/test_genomewide_budget.py", "tests/test_mechanism.py",
    "tests/test_toy_enumeration.py",
}


def test_every_production_beta_v_caller_passes_weights_explicitly():
    """THE GATE. The defect was an omitted keyword in a script that fed a table."""
    bad = [(str(f), ln) for f, ln, ok in _beta_v_calls()
           if not ok and str(f).replace("\\", "/") in PRODUCTION]
    assert not bad, (
        "these PRODUCTION beta_v call sites do not pass `weights=` and so silently "
        f"take the dataclass default: {bad}. Pass it explicitly. Two weightings are "
        "DIFFERENT utility functions (constraint 5) and sit at different "
        "`acc_ceiling` for one `eta` (constraint 7), so a default here is a silent "
        "wrong answer, not a convenience.")


def test_no_new_file_starts_omitting_weights():
    """The defect class, not the two instances: a NEW offender anywhere fails here.

    ⚠ This is what makes the gate more than a fixed-two-scripts check. `LEGACY` is a
    frozen census; anything outside it must pass `weights=`.
    """
    offenders = {str(f).replace("\\", "/") for f, _, ok in _beta_v_calls() if not ok}
    new = offenders - LEGACY - PRODUCTION
    assert not new, (
        f"new file(s) omit `weights=` on a beta_v call: {sorted(new)}. Either pass it "
        f"or, if this is genuinely legacy, say so by adding it to LEGACY with a "
        f"reason -- do not let it in silently.")


def test_production_scripts_still_exist_and_are_scanned():
    """⚠ If a PRODUCTION path were renamed, the gate above would pass vacuously."""
    scanned = {str(f).replace("\\", "/") for f, _, _ in _beta_v_calls()}
    # arm2_attack_a3a4 reaches `beta_v` only through A1's `load_state`, so it has no
    # call of its own; it is covered by the CLI test instead.
    expect = PRODUCTION - {"scripts/arm2_attack_a3a4.py"}
    missing = {p for p in expect if not (REPO / p).exists()}
    assert not missing, f"PRODUCTION paths no longer exist: {sorted(missing)}"
    unscanned = expect - scanned
    assert not unscanned, (
        f"these PRODUCTION scripts exist but the walk found no beta_v call in them: "
        f"{sorted(unscanned)}. The call moved or the walk broke; the gate is blind.")


def test_the_gate_can_fail():
    """⚠ Inject the exact omission and assert the checker reports it.

    Without this, `test_every_beta_v_caller_passes_weights_explicitly` passing tells us
    nothing about whether it CAN fail.
    """
    src = "import utility\nu = utility.beta_v(sup, in_chain=ic)\n"
    tree = ast.parse(src)
    found = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and getattr(node.func, "attr", None) == "beta_v":
            found.append({k.arg for k in node.keywords if k.arg})
    assert found == [{"in_chain"}], f"the detector did not see the call: {found}"
    assert "weights" not in found[0], "the injected omission was not detected"


@pytest.mark.parametrize("script", ["arm2_attack.py", "arm2_attack_a3a4.py"])
def test_attack_cli_defaults_to_the_released_weighting(script):
    """`--weights` must EXIST and default to `uniform` on both attack entry points.

    A1 and A3/A4 share `load_state`, so a flag on one and not the other would leave
    half the suite on the stale default -- which is how the defect reached a table.
    """
    r = subprocess.run([sys.executable, str(REPO / "scripts" / script), "--help"],
                       capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, f"{script} --help failed: {r.stderr[-2000:]}"
    assert "--weights" in r.stdout, (
        f"{script} has no --weights flag, so its weighting cannot be stated at the "
        f"command line and will be whatever the dataclass default happens to be")

    # THE DEFAULT, read from the parser rather than from prose. `--weights BOGUS` must
    # be refused by `choices`, which also proves the flag is wired to a real argument
    # and not just mentioned in a help string.
    bad = subprocess.run(
        [sys.executable, str(REPO / "scripts" / script), "--weights", "BOGUS"],
        capture_output=True, text=True, timeout=300)
    assert bad.returncode != 0 and "invalid choice" in bad.stderr, (
        f"{script} accepted --weights BOGUS, so the flag is not constrained by "
        f"choices and a typo would silently select no weighting: {bad.stderr[-500:]}")

    src = (REPO / "scripts" / script).read_text()
    assert 'choices=("uniform", "inv_support"), default="uniform"' in src, (
        f"{script} does not declare --weights with default 'uniform'; the released "
        f"tool's weighting must be what you get when the flag is omitted")


def test_the_two_weightings_are_actually_different():
    """⚠ If `uniform` and `inv_support` gave the same `beta`, the defect would be
    harmless and this whole gate pointless. They do not -- and the DIFFERENCE is what
    makes the mismatch a wrong answer rather than a naming quibble."""
    sys.path.insert(0, str(REPO / "src"))
    import numpy as np
    from mechanism import utility

    rng = np.random.default_rng(20261025)
    n = 5000
    # a support vector with real spread; 1/support(v) is flat only if support is
    support = rng.integers(1, 89, size=n).astype(float)
    ic = np.ones(n, dtype=bool)

    bu = utility.beta_v(support, in_chain=ic, weights="uniform").beta
    bi = utility.beta_v(support, in_chain=ic, weights="inv_support").beta

    assert np.isclose(bu.sum(), 1.0) and np.isclose(bi.sum(), 1.0)
    assert not np.allclose(bu, bi), (
        "uniform and inv_support produced the same beta; either the weighting branch "
        "is dead or this test's support vector has no spread")
    # `var_q_u` scales with sum(beta^2) at fixed structure, and that is the quantity
    # that moved `eta*` from 853 to 756 on chr21.
    su, si = float((bu ** 2).sum()), float((bi ** 2).sum())
    assert not np.isclose(su, si, rtol=1e-6), (
        f"sum(beta^2) identical ({su} vs {si}); the weightings would then sit at the "
        f"same acc_ceiling and the mismatch could not have changed a table")
