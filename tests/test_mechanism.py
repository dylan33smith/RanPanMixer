"""Unit and structural tests for the mechanism.

`tests/test_toy_enumeration.py` is the gate that can falsify Theorem 1. This file
covers the things that gate cannot see: the rules the code must obey structurally
(no target in `baseline_model`, no MAP decoding anywhere, no allele-frequency
branch), the guards that should fire on bad input, and the claims made in prose
that would otherwise go unchecked.

Runs two ways, matching `tests/test_preprocess.py`:

    python3 tests/test_mechanism.py
    python3 -m pytest -q tests/test_mechanism.py
"""

from __future__ import annotations

import ast
import math
import sys
import tracemalloc
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "src"
sys.path.insert(0, str(SRC))

from mechanism import calibration as C   # noqa: E402
from mechanism import exact as E         # noqa: E402
from mechanism import ffbs as F          # noqa: E402
from mechanism import hmm as H           # noqa: E402
from mechanism import toy as TOY         # noqa: E402
from mechanism import utility as U       # noqa: E402

MECH = SRC / "mechanism"


def _toy_bundle():
    t = TOY.toy_cohort()
    return t, U.beta_v(t.support), H.from_arrays(t.haplotypes, t.genetic_pos,
                                                 t.positions, cohort_id=t.name)


# ---------------------------------------------------------------------------
# Structural constraints. These read the SOURCE, so the rule breaks at CI time
# rather than silently at release time -- the same device `cohort/` uses, scoped
# here to one module because `utility.py` and `ffbs.py` must see the target.
# ---------------------------------------------------------------------------

def test_the_baseline_model_never_imports_or_names_a_target():
    """Standing constraint 1, enforced by construction rather than discipline.

    `hmm.py` builds `R_D` and `Y_D`, which must be functions of `G` and `D`
    alone. A target-derived quantity reaching the constructor is the failure
    `CLAUDE.md` calls the easiest way to buy utility and silently lose the
    guarantee, so it is checked from the AST and not by grep.
    """
    tree = ast.parse((MECH / "hmm.py").read_text())
    offenders = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            offenders += [f"import {n.name}" for n in node.names
                          if n.name.split(".")[0] in ("target", "utility", "ffbs")]
        elif isinstance(node, ast.ImportFrom):
            mod = (node.module or "").lstrip(".").split(".")[0]
            if mod in ("target", "utility", "ffbs"):
                offenders.append(f"from {node.module} import ...")
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            args = [a.arg for a in node.args.args + node.args.kwonlyargs]
            bad = [a for a in args if "target" in a.lower()]
            if bad:
                offenders.append(f"{node.name}({', '.join(bad)})")
    assert not offenders, (
        "hmm.py has acquired a target dependency:\n  " + "\n  ".join(offenders) +
        "\nbaseline_model is a function of G and D only.")


def test_the_enumeration_reference_is_independent_of_the_implementation():
    """A check that shares code with what it checks cannot fail.

    `exact.py` is the falsification target, so it must not import `hmm.py`,
    `ffbs.py` or `utility.py`. It restates the published transition formula and
    the utility instead. This project has already documented what the alternative
    looks like: PanMixer's own shape assertion compares two expressions that are
    equal by construction, so it can never fire.
    """
    for mod in ("exact.py", "toy.py"):
        tree = ast.parse((MECH / mod).read_text())
        for node in ast.walk(tree):
            names = []
            if isinstance(node, ast.Import):
                names = [n.name for n in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [(node.module or "").lstrip(".")]
            for n in names:
                assert n.split(".")[0] not in ("hmm", "ffbs", "utility", "mechanism"), (
                    f"{mod} imports {n}; the reference must be independent")


def test_the_toy_and_the_reference_read_nothing_from_disk():
    """The falsification gate must not depend on `data/`, or on anything at all.

    `docs/memory.md`, `docs/plan.md` and `docs/data.md` all state that the toy is
    deliberately independent of every artifact -- no graph, no VCF, no cohort
    output, no PanMixer path. That is load-bearing: a defect in preprocessing
    must not be able to mask a defect in the mechanism, and the gate has to keep
    working on a checkout with no `data/` symlink. So the claim gets a test
    rather than three restatements in prose.
    """
    banned = {"load", "loadtxt", "open", "read_text", "read_bytes", "save",
              "savetxt", "write_text", "genfromtxt", "memmap"}
    offenders = []
    for mod in ("toy.py", "exact.py"):
        tree = ast.parse((MECH / mod).read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                fn = node.func
                name = fn.attr if isinstance(fn, ast.Attribute) else getattr(fn, "id", None)
                if name in banned:
                    offenders.append(f"{mod}:{node.lineno}: {name}()")
    assert not offenders, (
        "the toy or the reference touches the filesystem:\n  " +
        "\n  ".join(offenders))


def test_the_reference_and_the_model_state_the_constants_separately():
    """Both must carry Ne and r, and they must currently agree.

    Stating them twice is the point -- it is what lets the two implementations
    disagree. This test is the tripwire for someone "removing the duplication",
    which would quietly turn the cross-check into a tautology.
    """
    assert (H.NE, H.R) == (10_000.0, 1.26)
    assert (E.NE, E.R) == (10_000.0, 1.26)
    for mod in ("hmm.py", "exact.py"):
        src = (MECH / mod).read_text()
        assert "NE = 10_000.0" in src and "R = 1.26" in src, (
            f"{mod} no longer states the published constants itself")


def test_no_map_decoding_anywhere_in_the_mechanism():
    """`viterbi_release` was RETIRED before use; it must not return by accident.

    MAP decoding returns the single highest-probability state sequence. It is not
    a draw from `Q^tau_p` and carries no guarantee, while producing output that
    looks entirely plausible -- the failure is completely silent, and in an HMM
    codebase it is the natural thing to reach for. Standing constraint 3.

    ⚠ Parsed from the AST, not grepped. The first version searched the source
    TEXT and fired on the docstrings explaining why the code does not do this --
    the identical mistake `test_preprocess.py` records about the shared-parser
    check: a test that reads prose rather than behaviour is a test that will keep
    lying. Names and attributes are real uses; docstrings are string constants
    and comments are not in the tree at all.
    """
    banned = {"argmax", "argmin", "amax", "viterbi", "argpartition"}
    offenders = []
    for py in sorted(MECH.glob("*.py")):
        tree = ast.parse(py.read_text())
        for node in ast.walk(tree):
            name = None
            if isinstance(node, ast.Attribute):
                name = node.attr
            elif isinstance(node, ast.Name):
                name = node.id
            if name and name.lower() in banned:
                offenders.append(f"{py.name}:{node.lineno}: {name}")
    assert not offenders, (
        "MAP/argmax decoding in the mechanism:\n  " + "\n  ".join(offenders) +
        "\nSample; never argmax (standing constraint 3).")


def test_no_allele_frequency_branch_exists():
    """Exit gate 2 of the sampler spec, as a source check.

    Model C has one code path for every chain position: no anchors, no singleton
    special case, and no allele-frequency fallback. PanMixer draws 89.6% of its
    blocks from independent per-variant frequencies instead of from the cohort
    chain, which is the behaviour this gate exists to keep out.
    """
    offenders = []
    for py in sorted(MECH.glob("*.py")):
        for line_no, line in enumerate(py.read_text().splitlines(), 1):
            code = line.split("#", 1)[0].lower()
            if "allele_freq" in code or "allele_frequenc" in code:
                offenders.append(f"{py.name}:{line_no}")
    assert not offenders, (
        "an allele-frequency branch appeared in the mechanism: " + ", ".join(offenders))


# ---------------------------------------------------------------------------
# calibration
# ---------------------------------------------------------------------------

def test_calibration_matches_the_table_in_the_glossary():
    """The code and `docs/terms.md`'s calibration table must not drift apart."""
    table = {0.00: (0.000000, 0.000000, 0.5000), 0.01: (0.010000, 0.020001, 0.5050),
             0.05: (0.050042, 0.100083, 0.5250), 0.10: (0.100335, 0.200671, 0.5500),
             0.25: (0.255413, 0.510826, 0.6250), 0.50: (0.549306, 1.098612, 0.7500),
             0.90: (1.472219, 2.944439, 0.9500)}
    for tau, (eta, eps, ps) in table.items():
        assert abs(C.eta_tau(tau) - eta) < 5e-7, tau
        assert abs(C.epsilon_tau(tau) - eps) < 5e-7, tau
        assert abs(C.p_succ_bound(tau) - ps) < 5e-7, tau


def test_epsilon_tau_is_exactly_twice_eta_tau():
    """⚠ Never `arctanh(tau)` alone -- that halves the stated budget."""
    for tau in (0.0, 0.01, 0.1, 0.25, 0.5, 0.9, 0.999):
        assert abs(C.epsilon_tau(tau) - 2 * C.eta_tau(tau)) < 1e-15


def test_an_unnormalised_utility_changes_the_calibration():
    """Standing constraint 2: range `delta_u` means `arctanh(tau)/delta_u`.

    An unnormalized utility does not weaken the bound -- it removes it. So
    `delta_u` is an argument here rather than an assumption made elsewhere.
    """
    assert abs(C.eta_tau(0.5, delta_u=2.0) - C.eta_tau(0.5) / 2) < 1e-15
    for bad in (-0.1, 1.0, 1.5):
        try:
            C.eta_tau(bad)
        except ValueError:
            continue
        raise AssertionError(f"tau={bad} was accepted")


# ---------------------------------------------------------------------------
# hmm -- baseline_model
# ---------------------------------------------------------------------------

def test_transition_constants_reproduce_the_published_closed_form():
    """`p_stay = e^(-d/n) + (1 - e^(-d/n))/n` with `d = Delta_x * Ne * r`.

    Computed here from the paper's formula directly and compared against the
    rank-one decomposition the sampler uses. ⚠ NOT PanMixer's constants: its
    `EFFECTIVE_N` is 1/10,000 and its distance is a difference of VCF row
    indices, giving a total switch mass of 5.8e-12 where this gives 0.13-0.99.
    """
    cm = np.array([0.0, 0.00002, 0.00012, 0.02])
    for K in (4, 44, 86, 88):
        s = H.switch_mass_from_cm(cm, K)
        for t, dx in enumerate(np.diff(cm)):
            d = dx * 10_000.0 * 1.26
            want = math.exp(-d / K) + (1 - math.exp(-d / K)) / K
            assert abs((1 - s[t] + s[t] / K) - want) < 1e-15, (K, t)


def test_switch_mass_is_the_total_not_the_per_state_value():
    """`stay_switch` warns that confusing them rescales every transition by K.

    `s_t` is the coefficient of `J` in `A = (1-s) I + (s/K) J`, so the chance of
    actually LEAVING a state is `s (K-1)/K` and the per-state switch is `s/K`.
    Pinned by a test because the two differ by a factor of K and nothing in the
    output makes the error visible.
    """
    t, util, m = _toy_bundle()
    for i, s in enumerate(m.switch_mass):
        A = E.dense_transition(float(np.diff(t.genetic_pos)[i]), t.K)
        assert abs(A[0, 1] - s / t.K) < 1e-15, "off-diagonal is not s/K"
        assert abs(A[0, 0] - m.p_stay[i]) < 1e-15, "diagonal is not p_stay"


def test_a_zero_distance_step_forbids_switching_and_a_large_one_goes_uniform():
    """The two ends of the transition model, as data rather than as prose."""
    t, util, m = _toy_bundle()
    flat = int(np.flatnonzero(np.diff(t.genetic_pos) == 0)[0])
    assert m.switch_mass[flat] == 0.0 and m.p_stay[flat] == 1.0
    wide = int(np.flatnonzero(np.diff(t.genetic_pos) > 0.01)[0])
    assert abs(m.p_stay[wide] - 1.0 / t.K) < 1e-9, (
        "a far-apart pair should give a near-uniform transition")


def test_a_clamped_or_missing_cm_is_refused_loudly():
    """`np.interp`'s default clamp is the defect the coordinate cut exists for.

    Interpolating with `left=nan, right=nan` and asserting here is what turns a
    silent 34,371-position failure -- one donor's real haplotype emitted verbatim
    across 7.25 Mb, a privacy failure -- into a crash.
    """
    for bad in (np.array([0.0, np.nan, 0.1]), np.array([0.0, np.inf, 0.1])):
        try:
            H.switch_mass_from_cm(bad, 4)
        except ValueError as e:
            assert "interpolated cM" in str(e)
            continue
        raise AssertionError(f"{bad} was accepted")
    try:
        H.switch_mass_from_cm(np.array([0.0, 0.2, 0.1]), 4)
    except ValueError as e:
        assert "monotonic" in str(e)
    else:
        raise AssertionError("a decreasing map was accepted")


def test_a_position_with_no_called_donor_is_refused_not_patched():
    """`support(v) = 0` means RENORMALISE has an empty pool: unusable, not awkward.

    It is excluded from the chain upstream by `genmap.chain_mask`, so reaching
    here means the mask and the matrix have drifted apart. Falling back would
    hide a position that cannot produce output.
    """
    hap = np.array([[0, -1], [1, -1]], dtype=np.int16)
    try:
        H.from_arrays(hap, np.array([0.0, 0.001]), np.array([1, 2]))
    except ValueError as e:
        assert "NO called donor" in str(e)
    else:
        raise AssertionError("an all-missing position was accepted")


def test_mismatched_artifact_lengths_are_refused():
    t, util, m = _toy_bundle()
    try:
        H.from_arrays(t.haplotypes, t.genetic_pos[:-1], t.positions)
    except ValueError as e:
        assert "different VCFs" in str(e)
    else:
        raise AssertionError("a short genetic_pos was accepted")


def test_the_chain_span_mask_selects_rows_without_renumbering_the_axis():
    """The cut is stored as a mask over the full axis, never applied as a filter.

    It is provisional -- Dylan is settling it with his PI -- so moving it must
    rewrite one small array and leave every other artifact's row numbering, and
    the site-axis digest, untouched. Filtering rows instead would renumber the
    whole tree.
    """
    t, util, m = _toy_bundle()
    mask = np.ones(t.T, dtype=bool)
    mask[0] = mask[-1] = False
    cut = H.from_arrays(t.haplotypes, t.genetic_pos, t.positions, in_chain=mask)
    assert cut.T == t.T - 2 and cut.K == t.K
    assert cut.positions[0] == t.positions[1]
    assert np.allclose(cut.switch_mass, H.switch_mass_from_cm(t.genetic_pos[mask], t.K))


# ---------------------------------------------------------------------------
# utility
# ---------------------------------------------------------------------------

def test_beta_sums_to_one_and_is_inverse_support_weighted():
    t, util, m = _toy_bundle()
    assert abs(util.beta.sum() - 1.0) < 1e-15
    w = 1.0 / t.support
    assert np.allclose(util.beta, w / w.sum())
    assert util.scope == "chain" and util.n_scope == util.n_chain == t.T


def test_beta_scope_chain_readable_changes_the_ceiling_not_the_bound():
    """Deferred decision R1, as a test rather than a claim.

    Normalising over `in_chain & readable` moves the achievable ceiling from
    below 1 to 1 and invents nothing -- but the release at the excluded positions
    is a pure prior draw either way, and `sum beta = 1` holds under both, so the
    calibration is untouched. What changes is the DENOMINATOR a utility number
    is quoted against, which is why the scope is stamped on every release.

    ⚠ It is legitimate only because `readable` is a property of the target PANEL.
    For a genuinely novel genome it would differ per target and normalising over
    it would be target-dependent -- see R4 in `docs/plan.md`.
    """
    t, _, _ = _toy_bundle()
    readable = np.ones(t.T, dtype=bool)
    readable[2] = False
    narrow = U.beta_v(t.support, readable=readable, scope="chain_readable")
    assert abs(narrow.beta.sum() - 1.0) < 1e-15
    assert narrow.beta[2] == 0.0 and narrow.n_scope == t.T - 1

    target = t.haplotypes[0].copy()
    target[2] = -1                       # no call exactly where it is unreadable
    wide = U.beta_v(t.support)
    assert U.u_ceiling(target, wide.beta) < 1.0 - 1e-9
    assert abs(U.u_ceiling(target, narrow.beta) - 1.0) < 1e-12

    try:
        U.beta_v(t.support, scope="chain_readable")
    except ValueError as e:
        assert "readable mask" in str(e)
    else:
        raise AssertionError("chain_readable was accepted with no readable mask")


def test_zero_support_is_refused_rather_than_weighted():
    """⚠ Why exclusion and not `1/(support + 1)`, which was the obvious patch.

    That patch keeps an unusable position in the chain and hands it a finite
    `beta_v` share it can never earn -- re-creating the inert-budget problem one
    position at a time, and hiding a position that cannot produce output.
    """
    try:
        U.support_weights(np.array([88, 0, 88]))
    except ValueError as e:
        assert "support 0" in str(e)
    else:
        raise AssertionError("support 0 was weighted")


def test_phi_v_is_an_indicator_and_a_missing_target_is_not_a_wildcard():
    """A target with no call scores ZERO for every donor, never 1.

    A wildcard match would make every donor equally attractive there and reward
    agreement with an unknown value. It is also not `absent_policy =
    "reference"`, which fills allele 0 and would reward donors for matching a
    fabricated call -- absence from a callset means the caller reported no
    variant there, not that this person is reference.
    """
    donors = np.array([0, 1, 2, -1])
    assert list(U.phi_v(donors, 1)) == [0.0, 1.0, 0.0, 0.0]
    assert list(U.phi_v(donors, -1)) == [0.0, 0.0, 0.0, 0.0]
    assert set(np.unique(U.phi_v(donors, 0)).tolist()) <= {0.0, 1.0}


def test_u_path_is_asserted_in_range_not_documented_in_range():
    t, util, m = _toy_bundle()
    assert abs(U.u_of_output(t.haplotypes[0], t.haplotypes[0], util.beta)[0] - 1.0) < 1e-12
    assert U.u_of_output(t.haplotypes[1], t.haplotypes[0], util.beta)[0] == 0.0
    try:
        U.u_of_output(t.haplotypes[0], t.haplotypes[0], util.beta * 3.0)
    except AssertionError as e:
        assert "left [0,1]" in str(e)
    else:
        raise AssertionError("an out-of-range utility was returned rather than raised")


def test_diploid_utility_is_the_mean_and_a_sum_is_caught():
    """⚠ A SUM leaves [0,1] and doubles the effective tilt; `tau` stops covering it."""
    assert abs(U.diploid_utility(1.0, 0.0) - 0.5) < 1e-15
    assert abs(U.diploid_utility(1.0, 1.0) - 1.0) < 1e-15
    try:
        U._assert_bounded(np.array([1.0 + 1.0]))      # what a SUM would produce
    except AssertionError:
        pass
    else:
        raise AssertionError("a summed diploid utility passed the range check")


def test_tilt_potentials_carry_both_the_mask_and_the_tilt():
    """One array, because applying one and forgetting the other must not be expressible."""
    t, util, m = _toy_bundle()
    eta = C.eta_tau(0.5)
    psi = U.tilt_potentials(m.allowed, t.haplotypes, t.haplotypes[0], util.beta, eta)
    assert psi.shape == (t.T, t.K) and psi.dtype == np.float64
    assert (psi.T[~m.allowed] == 0.0).all(), "an unavailable donor has a nonzero potential"
    match = (t.haplotypes == t.haplotypes[0][None, :]) & m.allowed
    assert np.allclose(psi.T[match], np.exp(eta * np.broadcast_to(
        util.beta[None, :], match.shape)[match]))

    flat = U.tilt_potentials(m.allowed, t.haplotypes, t.haplotypes[0], util.beta, 0.0)
    assert set(np.unique(flat).tolist()) == {0.0, 1.0}, (
        "at eta = 0 every potential must be the bare allowed-state mask")


# ---------------------------------------------------------------------------
# ffbs
# ---------------------------------------------------------------------------

def test_the_forward_pass_never_materialises_a_transition_matrix():
    """A dense `A` per position would be about 21 TB on chr21. Measured, not asserted.

    Peak allocation is compared against what a single K x K matrix would cost, at
    a K where the difference is three orders of magnitude. This is the one claim
    in `ffbs.py` that a correctness test cannot reach: a materialised matrix
    would give identical numbers and simply not run at scale.
    """
    K, T = 3_000, 30
    rng = np.random.default_rng(7)
    hap = rng.integers(0, 2, size=(K, T)).astype(np.int16)
    m = H.from_arrays(hap, np.linspace(0, 0.001, T), np.arange(1, T + 1) * 100)
    psi = np.ascontiguousarray(m.allowed.T.astype(np.float64))

    tracemalloc.start()
    F.forward(m, psi)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    one_matrix = K * K * 8
    alphas = T * K * 8
    assert peak < one_matrix / 10, (
        f"peak allocation {peak/1e6:.1f} MB against {one_matrix/1e6:.1f} MB for a "
        f"single {K}x{K} matrix and {alphas/1e6:.1f} MB of alphas -- A is being "
        f"materialised")


def test_every_alpha_column_is_stored_in_float64():
    """PanMixer keeps only the running vector; backward sampling needs all of them.

    float64 is load-bearing under model C, not a preference: the per-position
    tilt at `tau` = 0.5 is 1.61e-06 against float32's resolution of 1.19e-07, so
    only 13.5x the noise floor. Rounding would consume a meaningful fraction of
    the target's entire influence on its own release.
    """
    t, util, m = _toy_bundle()
    psi = U.tilt_potentials(m.allowed, t.haplotypes, t.haplotypes[0], util.beta,
                            C.eta_tau(0.5))
    fwd = F.forward(m, psi)
    assert fwd.alpha.shape == (t.T, t.K) and fwd.alpha.dtype == np.float64
    assert np.allclose(fwd.alpha.sum(axis=1), 1.0)
    try:
        F.forward(m, psi.astype(np.float32))
    except TypeError as e:
        assert "float64" in str(e)
    else:
        raise AssertionError("a float32 potential array was accepted")

    tilt = math.exp(C.eta_tau(0.5) * (1.0 / 305_886))
    assert (tilt - 1.0) / np.finfo(np.float32).eps < 20, (
        "the float32 headroom claim in ffbs.py is no longer the right order")


def test_log_z_p_needs_the_baseline_subtraction_and_says_so_when_it_is_missing():
    """The trap in point 4 of `ffbs.py`: the conditioning constant is not the tilt.

    RENORMALISE enters as a 0/1 potential, so `sum_t log c_t` carries `R_D`'s own
    normaliser as well as the tilt. Reporting that sum as `log_z_p` puts it far
    outside `[0, eta_tau]`, and the tempting fix -- widening the bracket -- would
    discard the cheapest correctness check in the project.
    """
    t, util, m = _toy_bundle()
    eta = C.eta_tau(0.5)
    psi = U.tilt_potentials(m.allowed, t.haplotypes, t.haplotypes[0], util.beta, eta)
    fwd = F.forward(m, psi)
    base = F.baseline_log_norm(m)
    assert base < 0, "the conditioning constant should be negative here"
    assert 0.0 <= F.log_z_p(fwd, base, eta) <= eta
    try:
        F.log_z_p(fwd, 0.0, eta)          # the forgotten subtraction
    except AssertionError as e:
        assert "outside [0, eta_tau]" in str(e)
    else:
        raise AssertionError("log_z_p without the baseline subtraction passed")


def test_the_baseline_normaliser_is_target_independent():
    """It is cached across every target and every `tau`, so this has to hold."""
    t, util, m = _toy_bundle()
    base = F.baseline_log_norm(m)
    for _ in range(3):
        assert F.baseline_log_norm(m) == base
    # and it equals the tau = 0 pass through the target-facing code path
    for target in TOY.toy_inputs(t, kind="called")[::23]:
        psi = U.tilt_potentials(m.allowed, t.haplotypes, target, util.beta, 0.0)
        assert abs(F.forward(m, psi).log_norm - base) < 1e-15


def test_a_release_is_reproducible_from_its_recorded_seed():
    """`docs/data.md`'s reproducibility contract, checked rather than promised."""
    t, util, m = _toy_bundle()
    a = F.release(m, t.haplotypes, t.haplotypes[0], util, 0.25, seed=99,
                  target_id="toyD1", graph_id=t.name)
    b = F.release(m, t.haplotypes, t.haplotypes[0], util, 0.25, seed=99,
                  target_id="toyD1", graph_id=t.name)
    assert np.array_equal(a.state_path, b.state_path)
    assert a.log_z_p == b.log_z_p and a.u_achieved == b.u_achieved
    c = F.release(m, t.haplotypes, t.haplotypes[0], util, 0.25, seed=100,
                  target_id="toyD1", graph_id=t.name)
    assert not np.array_equal(a.state_path, c.state_path), "the seed does nothing"


def test_a_release_record_is_self_describing():
    """A release that cannot answer "at what tau, under what utility" is not quotable."""
    t, util, m = _toy_bundle()
    r = F.release(m, t.haplotypes, t.haplotypes[0], util, 0.25, seed=1,
                  target_id="toyD1", graph_id=t.name)
    rec = r.record()
    for field in ("release_id", "target_id", "tau", "eta_tau", "utility_id",
                  "graph_id", "loo_id", "seed", "log_z_p", "released_at", "ploidy"):
        assert field in rec, field
    assert rec["loo_id"] is None, "a genuine external-target release holds out nobody"
    assert abs(rec["eta_tau"] - C.eta_tau(rec["tau"])) < 1e-15
    assert rec["utility_id"] == "indicator/w=inv_support/beta=chain"
    assert "tau0.25" in rec["release_id"]


def test_a_diploid_release_is_one_release_at_half_the_exponent_each():
    """Two haplotypes, ONE private input, one `tau`, one `log Z_p` in `[0, eta_tau]`.

    `exp(eta (u1+u2)/2) = exp((eta/2) u1) exp((eta/2) u2)`, so the chains are
    independent at half the exponent and `log Z_p` is the sum of the two -- which
    lands in `[0, eta/2] + [0, eta/2] = [0, eta]`. The bracket therefore holds for
    the diploid release as a whole, at the real `tau`.
    """
    t, util, m = _toy_bundle()
    tau = 0.5
    diploid_target = np.stack([t.haplotypes[0], t.haplotypes[1]], axis=1)
    r = F.release(m, t.haplotypes, diploid_target, util, tau, seed=5,
                  target_id="toyDip", graph_id=t.name)
    assert r.ploidy == 2
    assert r.state_path.shape == (t.T, 2) and r.released_path.shape == (t.T, 2)
    assert r.tau == tau and abs(r.eta_tau - C.eta_tau(tau)) < 1e-15
    assert 0.0 <= r.log_z_p <= C.eta_tau(tau) + 1e-9
    assert 0.0 <= r.u_achieved <= 1.0 + 1e-12
    assert not np.array_equal(r.state_path[:, 0], r.state_path[:, 1]), (
        "both haplotypes drew the identical path; the seeds are not distinct")

    # the halving is the guarantee: at the FULL exponent per haplotype the sum
    # can exceed the bracket, which is the SUM aggregator in disguise
    half = [F.log_z_p(F.forward(m, U.tilt_potentials(
        m.allowed, t.haplotypes, diploid_target[:, h], util.beta, C.eta_tau(tau))),
        F.baseline_log_norm(m), C.eta_tau(tau)) for h in (0, 1)]
    assert sum(half) > C.eta_tau(tau), (
        "running each haplotype at the full eta no longer overshoots, so this "
        "test can no longer tell the mean from the sum")


def test_a_dead_end_chain_fails_loudly_rather_than_returning_a_path():
    """The case that rules out reading RENORMALISE as a renormalised transition row.

    Across a `Delta_x = 0` step the transition is the identity, so a donor
    available at `t-1` and unavailable at `t` has nowhere to go. Under local row
    renormalisation that row sums to zero and there is no chain to normalise;
    here the potentials simply annihilate the mass and the forward pass says so.
    chr21 has 618 interior zero-cM intervals against 4.26% mean missingness, so
    this is reachable on the real data.
    """
    hap = np.array([[0, -1], [-1, 1]], dtype=np.int16)      # disjoint availability
    m = H.from_arrays(hap, np.array([0.0, 0.0]), np.array([1, 2]))   # dcM = 0
    assert m.switch_mass[0] == 0.0
    psi = np.ascontiguousarray(m.allowed.T.astype(np.float64))
    try:
        F.forward(m, psi)
    except ValueError as e:
        assert "zero total mass" in str(e) and "Delta_x = 0" in str(e)
    else:
        raise AssertionError("a dead-end chain returned a forward pass")


def test_the_fast_and_dense_forward_agree_away_from_the_toy():
    """The rank-one algebra at a K and T the toy does not reach, on random data."""
    rng = np.random.default_rng(2026)
    K, T = 60, 40
    hap = rng.integers(0, 3, size=(K, T)).astype(np.int16)
    hap[rng.random((K, T)) < 0.05] = -1
    hap[:, (hap != -1).sum(axis=0) == 0] = 0          # keep every position usable
    cm = np.sort(rng.random(T)) * 0.05
    m = H.from_arrays(hap, cm, np.arange(1, T + 1) * 1_000)
    util = U.beta_v((hap != -1).sum(axis=0).astype(np.int32))
    target = rng.integers(-1, 3, size=T).astype(np.int16)
    psi = U.tilt_potentials(m.allowed, hap, target, util.beta, C.eta_tau(0.5))
    fast = F.forward(m, psi)
    dense_alpha, dense_logc = E.dense_forward(hap, cm, psi)
    assert np.allclose(fast.alpha, dense_alpha, rtol=0, atol=1e-12)
    assert abs(fast.log_norm - dense_logc) < 1e-10


# ---------------------------------------------------------------------------
# Regressions from the 2026-09-30 independent audit. Each of these was a real
# defect that the 88-test suite passed over, so each gets a test that fails if it
# returns. Three were verified by mutation before the fix was accepted.
# ---------------------------------------------------------------------------

def test_calibration_is_checked_against_the_table_parsed_from_the_glossary():
    """Parse `docs/terms.md`'s table, then check the code against it.

    ⚠ THIS TEST EXISTS BECAUSE ITS ABSENCE HID A BUG. `calibration.py` claimed the
    docs contract already checked it against the glossary; it does not -- that test
    recomputes from `math.atanh` and never imports the module, and the test above
    compares against a hardcoded third copy. Three independent transcriptions of
    one table, no edge between the code and the document. Parsing the document is
    the only version of this check that can actually fail.
    """
    rows = []
    for line in (REPO / "docs" / "terms.md").read_text().splitlines():
        cells = [c.strip().strip("`") for c in line.strip().strip("|").split("|")]
        if len(cells) != 4:
            continue
        try:
            rows.append(tuple(float(c) for c in cells))
        except ValueError:
            continue
    assert len(rows) >= 7, f"parsed {len(rows)} calibration rows from terms.md, expected >= 7"
    for tau, eta, eps, ps in rows:
        assert abs(C.eta_tau(tau) - eta) < 5e-7, f"eta_tau({tau})"
        assert abs(C.epsilon_tau(tau) - eps) < 5e-7, f"epsilon_tau({tau})"
        assert abs(C.p_succ_bound(tau) - ps) < 5e-7, f"p_succ_bound({tau})"


def test_epsilon_tau_does_not_divide_by_the_utility_range():
    """REGRESSION 2026-09-30. `epsilon_tau` returned `2*arctanh(tau)/delta_u`.

    `delta_u` CANCELS in the budget: a utility spanning `delta_u` is tilted by an
    exponent whose range is `eta_tau * delta_u = arctanh(tau)`, so the pointwise
    ratio spans `exp(2*arctanh(tau))` whatever the range. The old code returned
    exactly `arctanh(tau)` at `delta_u` = 2 -- the expression `docs/terms.md`
    forbids in capitals as a silent understatement of the leakage. It was
    unreachable and untested at any non-unit range, so the only scenario the
    parameter exists for was the one that corrupted the headline privacy number.
    """
    for tau in (0.01, 0.1, 0.25, 0.5, 0.9):
        want = 2 * math.atanh(tau)
        for du in (0.5, 1.0, 2.0, 7.0):
            assert abs(C.epsilon_tau(tau, delta_u=du) - want) < 1e-15, (tau, du)
        assert abs(C.epsilon_tau(tau) - 2 * C.eta_tau(tau)) < 1e-15, (
            "the 2x relation must still hold at the unit range")


def test_a_diploid_release_brackets_the_SUM_of_its_two_log_z_values():
    """REGRESSION 2026-09-30, VERIFIED BY MUTATION. Nothing checked the sum.

    Each haplotype's `log_z` was asserted against `eta_tau/2`, and `release()`
    then returned their sum unchecked. Measured on chr21 before the fix: restoring
    the full `eta_tau` per haplotype -- the SUM aggregator in disguise -- was
    ACCEPTED at tau = 0.25, 0.5 and 0.9, returning a record whose published
    `log_z_p` (0.763 at tau = 0.5) violated its own stated invariant of
    [0, 0.549]. After the fix all three raise.

    Here the same hole is provoked without touching the source: two haplotype
    chains each legitimately bracketed at the FULL eta can sum past it.
    """
    t, util, m = _toy_bundle()
    tau = 0.5
    eta = C.eta_tau(tau)
    dip = np.stack([t.haplotypes[0], t.haplotypes[1]], axis=1)
    halves = []
    for h in (0, 1):
        psi = U.tilt_potentials(m.allowed, t.haplotypes, dip[:, h], util.beta, eta)
        halves.append(F.log_z_p(F.forward(m, psi), F.baseline_log_norm(m), eta))
    assert sum(halves) > eta, (
        f"the two full-exponent log_z values sum to {sum(halves)}, inside "
        f"[0, {eta}] -- this toy can no longer provoke the hole, so the test is "
        f"no longer testing anything")
    r = F.release(m, t.haplotypes, dip, util, tau, seed=3, target_id="dip")
    assert 0.0 <= r.log_z_p <= eta + 1e-9, "the correct diploid release must bracket"


def test_diploid_haplotype_streams_do_not_overlap_the_next_release():
    """REGRESSION 2026-09-30, MEASURED ON chr21. `seed + h` coupled the draws.

    Deriving the two haplotype streams as `seed` and `seed + 1` made
    `release(seed=s)`'s haplotype 1 BIT-IDENTICAL to `release(seed=s+1)`'s
    haplotype 0 -- 0 of 305,886 positions differed. So the first chr21 run's
    "3 independent draws per tau" were 4 correlated haplotype draws, which
    inflated the observed seed-to-seed spread 2.2x and made it coincidentally
    match the haploid `sd(u)`. `SeedSequence.spawn` gives provably disjoint
    streams from the one integer the record carries.
    """
    t, util, m = _toy_bundle()
    dip = np.stack([t.haplotypes[0], t.haplotypes[1]], axis=1)
    # A longer toy chain, so two independent draws agreeing everywhere by chance
    # is not the expected outcome.
    rng = np.random.default_rng(4)
    K, T = 12, 400
    hap = rng.integers(0, 2, size=(K, T)).astype(np.int16)
    cm = np.sort(rng.random(T)) * 0.5
    mm = H.from_arrays(hap, cm, np.arange(1, T + 1) * 500)
    uu = U.beta_v((hap != -1).sum(axis=0).astype(np.int32))
    tt = np.stack([hap[0], hap[1]], axis=1)
    a = F.release(mm, hap, tt, uu, 0.5, seed=500, target_id="a")
    b = F.release(mm, hap, tt, uu, 0.5, seed=501, target_id="b")
    for i in (0, 1):
        for j in (0, 1):
            assert not np.array_equal(a.state_path[:, i], b.state_path[:, j]), (
                f"release(500) haplotype {i} is identical to release(501) "
                f"haplotype {j}; the streams overlap again")


def test_the_allele_matrix_is_bound_to_the_model_it_was_built_from():
    """REGRESSION 2026-09-30. `release()` took the alleles as a free argument.

    `CohortHMM` carries only the availability mask, never the alleles, so this
    argument was the one place a caller could hand the sampler a DIFFERENT cohort
    of the same shape -- a target-pruned matrix above all, which is the standing
    constraint 1 violation the cohort/target directory split exists to prevent and
    the one shape that split cannot catch. The post-draw `-1` check only caught the
    subset that happened to be occupied.
    """
    t, util, m = _toy_bundle()
    pruned = t.haplotypes.copy()
    pruned[2, 0] = -1                      # prune one donor at one position
    try:
        F.release(m, pruned, t.haplotypes[0], util, 0.25, seed=1, target_id="x")
    except ValueError as e:
        assert "different cohorts" in str(e)
    else:
        raise AssertionError("a pruned allele matrix was accepted")
    try:
        F.release(m, t.haplotypes[:, :-1], t.haplotypes[0], util, 0.25,
                  seed=1, target_id="x")
    except ValueError as e:
        assert "different cohorts" in str(e) or "chain_span" in str(e)
    else:
        raise AssertionError("a mis-shaped allele matrix was accepted")


def test_phi_v_is_the_authority_for_the_indicator_the_fast_paths_inline():
    """REGRESSION 2026-09-30. `phi_v` is tested but never CALLED by the mechanism.

    `tilt_potentials` and `u_of_output` both re-inline the model-C indicator as a
    vectorised expression for speed, so the function carrying the documented
    semantics -- and the only one with a test -- was on no live code path. That is
    two implementations of one definition with no edge between them. This test is
    that edge: it checks the inlined comparisons agree with `phi_v` cell by cell,
    including the missing-target and missing-donor cases that are the whole reason
    `phi_v` has a docstring.
    """
    t, util, m = _toy_bundle()
    eta = C.eta_tau(0.5)
    for target in list(TOY.toy_inputs(t, kind="all")[::37]) + [t.haplotypes[2]]:
        psi = U.tilt_potentials(m.allowed, t.haplotypes, target, util.beta, eta)
        for v in range(t.T):
            want = U.phi_v(t.haplotypes[:, v], int(target[v]))
            got = np.where(m.allowed[:, v],
                           np.round(np.log(np.maximum(psi[v], 1e-300))
                                    / (eta * util.beta[v]), 6), want)
            assert np.allclose(got[m.allowed[:, v]], want[m.allowed[:, v]], atol=1e-5), (
                f"tilt_potentials disagrees with phi_v at position {v}")
        u_inline = float(U.u_of_output(
            t.haplotypes[0], target, util.beta)[0])
        u_phi = float(sum(util.beta[v] * U.phi_v(t.haplotypes[0:1, v], int(target[v]))[0]
                          for v in range(t.T)))
        assert abs(u_inline - u_phi) < 1e-12, "u_of_output disagrees with phi_v"


def test_the_two_independent_utility_implementations_agree():
    """`exact.utility_vector` and `utility.u_of_output` must compute one quantity.

    REGRESSION 2026-09-30. The gate's entire 3.5M-pair sweep scores with
    `exact.py`'s utility while every release scores with `utility.py`'s, and their
    null guards differ (`target != -1` against `target >= 0 & output >= 0`). If
    they diverged the gate would silently be measuring a different quantity than
    the mechanism optimises -- and nothing compared them.
    """
    t, util, m = _toy_bundle()
    base = E.enumerate_baseline(t.haplotypes, t.genetic_pos)
    worst = 0.0
    for target in TOY.toy_inputs(t, kind="all")[::11]:
        a = E.utility_vector(base, target, util.beta)
        b = U.u_of_output(base.paths, target, util.beta)
        worst = max(worst, float(np.abs(a - b).max()))
    assert worst < 1e-15, f"the two utility implementations differ by {worst}"


def test_state_marginals_match_brute_force_enumeration():
    """The exact posterior marginals, against the repo's own O(K^T) enumerator.

    ADDED 2026-10-03 after an independent audit found that 50 sampling draws had
    been spent estimating a quantity that is a LINEAR FUNCTIONAL of these
    marginals and therefore exact. The enumerator collapses the full joint over
    K^T state paths; `state_marginals` does it in O(TK) with the rank-one
    recursion and never materialises `A`. They must agree.
    """
    t, util, m = _toy_bundle()
    target, _ = TOY.sharp_pair(t)
    for tau in (0.0, 0.25, 0.5, 0.9):
        psi = U.tilt_potentials(m.allowed, t.haplotypes, target, util.beta, C.eta_tau(tau))
        gamma = F.state_marginals(m, F.forward(m, psi), psi)
        assert gamma.shape == (t.T, t.K) and gamma.dtype == np.float64
        assert np.allclose(gamma.sum(axis=1), 1.0, rtol=0, atol=1e-12)

        post = E.enumerate_state_posterior(t.haplotypes, t.genetic_pos, psi)
        ref = np.zeros((t.T, t.K))
        for i, pr in enumerate(post):
            if pr == 0.0:
                continue
            z, r = [], i
            for _ in range(t.T):
                z.append(r % t.K)
                r //= t.K
            for tt, st in enumerate(reversed(z)):
                ref[tt, st] += pr
        assert np.abs(gamma - ref).max() < 1e-12, (
            f"tau={tau}: marginals differ from enumeration by "
            f"{np.abs(gamma - ref).max():.3e}")

        # and the quantity they exist for, both ways
        known = target >= 0
        mism = (t.haplotypes != target[None, :]) & known[None, :]
        assert abs(F.expected_mismatches(m, t.haplotypes, target, gamma)
                   - float((ref.T * mism).sum())) < 1e-12


def test_expected_mismatches_is_exact_where_sampling_was_hopeless():
    """E[differ] is exact and LINEAR in `eta_tau`; no sample size was ever needed.

    REGRESSION 2026-10-03. The sampling route produced a signed monotonic pattern
    at n = 10 that reversed at n = 40, because the true effect on chr21 is about
    1.8 cells of 507,692 against a paired sd near 56. The exact route settles it
    in two forward and two backward passes. This test pins the linearity, which is
    what makes the slope quotable as a single number.
    """
    t, util, m = _toy_bundle()
    target, _ = TOY.sharp_pair(t)

    def differ(eta):
        psi = U.tilt_potentials(m.allowed, t.haplotypes, target, util.beta, eta)
        return F.expected_mismatches(
            m, t.haplotypes, target, F.state_marginals(m, F.forward(m, psi), psi))

    d0 = differ(0.0)
    for eta in (0.05, 0.1, 0.25):
        assert differ(eta) < d0, "the tilt must REDUCE expected mismatches"
    # a doubled exponent must double the effect -- the error the sampling script made
    one, two = d0 - differ(0.2), d0 - differ(0.4)
    assert abs(two / one - 2.0) < 0.05, (
        f"doubling eta changed the effect by {two/one:.3f}x, not ~2x; the effect is "
        f"no longer linear and the recorded slope is not quotable")


if __name__ == "__main__":
    tests = [(n, f) for n, f in sorted(globals().items())
             if n.startswith("test_") and callable(f)]
    failed = 0
    for name, fn in tests:
        try:
            fn()
            print(f"  PASS  {name}")
        except Exception as e:
            failed += 1
            print(f"  FAIL  {name}\n          {type(e).__name__}: {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
