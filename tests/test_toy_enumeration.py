"""[A-THY-toy-enumeration] The only test that can falsify our Theorem 1.

At chr21 scale `tv_empirical` is biased UPWARD by the plug-in estimator over a
large support, so a big measured value there is the estimator misbehaving and a
small one proves nothing at all. Falsification needs a support small enough to
enumerate exactly, which is what `mechanism/toy.py` is for. It reads nothing from
disk, so this gate is independent of the graph, the cohort, the genetic map, every
PanMixer artifact, and whether `data/` exists.

THE FIVE EXIT-GATE CLAUSES, and where each lives:

    (a) the sampler matches exact enumeration  -> test_ffbs_matches_*
    (b) exact TV <= tau over a tau grid and
        many input pairs                       -> test_exact_tv_never_exceeds_tau
    (c) log_z_p in [0, eta_tau] on every draw  -> test_log_z_p_*
    (d) tau = 0 reproduces baseline_model      -> test_tau_zero_*
    (e) the O(TK) path agrees with an O(TK^2)
        dense reference                        -> test_fast_forward_*

⚠ ONE CLAUSE WAS ADDED, BECAUSE (b) TURNS OUT TO BE A WEAK GATE. The worst exact
TV on this toy reaches only about a third of `tau`, so TV cannot detect a tilt
applied with `epsilon_tau` in place of `eta_tau` -- twice the intended exponent,
a one-token error, and the highest-severity error class in the project. The
POINTWISE RATIO SPREAD does detect it, and is attained exactly by a constructed
pair, so the gate is sharp rather than merely satisfied:
`test_pointwise_ratio_spread_*`, `test_the_sharp_pair_saturates_*`, and
`test_a_doubled_tilt_exponent_is_caught` which injects the bug and checks which
gates fire.

KILL CRITERION: any exact-TV value above `tau` on the toy support.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from mechanism import calibration as C   # noqa: E402
from mechanism import exact as E         # noqa: E402
from mechanism import ffbs as F          # noqa: E402
from mechanism import hmm as H           # noqa: E402
from mechanism import toy as TOY         # noqa: E402
from mechanism import utility as U       # noqa: E402

# The grid. 0 is the floor (an untilted draw from R_D) and 0.9 is far past
# anything defensible as a privacy claim -- included because a bound must hold
# where it is loose as well as where it is tight.
TAU_GRID = (0.0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.90)

_CACHE: dict = {}


def _toy():
    if "toy" not in _CACHE:
        t = TOY.toy_cohort()
        _CACHE["toy"] = t
        _CACHE["base"] = E.enumerate_baseline(t.haplotypes, t.genetic_pos)
        _CACHE["util"] = U.beta_v(t.support)
        _CACHE["hmm"] = H.from_arrays(t.haplotypes, t.genetic_pos, t.positions,
                                      cohort_id=t.name)
    return _CACHE["toy"], _CACHE["base"], _CACHE["util"], _CACHE["hmm"]


def _distinct_utility_vectors(base, beta, inputs):
    """Dedupe inputs by their utility VECTOR over the support.

    `Q_p` depends on `p` only through `u(p, .)`, so two inputs with the same
    utility vector induce the identical release distribution and contribute
    nothing new to an all-pairs sweep. Deduping is exact, not a sample: the
    sweep below still covers every pair of admissible inputs, at a fraction of
    the cost.
    """
    seen: dict[bytes, np.ndarray] = {}
    for p in inputs:
        u = E.utility_vector(base, p, beta)
        seen.setdefault(np.round(u, 15).tobytes(), p)
    return np.stack(list(seen.values())), len(seen)


def _worst_pair(mats, fn):
    """Sweep every unordered pair of rows, returning (worst value, i, j).

    Chunked row-against-remainder so nothing of size P^2 x M is ever allocated.
    """
    P = len(mats)
    worst, wi, wj = -np.inf, -1, -1
    for i in range(P - 1):
        vals = fn(mats[i], mats[i + 1:])
        j = int(np.argmax(vals))          # over PAIRS, not over a state space
        if vals[j] > worst:
            worst, wi, wj = float(vals[j]), i, i + 1 + j
    return worst, wi, wj


def _tv_row(qi, rest):
    return 0.5 * np.abs(rest - qi[None, :]).sum(axis=1)


def _spread_row(qi, rest):
    r = qi[None, :] / rest
    return r.max(axis=1) / r.min(axis=1)


# ---------------------------------------------------------------------------
# The fixture itself. A toy that has quietly stopped exercising a feature still
# reports green, so its structure is asserted before anything is measured on it.
# ---------------------------------------------------------------------------

def test_the_toy_is_small_enough_to_enumerate_and_rich_enough_to_matter():
    t, base, util, m = _toy()
    assert base.n_enumerated == t.K ** t.T == 16_384
    assert base.n_state_paths < base.n_enumerated, (
        "every state path has positive mass, so neither RENORMALISE nor the "
        "forced-stay step is constraining anything")
    assert base.M < base.n_state_paths, (
        "no two state paths collapse to one output path. The D1/D3 collision is "
        "gone, and TV is being measured over states rather than over releases.")
    assert abs(base.prob.sum() - 1.0) < 1e-12
    assert base.prob.min() > 0, "Y_D contains a path with zero probability"
    assert 0.0 in set(m.switch_mass.tolist()), "the forced-stay step is gone"
    assert len(set(util.beta.tolist())) > 1, "beta_v is uniform; R3 is untested"
    print(f"    K={t.K} T={t.T}: {base.n_enumerated} state paths -> "
          f"{base.n_state_paths} with mass -> {base.M} output paths in Y_D")


def test_output_support_is_target_independent():
    """Y_D is a function of G and D only -- standing constraint 1, measured.

    A target-dependent support is the §4.4 disjoint-support argument that voids
    the theorem, and it is the single rule `CLAUDE.md` says is most likely to be
    broken by accident. Here it is checked the only way that can actually fail:
    run the mechanism over every admissible input at every `tau` and confirm the
    set of outputs that can be emitted never moves.
    """
    t, base, util, m = _toy()
    reachable = {row.tobytes() for row in base.paths}
    for tau in TAU_GRID:
        eta = C.eta_tau(tau)
        for p in TOY.toy_inputs(t, kind="called"):
            u = E.utility_vector(base, p, util.beta)
            q, _ = E.release_distribution(base, u, eta)
            assert (q > 0).all(), (
                f"tau={tau} zeroed part of the support; the tilt must reweight "
                f"Y_D, never prune it")
            assert {row.tobytes() for row in base.paths[q > 0]} == reachable


def test_a_forced_stay_step_really_forbids_switching():
    """Delta_x = 0 gives P(switch) = 0, so Y_D is strictly smaller. Both halves matter.

    This is the behaviour `chain_span` in `docs/terms.md` calls CORRECT and
    deliberately leaves in place: 618 interior chr21 map intervals have dcM = 0,
    and between two variants with no genetic separation, not switching is the
    model rather than a defect. The mechanism must honour it at every `tau` --
    a tilt that could push mass onto an unreachable path would be emitting a
    path `R_D` assigns zero probability, which is outside the support the
    theorem is stated over.
    """
    t, base, util, m = _toy()
    flat = int(np.flatnonzero(m.switch_mass == 0)[0])
    assert np.diff(t.genetic_pos)[flat] == 0

    # D1 and D2 differ at every position, so "D1 then D2" across the flat step
    # is an output path only reachable by switching there.
    forbidden = t.haplotypes[0].copy()
    forbidden[flat + 1:] = t.haplotypes[1][flat + 1:]
    assert forbidden.tobytes() not in {row.tobytes() for row in base.paths}, (
        "a path that switches donors across a zero-distance step is in Y_D")

    fwd_seeds = range(400)
    m_ = m
    psi = U.tilt_potentials(m_.allowed, t.haplotypes, t.haplotypes[1], util.beta,
                            C.eta_tau(0.9))
    fwd = F.forward(m_, psi)
    for s in fwd_seeds:
        z = F.backward_sample(m_, fwd, np.random.default_rng(s))
        assert z[flat] == z[flat + 1], (
            f"seed {s} switched donors across a zero-distance step, from state "
            f"{z[flat]} to {z[flat+1]}")


# ---------------------------------------------------------------------------
# (c) and (d): the definition of the release distribution, with no sampler in it.
# ---------------------------------------------------------------------------

def test_log_z_p_is_inside_its_bracket_for_every_input_and_every_tau():
    """0 <= log Z_p <= eta_tau, by direct summation over the support.

    Computed from the enumeration rather than from a forward pass, so this tests
    the DEFINITION of `Q^tau_p` and not an implementation of it. `Z_p` is a
    convex combination of `exp(eta u)` with `u` in [0,1], so it lies in
    [1, e^eta] -- if it does not, the release distribution is not the object the
    theorem is about.
    """
    t, base, util, m = _toy()
    worst = 0.0
    for tau in TAU_GRID:
        eta = C.eta_tau(tau)
        for p in TOY.toy_inputs(t, kind="all"):
            u = E.utility_vector(base, p, util.beta)
            _, lz = E.release_distribution(base, u, eta)
            assert -1e-12 <= lz <= eta + 1e-12, (
                f"log Z_p = {lz} outside [0, {eta}] at tau={tau}, input={p}")
            if eta > 0:
                worst = max(worst, lz / eta)
    print(f"    worst log_z_p reached {worst:.4f} of eta_tau")


def test_tau_zero_reproduces_the_baseline_exactly():
    """The floor column. `tau` = 0 means the release ignores the target entirely.

    Not "approximately R_D" -- exactly, to floating point, for every input,
    because `eta_tau` = 0 makes every tilt potential exactly 1. A floor that is
    merely close to `R_D` would mean the mechanism has a target-dependent
    component outside the tilt.
    """
    t, base, util, m = _toy()
    for p in TOY.toy_inputs(t, kind="all"):
        u = E.utility_vector(base, p, util.beta)
        q, lz = E.release_distribution(base, u, 0.0)
        assert np.allclose(q, base.prob, rtol=0, atol=1e-15)
        assert abs(lz) < 1e-15


def test_the_floor_carries_no_target_information_but_looks_like_a_good_release():
    """Why `target_fidelity` has to be its own axis, demonstrated.

    A `tau` = 0 draw is a genuine mosaic of cohort haplotypes. It scores well on
    any cohort-level utility metric while carrying zero information about the
    target -- so quoting the floor without a fidelity axis makes it look like a
    free lunch. `docs/plan.md` records this as a MUST for the head-to-head; the
    toy shows it is a property of the mechanism, not of that comparison.
    """
    t, base, util, m = _toy()
    p, _ = TOY.sharp_pair(t)
    u = E.utility_vector(base, p, util.beta)
    q0, _ = E.release_distribution(base, u, 0.0)
    expected_u_at_floor = float(q0 @ u)
    q9, _ = E.release_distribution(base, u, C.eta_tau(0.9))
    assert float(q9 @ u) > expected_u_at_floor, "the tilt does not raise utility"
    assert 0.0 < expected_u_at_floor < 1.0
    print(f"    E[u] at the floor {expected_u_at_floor:.4f} -> "
          f"{float(q9 @ u):.4f} at tau=0.9 (ceiling {u.max():.4f})")


# ---------------------------------------------------------------------------
# (b) the TV bound, and the sharper pointwise statement it follows from.
# ---------------------------------------------------------------------------

def test_exact_tv_never_exceeds_tau():
    """THE KILL CRITERION. Exact TV over every pair of admissible inputs, per tau.

    Exhaustive, not sampled: inputs are deduped by their utility vector, which is
    exact because `Q_p` depends on `p` only through that vector, so every pair of
    admissible inputs is covered.
    """
    t, base, util, m = _toy()
    inputs, n_dedup = _distinct_utility_vectors(base, util.beta,
                                                TOY.toy_inputs(t, kind="all"))
    print(f"    {len(TOY.toy_inputs(t, kind='all'))} admissible inputs -> "
          f"{n_dedup} distinct release distributions -> "
          f"{n_dedup*(n_dedup-1)//2} pairs per tau")
    for tau in TAU_GRID:
        eta = C.eta_tau(tau)
        Q = np.stack([E.release_distribution(base, E.utility_vector(base, p, util.beta),
                                             eta)[0] for p in inputs])
        worst, i, j = _worst_pair(Q, _tv_row)
        assert worst <= tau + 1e-12, (
            f"EXACT TV {worst} EXCEEDS tau={tau} for inputs\n"
            f"  {inputs[i]}\n  {inputs[j]}\nTheorem 1 is falsified as implemented.")
        print(f"    tau={tau:<5} worst exact TV {worst:.6f} "
              f"({0.0 if tau == 0 else worst/tau:.1%} of tau)")


def test_pointwise_ratio_spread_never_exceeds_exp_epsilon_tau():
    """The sharp gate: `max_y Q_p/Q_q  /  min_y Q_p/Q_q  <=  exp(epsilon_tau)`.

    The ratio is `exp(eta_tau (u_p - u_q)) * Z_q/Z_p`, and the `Z` factor is
    common to every `y`, so it cancels out of the spread entirely -- leaving
    `exp(eta_tau * range(u_p - u_q))`, bounded by `exp(2 eta_tau)` because each
    utility is in [0,1]. Nothing is estimated and nothing cancels approximately,
    which is why this is the gate with power.
    """
    t, base, util, m = _toy()
    inputs, _ = _distinct_utility_vectors(base, util.beta,
                                          TOY.toy_inputs(t, kind="all"))
    for tau in TAU_GRID:
        eta, eps = C.eta_tau(tau), C.epsilon_tau(tau)
        Q = np.stack([E.release_distribution(base, E.utility_vector(base, p, util.beta),
                                             eta)[0] for p in inputs])
        worst, i, j = _worst_pair(Q, _spread_row)
        assert worst <= math.exp(eps) * (1 + 1e-12), (
            f"ratio spread {worst} exceeds exp(epsilon_tau)={math.exp(eps)} at "
            f"tau={tau} for inputs\n  {inputs[i]}\n  {inputs[j]}")
        print(f"    tau={tau:<5} worst spread {worst:.6f} vs "
              f"exp(eps_tau) {math.exp(eps):.6f}")


def test_the_sharp_pair_saturates_the_ratio_bound():
    """The bound is ATTAINED, not approached -- so any inflated tilt fails it.

    `p` = D1's alleles and `q` = D2's differ at every position, and copying
    either donor throughout is in the support, so `u_p - u_q` is exactly +1 at
    one output and -1 at another and the spread is exactly `exp(2 eta_tau)`.
    A gate that is never approached cannot detect anything; this one sits on the
    boundary at every `tau`.
    """
    t, base, util, m = _toy()
    p, q = TOY.sharp_pair(t)
    for tau in TAU_GRID[1:]:
        eta, eps = C.eta_tau(tau), C.epsilon_tau(tau)
        up = E.utility_vector(base, p, util.beta)
        uq = E.utility_vector(base, q, util.beta)
        assert abs(up.max() - 1.0) < 1e-12 and abs(uq.max() - 1.0) < 1e-12
        assert up.min() < 1e-12 and uq.min() < 1e-12
        qp, _ = E.release_distribution(base, up, eta)
        qq, _ = E.release_distribution(base, uq, eta)
        got = E.ratio_spread(qp, qq)
        assert abs(got - math.exp(eps)) < 1e-9 * math.exp(eps), (
            f"tau={tau}: spread {got} != exp(epsilon_tau) {math.exp(eps)}; the "
            f"sharp pair no longer saturates, so this gate has lost its power")


def test_a_doubled_tilt_exponent_is_caught_by_exactly_one_gate():
    """A MUTATION TEST, and the reason the gate list needed a sixth clause.

    Inject the likeliest calibration error -- `epsilon_tau = 2 * eta_tau`, so
    reaching for the wrong one doubles the tilt exponent, a one-token mistake and
    the highest-severity error class in the project -- then record which gates
    notice. A suite whose gates all pass a known-broken mechanism is worse than
    no suite, so what each one does is asserted rather than assumed.

    MEASURED ON THIS TOY, and it is not what was expected:

      | gate                      | tau 0.01-0.25 | tau 0.5-0.9 |
      |---------------------------|---------------|-------------|
      | exact TV <= tau           | PASSES        | PASSES      |
      | log_z_p in [0, eta_tau]   | PASSES        | fails       |
      | ratio spread <= exp(eps)  | fails         | fails       |

    So the pre-registered exit-gate clause (b) has NO power against this error at
    any `tau`, and clause (c) has none below about 0.5 -- which is the whole range
    anyone would defend as a privacy claim, since `tau` = 0.10 already bounds
    equal-prior identification at 0.55. The reason is arithmetic:
    `log Z_p` is about `eta_used * E_R[u]`, and `E_R[u]` is 0.455 here, so a
    doubled exponent still lands at 0.91 of the bracket. Only the pointwise ratio
    spread fires everywhere, because the normaliser cancels out of it and there
    is nothing left for a mis-scaled exponent to hide behind.

    ⚠ If the first column ever starts failing, the finding recorded in
    `docs/memory.md` for 2026-09-30 is what has gone stale, not this test.
    """
    t, base, util, m = _toy()
    p, q = TOY.sharp_pair(t)
    for tau in (0.01, 0.05, 0.10, 0.25, 0.50, 0.90):
        eta, eps = C.eta_tau(tau), C.epsilon_tau(tau)
        up = E.utility_vector(base, p, util.beta)
        uq = E.utility_vector(base, q, util.beta)
        bad_p, bad_lz = E.release_distribution(base, up, eps)   # the bug: eps for eta
        bad_q, _ = E.release_distribution(base, uq, eps)
        tv, spread = E.tv(bad_p, bad_q), E.ratio_spread(bad_p, bad_q)

        assert tv <= tau, (
            f"tau={tau}: TV now catches a doubled exponent ({tv} > {tau}). Good "
            f"news, but the recorded finding says otherwise -- update it.")
        assert spread > math.exp(eps) * (1 + 1e-9), (
            f"tau={tau}: the ratio spread {spread} did NOT catch a doubled "
            f"exponent. That was the only gate with power; the suite is now blind.")
        caught_by_logz = bad_lz > eta
        assert caught_by_logz == (tau >= 0.5), (
            f"tau={tau}: log_z_p caught={caught_by_logz}, expected "
            f"{tau >= 0.5}. The bracket's power against this error changed.")
        print(f"    tau={tau:<5} doubled exponent: TV {tv:.6f} (<= {tau}, PASSES) "
              f"| log_z_p {bad_lz:.6f} vs eta {eta:.6f} "
              f"({'FAILS' if caught_by_logz else 'PASSES'}) "
              f"| spread {spread:.4f} vs {math.exp(eps):.4f} (FAILS)")


# ---------------------------------------------------------------------------
# (a) and (e): the sampler, against the object it claims to sample from.
# ---------------------------------------------------------------------------

def test_fast_forward_matches_the_dense_reference():
    """O(TK) rank-one stepping against an O(TK^2) materialised matrix.

    `A_t = (1 - s_t) I + (s_t / K) J`, so the forward step needs one scalar sum
    instead of a matrix-vector product. The dense reference in `exact.py` writes
    the transition from the published closed form and never imports `hmm.py`, so
    the two can disagree -- which is the property that makes agreement evidence.
    """
    t, base, util, m = _toy()
    for tau in TAU_GRID:
        for p in TOY.toy_inputs(t, kind="called")[::17]:
            psi = U.tilt_potentials(m.allowed, t.haplotypes, p, util.beta,
                                    C.eta_tau(tau))
            fast = F.forward(m, psi)
            dense_alpha, dense_logc = E.dense_forward(t.haplotypes, t.genetic_pos, psi)
            assert np.allclose(fast.alpha, dense_alpha, rtol=0, atol=1e-13)
            assert abs(fast.log_norm - dense_logc) < 1e-12


def test_forward_pass_log_z_p_matches_the_enumerated_value():
    """The forward pass and the brute-force sum agree on `log Z_p`.

    Two independent routes to one number: `sum_t log c_t` minus the baseline's
    own conditioning constant, against `log sum_y R_D(y) exp(eta u)` over the
    enumerated support. This is the check that the subtraction in `ffbs.log_z_p`
    is the right one -- without it the forward value is 0.575 out on the toy and
    584.81 out on chr21 (measured 2026-09-30), against a bracket 0.549 wide at
    tau = 0.5.
    """
    t, base, util, m = _toy()
    baseline = F.baseline_log_norm(m)
    for tau in TAU_GRID:
        eta = C.eta_tau(tau)
        for p in TOY.toy_inputs(t, kind="all")[::7]:
            psi = U.tilt_potentials(m.allowed, t.haplotypes, p, util.beta, eta)
            got = F.log_z_p(F.forward(m, psi), baseline, eta)
            want = E.release_distribution(base, E.utility_vector(base, p, util.beta),
                                          eta)[1]
            assert abs(got - want) < 1e-12, (
                f"forward log_z_p {got} != enumerated {want} at tau={tau}")


def test_ffbs_matches_the_exact_release_distribution():
    """(a) The sampler draws from `Q^tau_p`, checked against enumeration.

    Deterministic: one fixed seed per draw, so the numbers below are reproducible
    and this is not a flaky test. Two statistics, because either alone can be
    fooled -- a per-path 5-sigma z-score catches a single mis-weighted output,
    and the total variation against its own analytic expectation catches a
    diffuse bias that no single path makes significant.
    """
    t, base, util, m = _toy()
    n_draws = 40_000
    index = {row.tobytes(): i for i, row in enumerate(base.paths)}

    for tau in (0.0, 0.25, 0.9):
        eta = C.eta_tau(tau)
        p, _ = TOY.sharp_pair(t)
        psi = U.tilt_potentials(m.allowed, t.haplotypes, p, util.beta, eta)
        fwd = F.forward(m, psi)
        want, _ = E.release_distribution(base, E.utility_vector(base, p, util.beta), eta)

        counts = np.zeros(base.M, dtype=np.int64)
        rng = np.random.default_rng(20260930)
        for _ in range(n_draws):
            z = F.backward_sample(m, fwd, rng)
            y = t.haplotypes[z, np.arange(t.T)].astype(np.int16)
            k = y.tobytes()
            assert k in index, f"the sampler emitted a path outside Y_D: {y}"
            counts[index[k]] += 1
        got = counts / n_draws

        big = want * n_draws >= 10           # where the normal approximation holds
        se = np.sqrt(want[big] * (1 - want[big]) / n_draws)
        # ⚠ PER-PATH, then take the max. Dividing the largest deviation by the
        # SMALLEST standard error mixes two different paths and inflates the
        # z-score by an order of magnitude -- it reported 19.9 sigma on a sampler
        # the joint-posterior test had already shown to be correct.
        z_max = float((np.abs(got[big] - want[big]) / se).max())
        # E[TV] for a multinomial plug-in estimate, from the half-normal mean.
        tv_expected = 0.5 * np.sqrt(2 * want * (1 - want) / (math.pi * n_draws)).sum()
        tv_got = E.tv(got, want)

        assert z_max < 5.0, (
            f"tau={tau}: a released path is {z_max:.1f} sigma from its exact "
            f"probability over {n_draws} draws")
        assert tv_got < 3 * tv_expected, (
            f"tau={tau}: TV(empirical, exact) = {tv_got:.6f} against an expected "
            f"{tv_expected:.6f}; the sampler is biased, not noisy")
        print(f"    tau={tau:<5} {n_draws} draws: max z {z_max:.2f}, "
              f"TV(emp, exact) {tv_got:.6f} vs expected {tv_expected:.6f}")


def test_ffbs_matches_the_exact_state_path_posterior():
    """The JOINT over donor states, not just the output marginals.

    Output marginals can agree while the joint is wrong -- a sampler that drew
    each position independently from its posterior marginal would pass a marginal
    check and produce mosaics with the wrong recombination structure, which is
    precisely PanMixer's defect in the other direction. So the state-path joint
    is checked directly.
    """
    t, base, util, m = _toy()
    tau, n_draws = 0.5, 30_000
    eta = C.eta_tau(tau)
    p, _ = TOY.sharp_pair(t)
    psi = U.tilt_potentials(m.allowed, t.haplotypes, p, util.beta, eta)
    fwd = F.forward(m, psi)
    want = E.enumerate_state_posterior(t.haplotypes, t.genetic_pos, psi)

    counts = np.zeros(t.K ** t.T, dtype=np.int64)
    rng = np.random.default_rng(31415)
    powers = t.K ** np.arange(t.T - 1, -1, -1)
    for _ in range(n_draws):
        counts[int(F.backward_sample(m, fwd, rng) @ powers)] += 1
    got = counts / n_draws

    tv_expected = 0.5 * np.sqrt(2 * want * (1 - want) / (math.pi * n_draws)).sum()
    tv_got = E.tv(got, want)
    assert (got[want == 0] == 0).all(), (
        "the sampler visited a state path with zero posterior probability")
    assert tv_got < 3 * tv_expected, (
        f"TV(empirical, exact) over state paths = {tv_got:.6f} against an "
        f"expected {tv_expected:.6f}")
    print(f"    {int((want > 0).sum())} state paths with mass, {n_draws} draws: "
          f"TV {tv_got:.6f} vs expected {tv_expected:.6f}")


def test_the_release_never_carries_a_missing_allele():
    """RENORMALISE, end to end: the chain cannot occupy an unavailable donor.

    The toy has donors missing at two positions, so a sampler that treated `-1`
    as an ordinary allele would emit it. `ffbs.release` asserts this too; here it
    is checked over many seeds and every `tau`, including `tau` = 0 where the
    tilt cannot be masking the error.
    """
    t, base, util, m = _toy()
    missing_at = np.flatnonzero((t.haplotypes == -1).any(axis=0))
    assert len(missing_at) >= 2
    for tau in (0.0, 0.5):
        for seed in range(50):
            r = F.release(m, t.haplotypes, t.haplotypes[0], util, tau,
                          seed=seed, target_id="toyD1", graph_id=t.name)
            assert (r.released_path != -1).all()
            for v in missing_at:
                assert t.haplotypes[r.state_path[v], v] != -1


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
