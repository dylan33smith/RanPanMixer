"""Algorithm 1: forward filtering, backward sampling. The mechanism itself.

Draws EXACTLY from `Q^tau_p(y) = R_D(y) exp(eta_tau * u(p,y)) / Z_p`. Exact, not
approximate -- which is the only reason the released path carries a guarantee at
all. Nothing upstream has a backward pass; PanMixer contains none, so this is new
work rather than something inherited.

FOUR THINGS THIS FILE IS CAREFUL ABOUT, each a real failure it would otherwise have.

1. ⚠ SAMPLE; NEVER ARGMAX. Viterbi/MAP decoding returns the single
   highest-probability state sequence. It is not a draw from `Q^tau_p` and carries
   NO guarantee whatsoever, while producing output that looks entirely plausible.
   `viterbi_release` was retired in `docs/terms.md` before use for that reason. A
   test walks the AST of this package and fails on any use of `argmax`, because the
   failure is silent and the temptation is real in an HMM codebase.

2. ⚠ EVERY ALPHA COLUMN IS STORED, in float64. PanMixer keeps only the running
   vector, which is all a forward-only ancestral simulation needs; backward
   sampling needs the whole (T, K) array. float64 is LOAD-BEARING here, not a
   preference: under model C the per-position tilt at tau = 0.5 is 1.61e-06
   against float32's resolution of 1.19e-07, so only 13.5x the noise floor
   (a block partition would have had 45.7x). Rounding would consume a meaningful
   fraction of the target's entire influence on its own release. Storage is
   about 215 MB at chr21's T = 305,886 and K = 88.

3. ⚠ `A` IS NEVER MATERIALISED. `A_t = (1 - s_t) I + (s_t / K) J` is rank-one plus
   diagonal, so the forward step is one scalar sum and the backward draw is one
   vector update -- O(K) each. A dense matrix per position would be about 21 TB on
   chr21.

4. ⚠ `log_z_p` IS THE TILTED NORMALISER MINUS THE BASELINE'S OWN, and this is the
   trap. `missing_policy` = RENORMALISE enters as a 0/1 potential (see `hmm.py`
   for why it cannot be a renormalised transition row), so `sum_t log c_t` from
   the tilted pass contains the conditioning constant of `R_D` as well as the
   tilt. Reporting that sum as `log_z_p` puts it far outside `[0, eta_tau]` --
   MEASURED -584.81 on chr21 against a bracket of [0, 0.549] at tau = 0.5, so it
   is out by a factor of about a thousand -- and the natural "fix" of widening the
   assertion would throw away the cheapest correctness check in the project.
   Subtract the
   baseline pass instead. It is target-independent, so it is computed once per
   cohort and cached.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass, field
from datetime import date

import numpy as np

from . import calibration, utility
from .hmm import CohortHMM


@dataclass(frozen=True)
class Forward:
    """The stored forward pass: every alpha column, and the per-position normalisers."""

    alpha: np.ndarray      # (T, K) float64, each row summing to 1
    log_c: np.ndarray      # (T,) float64, per-position log normaliser

    @property
    def log_norm(self) -> float:
        """`sum_t log c_t`. ⚠ NOT `log_z_p` -- see point 4 in the module docstring."""
        return float(self.log_c.sum())


def forward(hmm: CohortHMM, potentials: np.ndarray) -> Forward:
    """The O(TK) forward pass, storing every column.

    `potentials` is (T, K) and must already carry BOTH the RENORMALISE indicator
    and the tilt; `utility.tilt_potentials` builds it that way deliberately, so
    that applying one and forgetting the other is not expressible.

    The recursion, with `e_t = 1 - s_t` and `unif_t = s_t / K`:

        pred_t(k) = e_t * alpha_{t-1}(k) + unif_t * sum_j alpha_{t-1}(j)
        r_t(k)    = pred_t(k) * psi_t(k)
        c_t       = sum_k r_t(k),   alpha_t = r_t / c_t

    The scalar sum is computed rather than assumed to be 1. It costs one
    reduction we are doing anyway, and it keeps the recursion correct rather than
    merely correct-given-an-invariant.
    """
    T, K = hmm.T, hmm.K
    if potentials.shape != (T, K):
        raise ValueError(f"potentials must be (T, K) = {(T, K)}, got {potentials.shape}")
    if potentials.dtype != np.float64:
        raise TypeError(
            f"potentials are {potentials.dtype}; float64 is load-bearing here "
            f"(the per-position tilt is only ~13.5x float32 resolution)")
    if (potentials < 0).any():
        raise ValueError("negative potential; psi = exp(...) times a 0/1 mask cannot be")

    alpha = np.zeros((T, K), dtype=np.float64)
    log_c = np.zeros(T, dtype=np.float64)

    r = hmm.rho * potentials[0]
    c = r.sum()
    if c <= 0:
        raise ValueError(
            "position 0 of the chain has zero total mass. Every donor is either "
            "excluded or annihilated by its potential -- the chain has nothing to "
            "start from.")
    alpha[0] = r / c
    log_c[0] = math.log(c)

    s = hmm.switch_mass
    for t in range(1, T):
        stay = 1.0 - s[t - 1]
        unif = s[t - 1] / K
        prev = alpha[t - 1]
        r = (stay * prev + unif * prev.sum()) * potentials[t]
        c = r.sum()
        if c <= 0:
            raise ValueError(
                f"chain position {t} (pos {hmm.positions[t]}) has zero total "
                f"mass: every state reachable from position {t-1} is excluded "
                f"there. With Delta_x = 0 the transition is the identity, so a "
                f"donor that becomes unavailable dead-ends -- which is why "
                f"RENORMALISE is a potential and not a renormalised row.")
        alpha[t] = r / c
        log_c[t] = math.log(c)

    return Forward(alpha=alpha, log_c=log_c)


def baseline_log_norm(hmm: CohortHMM) -> float:
    """`sum_t log c_t` for the UNTILTED chain: the conditioning constant of `R_D`.

    ⚠ Target-independent by construction -- the potentials are the allowed-state
    mask and nothing else -- so this is a property of the cohort and is cached by
    the caller across every target and every `tau`. It is also exactly what
    `tau = 0` must reproduce, which makes `log_z_p == 0` at `tau = 0` a free
    end-to-end check rather than a claim.
    """
    return forward(hmm, np.ascontiguousarray(hmm.allowed.T.astype(np.float64))).log_norm


def log_z_p(tilted: Forward, baseline_log_norm_: float, eta: float,
            *, tol: float = 1e-6) -> float:
    """`log Z_p`, range-asserted against `[0, eta_tau]`.

    Theory pins it: `Z_p = sum_y R_D(y) exp(eta_tau u(p,y))` with `u` in [0,1], so
    `Z_p` in `[1, e^eta]`. If this leaves the bracket the sampler is wrong -- it is
    the cheapest correctness assertion in the whole implementation, so it is
    asserted and not reported.

    The tolerance covers accumulated floating-point error in the difference of two
    large sums, not modelling slack. On chr21 each sum is about -1.2e4 over
    305,886 terms, so the difference carries roughly 1e-9 of error against a
    bracket of width 0.549 at `tau` = 0.5. A failure at 1e-6 is a real defect.
    """
    lz = tilted.log_norm - baseline_log_norm_
    if not (-tol <= lz <= eta + tol):
        raise AssertionError(
            f"log_z_p = {lz!r} is outside [0, eta_tau] = [0, {eta!r}].\n"
            f"  tilted sum_t log c_t   {tilted.log_norm!r}\n"
            f"  baseline sum_t log c_t {baseline_log_norm_!r}\n"
            f"The release distribution is not the one the theorem is about. Do "
            f"NOT widen this tolerance -- fix the sampler."
        )
    return float(lz)


def backward_sample(hmm: CohortHMM, fwd: Forward, rng: np.random.Generator
                    ) -> np.ndarray:
    """Draw one state path from the posterior. O(TK), and never an argmax.

    `z_T ~ Categorical(alpha_T)`, then `z_t | z_{t+1}` with weights
    `alpha_t(j) * A_{t+1}(j, z_{t+1})`. The rank-one structure makes that weight
    `alpha_t(j) * s/K`, with the single entry `j = z_{t+1}` getting an extra
    `alpha_t(j) * (1 - s)` -- so it is a scaled copy of `alpha_t` with one element
    bumped, never a matrix column lookup.
    """
    T, K = hmm.T, hmm.K
    z = np.empty(T, dtype=np.int32)
    z[T - 1] = _categorical(fwd.alpha[T - 1], rng)
    s = hmm.switch_mass
    for t in range(T - 2, -1, -1):
        unif = s[t] / K
        w = fwd.alpha[t] * unif
        w[z[t + 1]] += fwd.alpha[t][z[t + 1]] * (1.0 - s[t])
        total = w.sum()
        if total <= 0:
            raise AssertionError(
                f"backward weights vanish at chain position {t}; alpha_t carries "
                f"no mass that can reach state {z[t+1]} at {t+1}")
        z[t] = _categorical(w / total, rng)
    return z


def state_marginals(hmm: CohortHMM, fwd: Forward, potentials: np.ndarray) -> np.ndarray:
    """(T, K) posterior state marginals `gamma_t(k) = P(Z_t = k | the tilt)`. O(TK).

    ⚠ WHY THIS EXISTS, AND IT IS A LESSON THE PROJECT ALREADY WROTE DOWN ONCE.
    `docs/bugs.md` records, for expected utility, that the fix for a sampling
    experiment with no power is to "report the exact quantity instead:
    `E_Q[u] = d log Z_p / d eta`, three forward passes and no sampling". That
    lesson was then NOT applied to the metric that replaced it -- the count of
    released cells differing from the target -- and 50 draws were spent estimating
    something computable exactly. Found by an independent audit, 2026-10-03.

    The general fact: ANY quantity that is a sum over positions of a function of
    the state at that position is a linear functional of these marginals, so its
    expectation under the release distribution is exact with zero Monte-Carlo
    error. `expected_mismatches` below is one such quantity; so is expected
    `u_path`, and so is any per-position fidelity or rewrite accounting.

    The backward recursion uses the same rank-one decomposition as the forward
    pass, so `A` is never materialised here either:

        b_t(j) = (1 - s_t) * psi_{t+1}(j) * b_{t+1}(j)
                 + (s_t / K) * sum_k psi_{t+1}(k) * b_{t+1}(k)

    ⚠ `b_t` is rescaled to sum 1 at every position. Without it the product
    underflows float64 long before T = 305,886 -- the same reason the forward pass
    normalises. The scaling cancels when `gamma` is normalised, so it is free.
    """
    T, K = hmm.T, hmm.K
    if potentials.shape != (T, K):
        raise ValueError(f"potentials must be (T, K) = {(T, K)}, got {potentials.shape}")
    b = np.zeros((T, K), dtype=np.float64)
    b[T - 1] = 1.0 / K
    s = hmm.switch_mass
    for t in range(T - 2, -1, -1):
        v = potentials[t + 1] * b[t + 1]
        r = (1.0 - s[t]) * v + (s[t] / K) * v.sum()
        tot = r.sum()
        if tot <= 0:
            raise AssertionError(
                f"backward message vanishes at chain position {t}; no state there "
                f"can reach the end of the chain")
        b[t] = r / tot

    g = fwd.alpha * b
    z = g.sum(axis=1, keepdims=True)
    if (z <= 0).any():
        raise AssertionError("a position has zero posterior mass across all states")
    return g / z


def expected_mismatches(hmm: CohortHMM, haplotypes_chain: np.ndarray,
                        target_1d: np.ndarray, gamma: np.ndarray) -> float:
    """`E[# positions where the released allele differs from the target]`. EXACT.

    Positions where the target has no call are excluded, matching the denominator
    `target_fidelity` is quoted against. Zero Monte-Carlo error: this is a linear
    functional of the marginals, not an estimate from draws.
    """
    known = target_1d >= 0
    mism = (haplotypes_chain != target_1d[None, :]) & known[None, :]
    return float((gamma.T * mism).sum())


def _categorical(p: np.ndarray, rng: np.random.Generator) -> int:
    """One draw from a categorical, by inverse CDF.

    Written out rather than `rng.choice(K, p=p)` for two reasons: `choice`
    validates and copies on every call, which matters 305,886 times per release;
    and it demands the probabilities sum to 1 within a fixed tolerance, so a
    legitimate 1e-15 drift from the forward normalisation raises instead of
    sampling. Clipping the search index is the standard guard for a uniform draw
    landing beyond the last cumulative bin through rounding.
    """
    cdf = np.cumsum(p)
    return int(min(int(np.searchsorted(cdf, rng.random() * cdf[-1], side="right")),
                   len(p) - 1))


def _digest(a: np.ndarray) -> str:
    """Content digest of an array, so a record identifies the object it describes."""
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


@dataclass(frozen=True)
class Release:
    """One released path, self-describing. Field names are `docs/data.md` section 4.

    ⚠ `release_id` IS NOT A UNIQUE KEY on its own: the key is
    (`target_id`, `tau`, `utility_id`, `seed`), and a collision must raise rather
    than overwrite. ⚠ ONE RELEASE PER GENOME, cached and reused -- independent
    re-releases of the same path compose and degrade the guarantee, so nothing
    here may be called per query (standing constraint 4).
    """

    target_id: str
    tau: float
    eta_tau: float
    utility_id: str
    graph_id: str
    seed: int
    state_path: np.ndarray        # (T,) or (T, 2) int32, the sampled Z_1:T
    released_path: np.ndarray     # (T,) or (T, 2) int16, the alleles released
    log_z_p: float
    u_achieved: float
    u_ceiling: float
    ploidy: int = 1               # 2 = one diploid release, both haplotypes
    loo_id: str | None = None     # comparison arm only; None for a real release
    released_at: str = field(default_factory=lambda: date.today().isoformat())

    @property
    def release_id(self) -> str:
        return (f"{self.target_id}_tau{self.tau:g}_{self.utility_id}"
                f"_seed{self.seed}".replace("/", "-").replace("=", ""))

    def record(self) -> dict:
        """The quotable record. The released PATH is identified; the target is not.

        ⚠ The arrays are carried by DIGEST, not inlined. `docs/data.md` section 4
        lists `state_path` and `released_path` as fields, and the first version of
        this method omitted BOTH -- so the "quotable object" appeared in no record
        at all. Inlining 305,886 integers per haplotype into every summary row is
        not what that schema wants either: the digest lets a record answer "is this
        the same release?", while the arrays go to disk through `store.save()`,
        which attaches their own provenance. `docs/data.md` records the deviation.
        """
        return {
            "release_id": self.release_id, "target_id": self.target_id,
            "tau": self.tau, "eta_tau": self.eta_tau,
            "utility_id": self.utility_id, "graph_id": self.graph_id,
            "loo_id": self.loo_id, "seed": self.seed,
            "log_z_p": self.log_z_p, "u_achieved": self.u_achieved,
            "u_ceiling": self.u_ceiling, "ploidy": self.ploidy,
            "released_at": self.released_at,
            "state_path_sha256": _digest(self.state_path),
            "released_path_sha256": _digest(self.released_path),
            "T": int(self.state_path.shape[0]),
            "donor_switches": int((np.diff(self.state_path, axis=0) != 0).sum()),
        }


def _draw_one(hmm: CohortHMM, haplotypes_chain: np.ndarray, target_1d: np.ndarray,
              beta: np.ndarray, eta: float, seed, baseline_norm: float
              ) -> tuple[np.ndarray, np.ndarray, float, float]:
    """One tilted chain: `(state_path, released_alleles, log_z, u_achieved)`.

    `seed` is a `numpy.random.SeedSequence`, not an integer, so the caller owns
    stream separation -- see the note in `release()` about the overlap the
    integer form caused.

    Takes `eta` rather than `tau` because the diploid release runs this at
    `eta_tau / 2` per haplotype while the RECORD must still carry the real `tau`.
    Keeping the split here means neither caller can quietly stamp a release with
    the wrong privacy parameter.
    """
    psi = utility.tilt_potentials(hmm.allowed, haplotypes_chain, target_1d, beta, eta)
    fwd = forward(hmm, psi)
    lz = log_z_p(fwd, baseline_norm, eta)

    z = backward_sample(hmm, fwd, np.random.default_rng(seed))
    released = haplotypes_chain[z, np.arange(hmm.T)].astype(np.int16)
    if (released == -1).any():
        raise AssertionError(
            "the released path carries a missing allele, so RENORMALISE did not "
            "hold: the chain occupied a donor with no call")

    u = float(utility.u_of_output(released, target_1d, beta)[0])
    return z, released, lz, u


def release(hmm: CohortHMM, haplotypes_chain: np.ndarray, target_chain: np.ndarray,
            util: utility.Utility, tau: float, *, seed: int, target_id: str,
            graph_id: str = "", loo_id: str | None = None,
            baseline_norm: float | None = None) -> Release:
    """Draw one release for one target at one `tau`. The whole mechanism, once.

    `haplotypes_chain` and `target_chain` are over CHAIN positions, in chain
    order, matching `hmm`. `target_chain` is (T,) for a haploid release or (T, 2)
    for a diploid one -- a diploid target's two haplotypes are ONE private input
    and get ONE release, because that is what the guarantee covers and what
    "one release per genome" means.

    `baseline_norm` is `baseline_log_norm(hmm)`; pass it in when releasing several
    targets or several `tau` values against one cohort, since it is
    target-independent and costs a full forward pass to recompute.

    ⚠ `tau` = 0 IS NOT A SPECIAL CASE. `eta_tau` = 0 makes every tilt potential 1
    and the draw is a sample from `R_D` alone -- the floor column of every result
    table, and a real cohort haplotype mosaic carrying zero target information.
    It is reached by the identical code path, so the floor cannot drift away from
    the mechanism it is the floor of.
    """
    eta = calibration.eta_tau(tau)
    base = baseline_log_norm(hmm) if baseline_norm is None else baseline_norm
    target_chain = np.asarray(target_chain)

    # ⚠ BIND THE ALLELES TO THE MODEL. `CohortHMM` carries only the availability
    # mask, never the alleles, so this argument is the one place a caller could
    # hand the sampler a DIFFERENT cohort of the same shape -- a target-pruned
    # matrix, say, which is the standing-constraint-1 violation the directory split
    # was built to prevent and the one shape it cannot catch. The `-1` check after
    # the draw only catches the subset that happens to be occupied, so the binding
    # is asserted here instead.
    if haplotypes_chain.shape != (hmm.K, hmm.T):
        raise ValueError(
            f"haplotypes_chain is {haplotypes_chain.shape}, but the model is "
            f"(K, T) = {(hmm.K, hmm.T)}. These are different cohorts or a "
            f"different chain_span.")
    if not np.array_equal(haplotypes_chain != -1, hmm.allowed):
        bad = int((( haplotypes_chain != -1) != hmm.allowed).sum())
        raise ValueError(
            f"haplotypes_chain disagrees with the model's availability mask at "
            f"{bad} cells. The mask and the alleles describe different cohorts, so "
            f"R_D is not the distribution these alleles came from -- which is how a "
            f"target-pruned state space would enter (standing constraint 1).")

    if target_chain.ndim == 1:
        z, released, lz, u = _draw_one(hmm, haplotypes_chain, target_chain,
                                       util.beta, eta,
                                       np.random.SeedSequence(seed).spawn(1)[0], base)
        ceiling = utility.u_ceiling(target_chain, util.beta)
        ploidy = 1
    elif target_chain.ndim == 2 and target_chain.shape[1] == 2:
        # The diploid mean factorises exactly:
        #     exp(eta * (u1 + u2)/2) = exp((eta/2) u1) * exp((eta/2) u2)
        # so the two chains are INDEPENDENT at half the exponent, `log Z_p` is the
        # SUM of the two, and that sum lands in [0, eta/2] + [0, eta/2] = [0, eta]
        # -- the bracket holds for the diploid release as a whole, at the real tau.
        # ⚠ Running each haplotype at the full `eta` is the SUM aggregator in
        # disguise: it doubles the effective tilt, puts the utility range at 2,
        # and `tau` stops covering the release.
        # ⚠ SeedSequence.spawn, NOT `seed + h`. The first version derived the two
        # haplotype streams as `seed` and `seed + 1`, so release(seed=s) and
        # release(seed=s+1) SHARED a draw -- measured, haplotype 1 of one release
        # was bit-identical to haplotype 0 of the next across all 305,886
        # positions. The three "independent" draws in the first chr21 run were
        # therefore four correlated haplotype draws, which inflated the observed
        # seed-to-seed spread 2.2x. spawn() gives provably disjoint streams from
        # one recorded integer, so a release stays reproducible from its `seed`.
        streams = np.random.SeedSequence(seed).spawn(2)
        draws = [_draw_one(hmm, haplotypes_chain, target_chain[:, h], util.beta,
                           eta / 2.0, streams[h], base) for h in (0, 1)]
        z = np.stack([d[0] for d in draws], axis=1)
        released = np.stack([d[1] for d in draws], axis=1)
        lz = draws[0][2] + draws[1][2]
        # ⚠ BRACKET THE SUM. Each haplotype's log_z was checked against eta/2, but
        # nothing checked their sum until 2026-09-30 -- so restoring the full `eta`
        # per haplotype (the SUM aggregator in disguise) was ACCEPTED at every tau
        # and returned a record whose published log_z_p violated its own stated
        # invariant. Measured at tau = 0.5: 0.763 against a bracket of [0, 0.549].
        if not (-1e-6 <= lz <= eta + 1e-6):
            raise AssertionError(
                f"diploid log_z_p = {lz!r} is outside [0, eta_tau] = [0, {eta!r}] "
                f"(haplotype values {draws[0][2]!r} and {draws[1][3]!r}). The two "
                f"chains must each run at eta_tau/2 so their sum lands in the "
                f"bracket; running each at the full eta_tau doubles the effective "
                f"tilt and tau no longer covers the release.")
        u = utility.diploid_utility(draws[0][3], draws[1][3])
        ceiling = utility.diploid_utility(
            utility.u_ceiling(target_chain[:, 0], util.beta),
            utility.u_ceiling(target_chain[:, 1], util.beta))
        ploidy = 2
    else:
        raise ValueError(
            f"target_chain must be (T,) or (T, 2), got {target_chain.shape}")

    return Release(
        target_id=target_id, tau=float(tau), eta_tau=eta,
        utility_id=utility_id(util), graph_id=graph_id, seed=int(seed),
        state_path=z, released_path=released, log_z_p=lz,
        u_achieved=u, u_ceiling=ceiling, ploidy=ploidy, loo_id=loo_id)


def utility_id(util: utility.Utility) -> str:
    """The utility's identity: `phi_v`, the weights, and the `beta_v` scope.

    Stamped on every release because two runs under different scopes are
    different quantities wearing one name, and `docs/plan.md` R1 makes the scope
    an open decision -- so a number that does not carry it cannot be compared to
    anything later.
    """
    return f"indicator/w=inv_support/beta={util.scope}"
