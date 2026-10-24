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
   preference: under model C the per-position tilt is eta_tau * mean(beta_v) =
   eta_tau / T. ⚠ **[CORRECTION - 2026-10-12]** This read "the CHAIN's T = 305,886
   ... 1.80e-06 ... only 15.1x the noise floor", which is the SUPERSEDED
   coordinate-cut chain. The live decision-K chain is T = 332,521, giving
   1.6519e-06 = 13.86x float32's 1.1921e-07. And the release runs each haplotype at
   `eta_tau / 2` (see `_draw_one`), so the REALISED per-position tilt is half that
   again: 8.2597e-07 = 6.93x. ⚠⚠ GENOME-WIDE IT GOES UNDER: at T = 22,134,562 the
   realised per-haplotype tilt at tau = 0.5 is 1.2408e-08, which is **0.104x float32
   epsilon** -- an order of magnitude BELOW the representable step, so in float32 the
   tilt would vanish entirely rather than merely lose precision. The guard below is
   load-bearing by a wider margin than it used to claim. (MODEL B's coarser
   T = 100,757 would have had 5.45e-06, 45.7x -- a fallback-model figure, quoted for
   contrast only, not a model-C number.) ⚠ The denominator is the chain, NOT the
   340,849-record site axis; quoting the site axis understates the headroom. Rounding
   would consume a meaningful fraction of the target's entire influence on its
   own release. Storage is about 215 MB at chr21's T = 305,886 and K = 88.

   ⚠ FORWARD-LOOKING, NOT IMPLEMENTED (the genome-wide work is tracked in
   `docs/plan.md`; nothing below describes a code path that exists). At genome
   scale float64 stops being a preference and becomes a CORRECTNESS requirement.
   `beta_v` sums to 1 genome-wide over the 22 autosomes (standing constraint 4),
   so T is the genome-wide chain length and the per-position tilt eta_tau / T
   falls BELOW float32's resolution once T exceeds eta_tau / 1.19e-07 -- about
   4.6e6 positions at tau = 0.5 and only 8.4e5 at tau = 0.10 (derived; the
   threshold tightens as tau falls). The 22-autosome chain length is UNMEASURED,
   but chr21 alone already contributes 305,886. Below that floor a dtype
   downgrade does not round the tilt, it ERASES it: the draw silently becomes an
   untilted sample from `R_D`, and `log_z_p` would sit at 0, inside its own
   bracket, so no assertion in this file would fire. The float64 check in
   `forward` is the only thing standing there.

3. ⚠ `A` IS NEVER MATERIALISED. `A_t = (1 - s_t) I + (s_t / K) J` is rank-one plus
   diagonal, so the forward step is one scalar sum and the backward draw is one
   vector update -- O(K) each. A dense matrix per position would be about
   21.1 GB on chr21 (340,849 * 88 * 88 * 8 B). ⚠ This read "21 TB" until
   2026-10-09 -- wrong by 1000x, in both this docstring and `hmm.py`.

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
        """`sum_t log c_t`. ⚠ NOT `log_z_p` -- see point 4 in the module docstring.

        ⚠ `math.fsum`, NOT `.sum()`, and not as a micro-optimisation. `log_z_p` is
        a DIFFERENCE of two of these whose magnitude grows with T while the result
        stays bounded by `eta_tau`, so the cancellation worsens with chromosome
        length. numpy's pairwise summation holds the error near 1.5e-10 at genome
        scale, but a sequential accumulation -- which is what someone reaches for
        when the `(T,)` array looks like wasted memory -- reaches ~9.9e-05 against
        a 1e-6 tolerance. `fsum` is exactly rounded, so the failure mode is closed
        rather than merely unlikely, and it costs nothing here. See `docs/bugs.md`.
        """
        return float(math.fsum(self.log_c))


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
            f"(the realised per-haplotype tilt is ~6.9x float32 resolution on chr21's "
            f"live chain and ~0.10x -- BELOW the representable step -- genome-wide)")
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


def expand_to_positions(hmm: CohortHMM, haplotypes_chain: np.ndarray,
                        z: np.ndarray) -> np.ndarray:
    """A sampled unit path -> the released allele at every CHAIN POSITION.

    Under model C a unit IS a position and this is a gather. Under model B the
    chosen donor is copied across the unit's whole span, which is what makes the
    unit atomic -- and is also why the release carries a donor's real haplotype
    verbatim over up to 618 chr21 variants. Departure (b), by construction.
    """
    if hmm.unit_spans is None:
        return haplotypes_chain[z, np.arange(hmm.T)].astype(np.int16)
    out = np.empty(hmm.n_positions, dtype=np.int16)
    for u, (a, b) in enumerate(hmm.unit_spans):
        out[a:b + 1] = haplotypes_chain[z[u], a:b + 1]
    return out


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
            *, budget: float | None = None, tol: float = 1e-6) -> float:
    """`log Z_p`, range-asserted against `[0, budget]`, which defaults to `eta_tau`.

    ⚠ `budget` EXISTS BECAUSE THE DEFAULT IS ONLY RIGHT FOR A WHOLE-BUDGET CHAIN.
    `Z_p = sum_y R_D(y) exp(eta * u)` reaches `e^eta` only if `u` can reach 1,
    i.e. only if this chain's `beta` sums to 1. Under decision B a chromosome's
    `beta` sums to its SHARE `w_c`, so its own `log Z` lives in
    `[0, eta * w_c]` -- roughly `eta/22` genome-wide. Bracketing that against the
    full `eta_tau` is a check that cannot fire, and the error it would wave
    through is running every chromosome at the full budget: the 22x composition
    violation of standing constraint 5. Callers pass
    `budget = exponent * beta.sum()`, which is tight by construction and cannot
    drift from the weights it is derived from.

    Theory pins it: `Z_p = sum_y R_D(y) exp(eta_tau u(p,y))` with `u` in [0,1], so
    `Z_p` in `[1, e^eta]`. If this leaves the bracket the sampler is wrong -- it is
    the cheapest correctness assertion in the whole implementation, so it is
    asserted and not reported.

    The tolerance covers accumulated floating-point error in the difference of two
    large sums, not modelling slack. On chr21 each sum is about -5.85e2 over
    305,886 terms -- MEASURED `baseline_log_norm` = -584.8135, the same figure
    point 4 of the module docstring quotes -- so the difference carries roughly
    5e-11 of error against a bracket of width 0.549 at `tau` = 0.5. (DERIVED, and
    unmeasured: the estimate scales with the magnitude of the sums, so correcting
    that magnitude from a stale -1.2e4 lowers the error figure by the same ~20x.)
    A failure at 1e-6 is a real defect.
    """
    cap = eta if budget is None else float(budget)
    if cap < -tol:
        raise ValueError(f"budget must be non-negative, got {cap!r}")
    lz = tilted.log_norm - baseline_log_norm_
    if not (-tol <= lz <= cap + tol):
        # Keep the whole-budget wording verbatim when no share is in play: it is
        # the form the existing regression test pins, and the only thing a share
        # changes is the cap, so say so additively rather than rewriting.
        where = ("[0, eta_tau]" if budget is None else
                 f"[0, eta_tau * share] = [0, {cap!r}] (eta_tau = {eta!r})")
        raise AssertionError(
            f"log_z_p = {lz!r} is outside {where}.\n"
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

    Under model B the marginals are per UNIT, so each donor's contribution is the
    number of mismatching positions across the whole unit -- the same linearity,
    just a coarser index.
    """
    known = target_1d >= 0
    mism = (haplotypes_chain != target_1d[None, :]) & known[None, :]
    if hmm.unit_spans is None:
        return float((gamma.T * mism).sum())
    cs = np.concatenate([np.zeros((mism.shape[0], 1)),
                         np.cumsum(mism, axis=1)], axis=1)
    per_unit = cs[:, hmm.unit_spans[:, 1] + 1] - cs[:, hmm.unit_spans[:, 0]]
    return float((gamma.T * per_unit).sum())


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
    here may be called per query (standing constraint 5).
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
    chrom: str = ""               # which chromosome this PART of the release covers
    share: float = 1.0            # its share of the one genome-wide budget
    released_at: str = field(default_factory=lambda: date.today().isoformat())

    @property
    def release_id(self) -> str:
        """⚠ THE CHROMOSOME IS PART OF THE KEY, because without it all 22 parts of
        one genome-wide release collapse to a single id -- and `docs/data.md` says
        a collision must raise rather than overwrite, so 21 of the 22 would be
        silently lost. Omitted when empty so every single-chromosome id on record
        stays byte-identical.
        """
        stem = (f"{self.target_id}_tau{self.tau:g}_{self.utility_id}"
                f"_seed{self.seed}")
        if self.chrom:
            stem += f"_{self.chrom}"
        return stem.replace("/", "-").replace("=", "")

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
            "chrom": self.chrom, "share": self.share,
            "released_at": self.released_at,
            "state_path_sha256": _digest(self.state_path),
            "released_path_sha256": _digest(self.released_path),
            "T": int(self.state_path.shape[0]),
            "donor_switches": int((np.diff(self.state_path, axis=0) != 0).sum()),
        }


def _draw_one(hmm: CohortHMM, haplotypes_chain: np.ndarray, target_1d: np.ndarray,
              beta: np.ndarray, eta: float, seed, baseline_norm: float,
              exact_out: dict | None = None
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
    psi = utility.tilt_potentials(hmm.allowed, haplotypes_chain, target_1d, beta, eta,
                                  unit_spans=hmm.unit_spans)
    fwd = forward(hmm, psi)
    # ⚠ The budget is this chain's OWN reach: `eta * sum(beta)`. With a
    # genome-wide normaliser `sum(beta)` is the chromosome's share, so this
    # is ~22x tighter than `eta` and is the only form that can fire.
    lz = log_z_p(fwd, baseline_norm, eta, budget=eta * float(beta.sum()))

    # ⚠⚠ THE EXACT ENDPOINT, FROM THE SAME `psi` AND THE SAME `fwd`. `E[differ]` and
    # `E_Q[u]` need the state marginals, which are a backward pass over the forward
    # messages this draw already computed. Computing them in a SEPARATE run would
    # recompute `psi` and `forward` from scratch: measured 2026-10-12, that doubling
    # was the difference between ~17.4 h and ~10.9 h for arm 1's 308 cells.
    # ⚠ It consumes NO randomness and is placed BEFORE the draw, so the released path
    # is bit-identical whether or not it was requested -- pinned by a test.
    if exact_out is not None:
        gam = state_marginals(hmm, fwd, psi)
        known = target_1d >= 0
        match = (haplotypes_chain == target_1d[None, :]) & known[None, :]
        exact_out["e_differ"] = exact_out.get("e_differ", 0.0) + expected_mismatches(
            hmm, haplotypes_chain, target_1d, gam)
        exact_out.setdefault("u_haps", []).append(
            float((gam * match.T).sum(axis=1) @ beta))
        exact_out["called"] = exact_out.get("called", 0) + int(known.sum())
        del gam

    z = backward_sample(hmm, fwd, np.random.default_rng(seed))
    released = expand_to_positions(hmm, haplotypes_chain, z)
    if (released == -1).any():
        raise AssertionError(
            "the released path carries a missing allele, so RENORMALISE did not "
            "hold: the chain occupied a donor with no call")

    u = float(utility.u_of_output(released, target_1d, beta)[0])
    return z, released, lz, u


def _streams(seed: int, stream_key: str, n: int):
    """`n` disjoint streams from one RECORDED integer, keyed by `stream_key`.

    ⚠⚠ WHY `stream_key` EXISTS (added 2026-10-12, before any genome release was written).
    `GenomeRelease` requires all 22 parts to agree on `seed`, because one genome release
    is one draw of one person. But the stream was derived from `seed` ALONE, so every
    chromosome of a genome release got the SAME two streams -- verified bit-identical --
    and `_pick` samples by inverse CDF on one uniform per position (see `_pick`). Every
    chromosome would therefore have chosen its donor at the SAME QUANTILE of its own
    posterior at every step: a systematic cross-chromosome coupling that nothing in the
    record would show. Keying the stream by the chromosome removes it while leaving the
    recorded `seed` the genome's, so a part stays reproducible from `(seed, chrom)` --
    both already fields of `Release`, so no record schema changed.

    ⚠ This is the FOURTH common-random-numbers defect in this project, and the first
    caught before it reached an artifact: no genome-scope release existed when it was
    found, so nothing on record is affected.

    ⚠ `stream_key=""` reproduces the previous behaviour EXACTLY, so every
    single-chromosome release already on disk stays reproducible from its seed.

    ⚠ hashlib, never `hash()`: Python's `hash()` is PYTHONHASHSEED-salted, which made a
    decoy choice irreproducible once already in this project.
    """
    if not stream_key:
        return np.random.SeedSequence(seed).spawn(n)
    h = int.from_bytes(hashlib.sha256(stream_key.encode()).digest()[:4], "big")
    return np.random.SeedSequence([int(seed), h]).spawn(n)


def release(hmm: CohortHMM, haplotypes_chain: np.ndarray, target_chain: np.ndarray,
            util: utility.Utility, tau: float | None = None, *,
            eta: float | None = None, seed: int, target_id: str,
            graph_id: str = "", loo_id: str | None = None,
            baseline_norm: float | None = None, stream_key: str = "",
            exact_out: dict | None = None) -> Release:
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
    # ⚠⚠ DRIVEN BY `tau` OR BY `eta`, NEVER BOTH. The `eta` entry point exists
    # because `tau` cannot ADDRESS the high end of the tilt: `tanh(eta)` is exactly
    # 1.0 in float64 past `eta` ~ 19 and `calibration.eta_tau` rightly refuses
    # `tau` = 1, so an eta-parameterised sweep (scripts/eta_sweep.py, and the
    # eta-driven attack cells) could not call this function at all. It is ONE
    # implementation with two entry points, never a second sampler.
    # ⚠ THE RECORD'S `tau` IS THEN `tanh(eta)`, WHICH SATURATES TO 1.0. That is the
    # true equivalent tau, and a release record carrying `tau` = 1 means the run was
    # eta-driven and its guarantee is VACUOUS -- read `eta_tau` instead, which is
    # exact at any magnitude. Nothing in this project persists an eta-driven release;
    # they exist for characterisation and for attack trials held in memory.
    if (tau is None) == (eta is None):
        raise ValueError("pass exactly one of tau= or eta=")
    if eta is None:
        eta = calibration.eta_tau(tau)
    else:
        if eta < 0.0:
            raise ValueError(f"eta must be non-negative, got {eta}")
        tau = math.tanh(eta)
    base = baseline_log_norm(hmm) if baseline_norm is None else baseline_norm
    target_chain = np.asarray(target_chain)

    # ⚠ BIND THE ALLELES TO THE MODEL. `CohortHMM` carries only the availability
    # mask, never the alleles, so this argument is the one place a caller could
    # hand the sampler a DIFFERENT cohort of the same shape -- a target-pruned
    # matrix, say, which is the standing-constraint-1 violation the directory split
    # was built to prevent and the one shape it cannot catch. The `-1` check after
    # the draw only catches the subset that happens to be occupied, so the binding
    # is asserted here instead.
    if haplotypes_chain.shape != (hmm.K, hmm.n_positions):
        raise ValueError(
            f"haplotypes_chain is {haplotypes_chain.shape}, but the model covers "
            f"(K, positions) = {(hmm.K, hmm.n_positions)} "
            f"({hmm.T} chain units under model {hmm.model}). These are different "
            f"cohorts or a different chain_span.")
    # ⚠ The mask is rebuilt from the alleles the SAME way `from_arrays` built it --
    # per position under model C, and ANDed across each unit under model B -- so the
    # guard stays exact in both and cannot be satisfied by a reshaped near-miss.
    called = haplotypes_chain != -1
    if hmm.unit_spans is None:
        rebuilt = called
    else:
        rebuilt = np.stack([called[:, a:b + 1].all(axis=1)
                            for a, b in hmm.unit_spans], axis=1)
    if not np.array_equal(rebuilt, hmm.allowed):
        bad = int((rebuilt != hmm.allowed).sum())
        raise ValueError(
            f"haplotypes_chain disagrees with the model's availability mask at "
            f"{bad} entries. The mask and the alleles describe different cohorts, so "
            f"R_D is not the distribution these alleles came from -- which is how a "
            f"target-pruned state space would enter (standing constraint 1).")

    if target_chain.ndim == 1:
        z, released, lz, u = _draw_one(hmm, haplotypes_chain, target_chain,
                                       util.beta, eta,
                                       _streams(seed, stream_key, 1)[0], base,
                                       exact_out=exact_out)
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
        streams = _streams(seed, stream_key, 2)
        draws = [_draw_one(hmm, haplotypes_chain, target_chain[:, h], util.beta,
                           eta / 2.0, streams[h], base, exact_out=exact_out)
                 for h in (0, 1)]
        z = np.stack([d[0] for d in draws], axis=1)
        released = np.stack([d[1] for d in draws], axis=1)
        lz = draws[0][2] + draws[1][2]
        # ⚠ BRACKET THE SUM. Each haplotype's log_z was checked against eta/2, but
        # nothing checked their sum until 2026-09-30 -- so restoring the full `eta`
        # per haplotype (the SUM aggregator in disguise) was ACCEPTED at every tau
        # and returned a record whose published log_z_p violated its own stated
        # invariant. Measured at tau = 0.5: 0.763 against a bracket of [0, 0.549].
        # ⚠ AND BRACKET IT AGAINST THIS CHROMOSOME'S SHARE, not the whole budget.
        # Each haplotype reaches eta/2 * share, so the diploid sum reaches
        # eta * share -- which equals eta only in the single-chromosome case where
        # share is 1. Under decision B a per-chromosome sum checked against the
        # full eta_tau is ~22x too loose and would accept every chromosome running
        # at the whole budget: standing constraint 5's composition violation, and
        # the error `test_diploid_gate.py` exists to catch.
        cap = eta * float(util.share)
        if not (-1e-6 <= lz <= cap + 1e-6):
            raise AssertionError(
                f"diploid log_z_p = {lz!r} is outside [0, eta_tau * share] = "
                f"[0, {cap!r}] with eta_tau = {eta!r} and share = {util.share!r} "
                f"(haplotype values {draws[0][2]!r} and {draws[1][2]!r}). The two "
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
        chrom=util.chrom, share=float(util.share),
        state_path=z, released_path=released, log_z_p=lz,
        u_achieved=u, u_ceiling=ceiling, ploidy=ploidy, loo_id=loo_id)


def utility_id(util: utility.Utility) -> str:
    """The utility's identity: `phi_v`, the weights, and the `beta_v` scope.

    Stamped on every release because two runs under different scopes are
    different quantities wearing one name, and `docs/plan.md` R1 makes the scope
    an open decision -- so a number that does not carry it cannot be compared to
    anything later.
    """
    base = f"indicator/w=inv_support/beta={util.scope}"
    # ⚠ The share belongs in the ID. A chromosome released at share 1 and the same
    # chromosome released at its genome-wide share are DIFFERENT utility functions
    # -- different weights, different ceiling, different reachable log_z -- and
    # standing constraint 8 forbids comparing two utility functions. Omitted when
    # it is 1 so every single-chromosome ID on record stays byte-identical.
    if abs(util.share - 1.0) > 1e-12:
        base += f"/share={util.share:.9g}"
    if util.chrom:
        base += f"/chrom={util.chrom}"
    return base


@dataclass(frozen=True)
class GenomeRelease:
    """One release of one genome, assembled from its per-chromosome parts.

    ⚠ WHY THIS TYPE EXISTS AT ALL. Under decision B there is ONE `tau` and ONE
    budget per person, spread across the 22 autosomes. Every invariant that makes
    that a release rather than 22 releases is a statement about the WHOLE set:
    `sum_c w_c == 1`, `sum_c log Z_c` in `[0, eta_tau]`, `sum_c u_c` in `[0,1]`,
    and one `(target, tau, utility, seed)` across all of them. A per-chromosome
    `Release` is individually well formed in all 22 cases where the budget has
    been spent 22 times over, so nothing that sees one chromosome can detect it.
    This is the object that can.

    ⚠ It also makes the released OBJECT match the guaranteed one. `tau` bounds the
    distinguishability of the whole released path; quoting a per-chromosome
    `log_z_p` or `u_achieved` as "the release's" is a category error that the
    reporting contract's arm/`tau`/denominator stamp would not catch.
    """

    parts: tuple[Release, ...]

    def __post_init__(self) -> None:
        if not self.parts:
            raise ValueError("a genome release needs at least one chromosome")
        first = self.parts[0]

        # One private input, one privacy level, one utility family, one draw.
        for field_name in ("target_id", "tau", "eta_tau", "seed", "ploidy",
                           "loo_id", "graph_id"):
            vals = {getattr(p, field_name) for p in self.parts}
            if len(vals) != 1:
                raise AssertionError(
                    f"the parts of one genome release disagree on {field_name!r}: "
                    f"{sorted(vals, key=repr)}. These are not parts of one "
                    f"release -- they are separate releases, and separate "
                    f"releases of one genome compose (standing constraint 5)."
                )

        # ⚠ The utility must be the same FAMILY but NOT the same id: each
        # chromosome's id carries its own share and name, so identical ids across
        # parts would mean the share was never applied.
        stems = {p.utility_id.split("/share=")[0].split("/chrom=")[0]
                 for p in self.parts}
        if len(stems) != 1:
            raise AssertionError(
                f"the parts use different utility functions: {sorted(stems)}. "
                f"Utility is comparable only at a matched utility id "
                f"(standing constraint 8).")

        chroms = [p.chrom for p in self.parts]
        if len(self.parts) > 1:
            if any(not c for c in chroms):
                raise AssertionError(
                    "a multi-chromosome release has a part with no `chrom`, so "
                    "the duplicate check below cannot do its job")
            if len(set(chroms)) != len(chroms):
                dup = sorted({c for c in chroms if chroms.count(c) > 1})
                raise AssertionError(
                    f"chromosome(s) {dup} appear more than once. Their share "
                    f"would be double counted, and the total could still look "
                    f"right with another chromosome missing.")

        total_share = math.fsum(p.share for p in self.parts)
        if abs(total_share - 1.0) > 1e-9:
            raise AssertionError(
                f"the parts' shares sum to {total_share!r}, not 1, over "
                f"{len(self.parts)} chromosomes: "
                f"{[(p.chrom or '?', p.share) for p in self.parts]}.\n"
                f"  A total of {len(self.parts)} means every chromosome was "
                f"normalised to the whole budget, which is {len(self.parts)} "
                f"independent releases of this person rather than one."
            )

        if not (-1e-6 <= self.log_z_p <= first.eta_tau + 1e-6):
            raise AssertionError(
                f"genome-wide log_z_p = {self.log_z_p!r} is outside "
                f"[0, eta_tau] = [0, {first.eta_tau!r}]. Per-chromosome values: "
                f"{[(p.chrom or '?', p.log_z_p) for p in self.parts]}."
            )
        if not (-1e-9 <= self.u_achieved <= 1.0 + 1e-9):
            raise AssertionError(
                f"genome-wide u_achieved = {self.u_achieved!r} is outside [0,1]; "
                f"an unnormalised utility removes the bound (standing "
                f"constraint 2)."
            )

    # `u` is the SUM over chromosomes, not the mean: each part's `u_c` is already
    # scaled by its own share, so the parts add to at most 1. Averaging would
    # divide the budget a second time and silently understate every release.
    @property
    def u_achieved(self) -> float:
        return float(math.fsum(p.u_achieved for p in self.parts))

    @property
    def u_ceiling(self) -> float:
        return float(math.fsum(p.u_ceiling for p in self.parts))

    @property
    def log_z_p(self) -> float:
        return float(math.fsum(p.log_z_p for p in self.parts))

    @property
    def donor_switches(self) -> int:
        return int(sum(int((np.diff(p.state_path, axis=0) != 0).sum())
                       for p in self.parts))

    @property
    def release_id(self) -> str:
        """No chromosome in it: this IS the release, and its parts are parts."""
        p = self.parts[0]
        return (f"{p.target_id}_tau{p.tau:g}_"
                f"{p.utility_id.split('/share=')[0].split('/chrom=')[0]}"
                f"_seed{p.seed}_genome").replace("/", "-").replace("=", "")

    def record(self) -> dict:
        p = self.parts[0]
        return {
            "release_id": self.release_id, "target_id": p.target_id,
            "tau": p.tau, "eta_tau": p.eta_tau, "seed": p.seed,
            "ploidy": p.ploidy, "loo_id": p.loo_id, "graph_id": p.graph_id,
            "scope": "genome", "n_chromosomes": len(self.parts),
            "chromosomes": [pp.chrom for pp in self.parts],
            "shares": [pp.share for pp in self.parts],
            "log_z_p": self.log_z_p, "u_achieved": self.u_achieved,
            "u_ceiling": self.u_ceiling,
            "donor_switches": self.donor_switches,
            "T": int(sum(pp.state_path.shape[0] for pp in self.parts)),
            "parts": [pp.release_id for pp in self.parts],
            "released_at": p.released_at,
        }
