"""`u_path`, and the two pieces it is built from. The only place the target enters.

Under model C the utility collapses to something very simple:

    u(p, y) = sum_v beta_v * phi_v(p_v, y_v),   phi_v = 1[y_v == p_v],
    beta_v  = w_v / sum_u w_u,                  w_v  = 1 / support(v)

`phi_v` is an INDICATOR -- there is no block to aggregate over, no per-block
normaliser and no aggregation choice to make. That is what choosing the variant
as the chain unit bought, and it is why model C is simpler than the block-based
design it replaced rather than harder.

TWO PROPERTIES ARE LOAD-BEARING AND ARE ASSERTED AT RUNTIME, NOT DOCUMENTED.

1. `sum_v beta_v == 1`. The normalisation is what keeps `u` inside [0,1]; weights
   that do not sum to 1 break the CALIBRATION, not merely the emphasis, because
   `eta_tau = arctanh(tau)` assumes a utility of range exactly 1.
2. `u` in [0,1]. An unbounded utility does not weaken the bound -- it removes it
   (standing constraint 2). A run-time assertion is the difference between a
   wrong number and a crash.

⚠ THE BUDGET DOES NOT THIN AS T GROWS, which is the intuition to discard here.
Since `sum beta = 1`, the total tilt between a perfectly-matching and a
non-matching path is `exp(eta_tau)` however the path is partitioned: each
position moves less and there are proportionally more of them. Measured at
tau = 0.5, the per-position tilt is 1.0000016 at T = 340,849 against 1.0000055 at
T = 100,757, and `exp(eta_tau) = 1.7321` either way. What DOES change with T is
numerical headroom -- see `ffbs.py`.

⚠ WHAT `beta_v` EMPHASISES IS INHERITED AND NOT ARGUED FOR (deferred decision R3).
`w_v = 1/support(v)` upweights POORLY-CALLED positions, which is not the same as
INFORMATIVE ones: a site where 87 of 88 donors are called but the target carries a
singleton is highly informative and gets low weight. It came from PanMixer, where
it was a cost weight doing a different job. It is survivable for now because on
this cohort it is nearly uniform -- median support 88, and 98.6% of chain
positions have support >= 80 and hold 96.8% of the weight -- so the 88x ratio
exists in principle and applies to almost nothing. Probably it should be rarity
(`-log f`) instead. That is a pre-registration item, not a code change to make
casually: standing constraint 5 says the utility function is frozen before any
attack outcome is inspected.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# `beta_v` normalisation scopes. The choice is deferred decision R1 and it must
# be stated on every number, because it moves the achievable ceiling.
BETA_SCOPES = ("chain", "chain_readable")


@dataclass(frozen=True)
class Utility:
    """`beta_v` over the chain, plus the provenance needed to quote a number.

    `scope` and `u_ceiling` are carried rather than recomputed because a utility
    figure is uninterpretable without them -- see the reporting contract. Two
    runs under different scopes are different quantities wearing one name.
    """

    beta: np.ndarray        # (T,) float64 over CHAIN positions, sums to 1
    scope: str
    n_scope: int            # positions the normaliser ran over
    n_chain: int            # chain positions in total

    @property
    def T(self) -> int:
        return len(self.beta)


def support_weights(support: np.ndarray) -> np.ndarray:
    """`w_v = 1 / support(v)`, refusing the zero it cannot handle.

    ⚠ `support(v) = 0` is EXCLUDED FROM THE CHAIN upstream rather than patched
    here. The obvious patch, `1/(support + 1)`, was rejected deliberately: a
    position where no donor has a called allele is unusable, not merely awkward
    -- the sampler has nothing to emit and RENORMALISE has an empty pool -- and
    the patch would keep it in the chain with a finite `beta_v` share it can
    never earn, re-creating the inert-budget problem one position at a time while
    hiding a position that cannot produce output. The exclusion is a function of
    the cohort alone, so it stays target-independent.
    """
    support = np.asarray(support)
    if (support < 0).any():
        raise ValueError("negative support; the matrix and the VCF disagree")
    zero = np.flatnonzero(support == 0)
    if len(zero):
        raise ValueError(
            f"{len(zero)} positions have support 0 (first at chain index "
            f"{zero[0]}). Those are excluded from the chain by "
            f"genmap.chain_mask, not weighted here. Rebuild in_chain.npy."
        )
    return 1.0 / support.astype(np.float64)


def beta_v(support: np.ndarray, *, in_chain: np.ndarray | None = None,
           readable: np.ndarray | None = None, scope: str = "chain") -> Utility:
    """`beta_v` over the chain positions, summing to 1 across the chosen scope.

    `support`, `in_chain` and `readable` are all over the FULL site axis; the
    returned `beta` is over chain positions, in chain order, so it lines up with
    `CohortHMM`'s arrays.

    ⚠ WHICH SCOPE `tau` COVERS IS DEFERRED DECISION R1, AND IT IS NOT COSMETIC.

    `scope="chain"` normalises over every chain position. Honest, and the current
    default. On chr21 that means 16.4% of the weight sits on positions where the
    target has no call at all, so `u` can never exceed 0.815: the full `tau` is
    paid for 81.5% of the leverage, and the release at the other positions is a
    pure prior draw with no reference to the target whatsoever.

    `scope="chain_readable"` normalises over `in_chain & readable`, taking the
    ceiling to 1.000 over 253,846 of 305,886 chain positions. Legitimate in the
    test harness because `readable` is derived from the target PANEL's site list,
    which is one file for all ten test targets -- verified 1 distinct set, and
    locked by a test.

    ⚠ AND IT DOES NOT SURVIVE CONTACT WITH A REAL TARGET AS WRITTEN (R4). A
    genuinely novel genome arrives from its own calling pipeline with its own
    site list, so `readable` would differ per target, and normalising over it
    would be target-dependent after all -- the exact standing-constraint-1
    violation the scope choice is meant to avoid. Before adopting this scope for
    anything but the harness, `readable` must be redefined as a function of
    something fixed and public: `G` itself, or a named reference callset pinned in
    advance. Never of the target's own VCF.

    ⚠ AND IT IS PER-CHROMOSOME EITHER WAY. Today `beta` sums to 1 over one
    chromosome, which silently means "`tau` per chromosome". A genome-wide `tau`
    needs a genome-wide normaliser, with each chromosome taking a SHARE of the
    budget; 22 chromosomes each normalised to 1 would be 22 independent releases
    of the same person and would compose, which standing constraint 4 forbids.
    Decide before any multi-chromosome run.
    """
    if scope not in BETA_SCOPES:
        raise ValueError(f"scope must be one of {BETA_SCOPES}, not {scope!r}")
    support = np.asarray(support)
    if in_chain is None:
        in_chain = np.ones(len(support), dtype=bool)
    if len(in_chain) != len(support):
        raise ValueError(
            f"in_chain has {len(in_chain)} entries, support has {len(support)}")

    w = support_weights(support[in_chain])

    if scope == "chain":
        keep = np.ones(len(w), dtype=bool)
    else:
        if readable is None:
            raise ValueError(
                "scope='chain_readable' needs the readable mask. It is emitted "
                "beside path.npy and is target-independent ONLY because our test "
                "targets share one panel -- see R4 before using it for real.")
        if len(readable) != len(support):
            raise ValueError(
                f"readable has {len(readable)} entries, support has {len(support)}")
        keep = readable[in_chain]

    if not keep.any():
        raise ValueError(f"scope {scope!r} selected no positions; beta is undefined")

    beta = np.where(keep, w, 0.0)
    total = beta.sum()
    if not np.isfinite(total) or total <= 0:
        raise AssertionError(f"beta weights sum to {total!r}; cannot normalise")
    beta = beta / total

    s = beta.sum()
    if abs(s - 1.0) > 1e-12:
        raise AssertionError(f"sum beta_v = {s!r}, not 1. The calibration is void.")
    if (beta < 0).any():
        raise AssertionError("negative beta_v")

    return Utility(beta=beta, scope=scope, n_scope=int(keep.sum()),
                   n_chain=int(in_chain.sum()))


def phi_v(donor_alleles: np.ndarray, target_allele: int) -> np.ndarray:
    """`phi_v` at one position: 1 where the donor carries the target's allele.

    ⚠ A TARGET WITH NO CALL SCORES ZERO FOR EVERY DONOR, which makes the tilt
    potential 1 everywhere and the position a pure prior draw. It is NOT a
    wildcard match -- rewarding donors for agreeing with an unknown value would
    invent information. It is also not the same thing as `absent_policy =
    "reference"`, which fills allele 0 and would reward donors for matching a
    FABRICATED call: absence from a callset means the caller reported no variant
    there, not that this person is reference.

    ⚠ A DONOR WITH NO CALL ALSO SCORES ZERO HERE, but that is not how it is
    handled -- `missing_policy` = RENORMALISE removes it from the state set
    entirely, via `CohortHMM.allowed`, so the chain cannot occupy it at all. The
    zero here is only so this function is total.
    """
    donor_alleles = np.asarray(donor_alleles)
    if target_allele < 0:
        return np.zeros(len(donor_alleles), dtype=np.float64)
    return ((donor_alleles == target_allele) & (donor_alleles >= 0)).astype(np.float64)


def u_of_output(output: np.ndarray, target: np.ndarray, beta: np.ndarray) -> np.ndarray:
    """`u(p, y)` for one or many released allele paths. Asserted in [0,1].

    `output` is (T,) or (n, T) of allele indices; `target` is (T,). The utility is
    a function of the RELEASED ALLELES, not of which donor they came from -- which
    is why the toy's D1/D3 collision is a real test and not decoration.
    """
    output = np.atleast_2d(np.asarray(output))
    target = np.asarray(target)
    if output.shape[1] != len(target) or len(beta) != len(target):
        raise ValueError(
            f"shape mismatch: output {output.shape}, target {target.shape}, "
            f"beta {beta.shape}")
    match = (output == target[None, :]) & (target[None, :] >= 0) & (output >= 0)
    u = match.astype(np.float64) @ beta
    _assert_bounded(u)
    return u


def u_ceiling(target: np.ndarray, beta: np.ndarray) -> float:
    """The largest `u` this target could possibly achieve: `sum beta_v` over its calls.

    ⚠ REPORT THIS WITH EVERY UTILITY FIGURE. It is the honest ceiling, and it is
    below 1 whenever the target has no call at some position carrying weight --
    which under `scope="chain"` is 16.4% of the chr21 weight. Quoting
    `utility_retained` against an implicit ceiling of 1 understates the mechanism
    by exactly that factor, and quoting it against this number without saying so
    overstates it. Both are the reporting error the two measurement paths in
    `docs/terms.md` are about.

    It is NOT the same as `target_fidelity`'s ceiling, which is about target
    variants with no record in the graph at all -- a different loss in the
    opposite direction. `docs/plan.md` item 1b keeps the three apart.
    """
    target = np.asarray(target)
    if len(target) != len(beta):
        raise ValueError(f"target has {len(target)} positions, beta {len(beta)}")
    return float(beta[target >= 0].sum())


def diploid_utility(u1: float, u2: float) -> float:
    """`(u1 + u2) / 2` for a diploid target, whose two haplotypes are ONE input.

    ⚠ THE MEAN, NEVER THE SUM. A sum leaves [0,1] and doubles the effective tilt,
    so the same `tau` no longer covers the release and the guarantee is gone. This
    is a one-character error with no visible symptom, which is why it has its own
    function and its own test.

    The pairing is positional: the target is phased, so haplotype 1 is scored
    against released haplotype 1. A maximum over the two pairings would be the
    more generous quantity and is deliberately not offered -- it does not
    decompose along the chain, so `ffbs` could not draw from the resulting tilted
    distribution exactly, and an approximate draw carries no bound at all.
    """
    u = 0.5 * (float(u1) + float(u2))
    _assert_bounded(np.array([u]))
    return u


def _assert_bounded(u: np.ndarray, *, tol: float = 1e-9) -> None:
    if u.size and (u.min() < -tol or u.max() > 1.0 + tol):
        raise AssertionError(
            f"u_path left [0,1]: [{u.min():.10g}, {u.max():.10g}]. The theorem "
            f"assumes range exactly 1; outside it, eta_tau = arctanh(tau) is the "
            f"wrong calibration and the bound does not hold (constraint 2)."
        )


def tilt_potentials(hmm_allowed: np.ndarray, haplotypes_chain: np.ndarray,
                    target_chain: np.ndarray, beta: np.ndarray, eta: float
                    ) -> np.ndarray:
    """(T, K) `psi_t(j) = exp(eta_tau * beta_v * phi_v)`, times the RENORMALISE mask.

    The ONLY place the private target enters the sampler; the cohort chain stays
    fixed. Returns the allowed-state indicator multiplied in, because the forward
    pass needs one potential array and keeping them separate invites applying one
    and forgetting the other.

    ⚠ Potentials cached across `tau` values are a correctness hazard -- `eta`
    multiplies every entry, so a stale array is a release at the wrong privacy
    level with no visible symptom. Nothing here caches.
    """
    K, T = haplotypes_chain.shape
    if hmm_allowed.shape != (K, T):
        raise ValueError(f"allowed is {hmm_allowed.shape}, expected {(K, T)}")
    for name, arr in (("target", target_chain), ("beta", beta)):
        if len(arr) != T:
            raise ValueError(f"{name} has {len(arr)} positions, chain has {T}")
    if eta < 0:
        raise ValueError(f"eta_tau must be non-negative, got {eta}")

    match = (haplotypes_chain == target_chain[None, :]) & (target_chain[None, :] >= 0)
    psi = np.exp(eta * (match * beta[None, :]))
    psi *= hmm_allowed
    return np.ascontiguousarray(psi.T)
