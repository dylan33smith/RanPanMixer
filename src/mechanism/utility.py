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
casually: standing constraint 6 says the utility function is frozen before any
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

    beta: np.ndarray        # (T,) float64 over CHAIN positions, sums to `share`
    scope: str
    n_scope: int            # positions the normaliser ran over
    n_chain: int            # chain positions in total
    share: float = 1.0      # this chromosome's share of the ONE genome-wide budget
    weights: str = "inv_support"   # which w_v; see beta_v. "uniform" is the default
                                   # for the released tool as of 2026-10-24
    chrom: str = ""         # which chromosome, for the genome-level completeness check

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
           readable: np.ndarray | None = None, scope: str = "chain",
           share: float = 1.0, chrom: str = "",
           weights: str = "inv_support") -> Utility:
    """`beta_v` over the chain positions, summing to `share` across the chosen scope.

    ⚠ `share` IS THIS CHROMOSOME'S PORTION OF ONE GENOME-WIDE BUDGET (decision B,
    standing constraint 4), and the default of 1.0 is the single-chromosome case
    -- the whole budget on one chromosome. For a genome-wide release the caller
    passes `share = w_c` with `sum_c w_c == 1`, and `assert_shares_complete`
    checks that sum where it can be seen. Passing 1.0 per chromosome across 22
    chromosomes is 22 independent releases of one person; nothing inside this
    function can detect that, which is why the completeness check is separate.

    `support`, `in_chain` and `readable` are all over the FULL site axis; the
    returned `beta` is over chain positions, in chain order, so it lines up with
    `CohortHMM`'s arrays.

    ⚠ WHICH SCOPE `tau` COVERS IS DEFERRED DECISION R1, AND IT IS NOT COSMETIC.

    `scope="chain"` normalises over every chain position. Honest, and the current
    default. ⚠ **[CORRECTION - 2026-10-12]** This read "16.4% of the weight ... `u`
    can never exceed 0.815 ... paid for 81.5% of the leverage", which was measured on
    the SUPERSEDED coordinate-cut chain of T = 305,886. On the live decision-K chain
    of T = 332,521 it is **30.749104%** of the weight on positions where HG00096 has
    no call at all, so `u` cannot exceed **0.692509** — the inert share very nearly
    DOUBLED, because the rebuild added 26,635 positions of which only 4,362 are
    called. So the full `tau` is paid for 69.3% of the leverage, and the release at
    the other positions is a pure prior draw with no reference to the target
    whatsoever. ⚠ The figure is per target; this one is HG00096, ARM 2.

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

    ✅ **RESOLVED BY DECISION J (2026-10-11); THIS USED TO SAY "DECIDE BEFORE ANY
    MULTI-CHROMOSOME RUN".** `beta` now sums to the declared `share`, not to 1, and
    the assertion below checks it against `share` rather than against 1. Each
    chromosome takes `share` = T_c / T_total by variant count, the 22 autosome
    shares sum to 1, and `tau` is therefore genome-wide. Run per chromosome with
    `share` = 1 it still means "`tau` on that chromosome alone" -- which is what
    every chr21 figure on record is, legitimately, as a single-chromosome release.
    22 chromosomes each normalised to 1 would be 22 independent releases of the
    same person and would compose, which standing constraint 5 forbids; passing a
    share is what prevents that, and the caller must pass one.
    """
    if scope not in BETA_SCOPES:
        raise ValueError(f"scope must be one of {BETA_SCOPES}, not {scope!r}")
    if not (0.0 < share <= 1.0):
        raise ValueError(
            f"share must be in (0, 1], got {share!r}. A share above 1 is more "
            f"than the whole budget; a share of 0 is a chromosome with no weight, "
            f"which should be left out of the release rather than given zero.")
    support = np.asarray(support)
    if in_chain is None:
        in_chain = np.ones(len(support), dtype=bool)
    if len(in_chain) != len(support):
        raise ValueError(
            f"in_chain has {len(in_chain)} entries, support has {len(support)}")

    # ⚠⚠ `weights="uniform"` IS THE RELEASED TOOL'S CHOICE AS OF 2026-10-24, AND IT
    # IS ARGUED FOR RATHER THAN INHERITED. The full argument is in docs/terms.md
    # under `beta_v`; the four load-bearing parts:
    #  (1) It gives `u` a one-sentence meaning -- the FRACTION of the target's
    #      called positions the release reproduces -- which is comparable across
    #      targets, cohorts and tools. Under 1/support(v), `u` is a weighted
    #      average whose weights encode OUR cohort's genotyping completeness and
    #      which nobody can state the meaning of in a sentence.
    #  (2) 1/support(v) weights a property of OUR DATA, not of the position. It
    #      directs the budget at the positions where the cohort knows LEAST --
    #      measured genome-wide: 75,998 positions (0.343% of the chain) with
    #      support < 40 carry 2.491% of the weight, 7.3x their positional share,
    #      and a support-1 position carries 85.5x the mean.
    #  (3) It matches how a release is CONSUMED. A path is read positionally; no
    #      downstream tool reweights by our call rate.
    #  (4) It is the NEUTRAL INSTRUMENT for a claim about the mechanism.
    #      Concentration is the lever that moves the exchange rate (as sqrt(Var)),
    #      so a non-uniform default would conflate "the tilt does X" with "we
    #      chose a narrow utility". Uniform is the control; the weighting
    #      comparison reports the headroom separately.
    # ⚠ AND IT COSTS SOMETHING, which is the right direction to err: uniform is the
    # LEAST concentrated weighting, sum(beta^2) 4.5178e-08 against inv_support's
    # 7.3440e-08, so the headline result is CONSERVATIVE by about sqrt(1.63) = 1.27x.
    if weights == "uniform":
        w = np.ones(int(np.asarray(in_chain).sum()), dtype=np.float64)
    elif weights == "inv_support":
        w = support_weights(support[in_chain])
    else:
        raise ValueError(
            f"unknown weights {weights!r}; expected 'uniform' or 'inv_support'")

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

    beta = beta * share

    s = beta.sum()
    if abs(s - share) > 1e-12:
        raise AssertionError(
            f"sum beta_v = {s!r}, not the declared share {share!r}. The "
            f"calibration is void.")
    if (beta < 0).any():
        raise AssertionError("negative beta_v")

    return Utility(beta=beta, scope=scope, n_scope=int(keep.sum()),
                   n_chain=int(in_chain.sum()), share=float(share),
                   chrom=str(chrom), weights=str(weights))


def assert_shares_complete(utils: "list[Utility]", *, tol: float = 1e-12) -> float:
    """The genome-level half of standing constraint 4: `sum_c w_c == 1`.

    ⚠ THIS CHECK CANNOT LIVE IN `beta_v`, and that is the whole reason it is a
    separate function. Each chromosome's `Utility` is well formed in isolation --
    its weights sum to its own declared share, exactly as asserted -- and a set of
    22 of them summing to 22 is also made of 22 individually valid objects. Only
    something that sees every chromosome at once can tell one release of a genome
    from 22 independent releases of the same person, and the latter compose and
    degrade the guarantee (standing constraint 5).

    Returns the total, so a caller can stamp it.
    """
    if not utils:
        raise ValueError("no utilities to check; a genome-wide release needs at least one")
    seen: dict[str, int] = {}
    for u in utils:
        if u.chrom in seen:
            raise AssertionError(
                f"chromosome {u.chrom!r} appears twice in the budget. Its share "
                f"would be counted twice and the genome total would still look "
                f"right if another chromosome were missing.")
        seen[u.chrom] = 1
    if len(utils) > 1 and any(not u.chrom for u in utils):
        raise AssertionError(
            "a multi-chromosome budget has a Utility with no `chrom`, so the "
            "duplicate check above cannot do its job. Name every chromosome.")
    total = float(sum(u.share for u in utils))
    if abs(total - 1.0) > tol:
        raise AssertionError(
            f"the per-chromosome shares sum to {total!r}, not 1, over "
            f"{len(utils)} chromosomes: "
            f"{[(u.chrom or '?', u.share) for u in utils]}.\n"
            f"  `beta_v` sums to 1 over the whole scope `tau` covers -- the 22 "
            f"autosomes TOGETHER -- and each chromosome takes a SHARE. A total of "
            f"{len(utils)} means each chromosome was normalised to 1 "
            f"independently, which is {len(utils)} releases of one genome."
        )
    return total


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
    which under `scope="chain"` is **30.749104%** of the live chr21 weight for
    HG00096 (this read 16.4%, measured on the superseded T = 305,886 chain; the live
    chain is T = 332,521 -- corrected 2026-10-12). Quoting
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
                    target_chain: np.ndarray, beta: np.ndarray, eta: float,
                    *, unit_spans: np.ndarray | None = None) -> np.ndarray:
    """(T, K) `psi_t(j) = exp(eta_tau * beta_v * phi_v)`, times the RENORMALISE mask.

    The ONLY place the private target enters the sampler; the cohort chain stays
    fixed. Returns the allowed-state indicator multiplied in, because the forward
    pass needs one potential array and keeping them separate invites applying one
    and forgetting the other.

    ⚠ Potentials cached across `tau` values are a correctness hazard -- `eta`
    multiplies every entry, so a stale array is a release at the wrong privacy
    level with no visible symptom. Nothing here caches.
    """
    K, n_pos = haplotypes_chain.shape
    for name, arr in (("target", target_chain), ("beta", beta)):
        if len(arr) != n_pos:
            raise ValueError(f"{name} has {len(arr)} positions, chain has {n_pos}")
    if eta < 0:
        raise ValueError(f"eta_tau must be non-negative, got {eta}")
    # ⚠⚠ THE ONLY STRUCTURAL CHECK ON `beta` ON THIS PATH, added 2026-10-12. `beta_v`
    # asserts sum == share, but NOTHING REQUIRES A CALLER TO GO THROUGH IT, and
    # `ffbs._draw_one` hands `log_z_p` a bracket of `eta * float(beta.sum())` -- which
    # WIDENS in lockstep with a mis-normalised `beta` instead of firing on it. So before
    # this, a hand-built `beta` summing to 2 ran the whole release path with no symptom.
    # Every legitimate `beta` sums to one chromosome's SHARE of a single genome-wide
    # budget and the shares sum to 1, so `sum(beta) <= 1` ALWAYS. Above 1 means
    # `delta_u > 1`, and standing constraint 2 is explicit that an unnormalised utility
    # does not weaken the bound -- it REMOVES it. Checked here because this is the one
    # function every tilted release must call.
    s_beta = float(beta.sum())
    if s_beta > 1.0 + 1e-9:
        raise ValueError(
            f"sum beta_v = {s_beta!r}, above 1, so delta_u > 1 and "
            f"eta_tau = arctanh(tau) is the WRONG calibration: the release would carry "
            f"a larger tau than the one it is stamped with. The calibration is void.")

    match = (haplotypes_chain == target_chain[None, :]) & (target_chain[None, :] >= 0)
    weighted = match * beta[None, :]                       # (K, n_pos)

    if unit_spans is None:
        if hmm_allowed.shape != (K, n_pos):
            raise ValueError(f"allowed is {hmm_allowed.shape}, expected {(K, n_pos)}")
        exponent = weighted
    else:
        # MODEL B: one donor covers a whole unit, so the unit's exponent is the SUM
        # of its positions' beta * phi. `sum beta` over units is still 1, because it
        # is the same weights regrouped -- the calibration is untouched.
        T = len(unit_spans)
        if hmm_allowed.shape != (K, T):
            raise ValueError(f"allowed is {hmm_allowed.shape}, expected {(K, T)}")
        cs = np.concatenate([np.zeros((K, 1)), np.cumsum(weighted, axis=1)], axis=1)
        exponent = cs[:, unit_spans[:, 1] + 1] - cs[:, unit_spans[:, 0]]

    psi = np.exp(eta * exponent)
    psi *= hmm_allowed
    return np.ascontiguousarray(psi.T)
