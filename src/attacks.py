"""The attacker set, and the arm guard that stops an Arm 1 attack being read as Arm 2.

⚠⚠ WHY THE GUARD EXISTS, AND IT IS THE POINT OF THIS MODULE. The released object
contains NONE of the target's own alleles. Every released allele is some cohort
donor's, because the output support holds only cohort-supported paths; the only
target-dependent quantity in the whole object is the TILT. That is not the shape the
standard genomic attacks assume -- they assume the release contains the target's own
genotypes, possibly perturbed -- and two of them do not merely weaken against us,
they become ILL-POSED or INVERTED in Arm 2:

  * DATABASE SEARCH ("which database row produced this release?") has NO CORRECT
    ANSWER in Arm 2. The target is outside G, outside D, outside every panel and
    outside the attack database BY CONSTRUCTION. So a 0% link rate is the
    experimental design, not privacy, and reporting it as privacy would be a number
    that cannot come out any other way.
  * MEMBERSHIP INFERENCE INVERTS. The release is a mosaic of cohort donors, so it
    matches cohort members well -- which identifies DONORS, not the target. In Arm 2
    there is no positive class at all, and a high match fraction is evidence about
    the public cohort rather than about our guarantee. Read as membership it would
    reverse the conclusion.

An Arm 1 attacker carried into Arm 2 unexamined would therefore produce a number
that LOOKS like a privacy result. This module makes that combination raise instead of
returning a figure, which is the same device `src/cli.py` uses for standing
constraint 1: make the error inexpressible rather than something to remember not to do.

WHAT IS WELL POSED IN BOTH ARMS is the theorem's own game: an equal-prior binary
hypothesis test between two input paths. The proposal specifies no attacker -- it is
seven pages of theory with no evaluation section -- but its algebra hands over the
PROVABLY OPTIMAL one. With

    Q^tau_p(y) = R_D(y) exp(eta_tau u(p,y)) / Z_p

the Neyman-Pearson statistic between inputs p and q is

    Lambda(y) = eta_tau [u(p,y) - u(q,y)] - [log Z_p - log Z_q]

with `R_D` cancelling. By Neyman-Pearson its advantage EQUALS TV(Q_p, Q_q), the
quantity Theorem 1 bounds, so no attack in any literature can beat it for a given
pair. That is worth more than it sounds: in the DP-auditing literature the trade-off
function is unknown and must be approximated by heuristics, which is why "was my
attacker strong enough" is an open problem there. Here it is closed form, so this
attacker provably attains the optimum -- and that is the whole answer to a reviewer
asking whether a null came from a weak attack. ⚠ It is only available if the
PRE-REGISTERED attacker is this exact likelihood ratio and not a heuristic.

⚠ SO STRENGTH IS NOT THE RISK; POWER IS. Even the optimal attacker returns a null at
any feasible n -- R5 measures signal-to-noise about 1:900 and standing constraint 5
permits one release per genome. A null without a power statement is what standing
constraint 7 forbids calling evidence. Every attack here therefore carries a
POSITIVE CONTROL: a stated condition under which it must demonstrably SUCCEED.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

ARMS = ("arm1", "arm2")

# How an attack is posed, which is what decides whether an arm admits it.
KINDS = ("lrt_pair", "db_search", "membership", "reconstruction")

# Why each (kind, arm) pairing is refused. Keyed by kind; the text is the whole
# point of the refusal, so it lives here rather than in a bare exception.
_ILL_POSED_IN_ARM2 = {
    "db_search": (
        "a database search has NO CORRECT ANSWER in arm 2. The target is external "
        "by construction -- outside G, outside D, outside every panel and outside "
        "the attack database -- so no database row produced this release. A 0% link "
        "rate is the experimental design, not privacy, and it cannot come out any "
        "other way. Re-pose it as a 1-of-2 discrimination with the target supplied "
        "as one of the two candidates (kind='lrt_pair'), which IS well posed and is "
        "bounded by p_succ_bound."
    ),
    "membership": (
        "membership inference INVERTS in arm 2. There is no positive class: the "
        "target is not a member. And the release is a mosaic of cohort donors, so "
        "it matches cohort members well -- a high match fraction identifies DONORS, "
        "not the target, and is evidence about the public cohort rather than about "
        "our guarantee. Reported as membership it would reverse the conclusion."
    ),
}


def assert_well_posed(kind: str, arm: str) -> None:
    """Refuse an attack whose question has no answer in this arm.

    ⚠ This is a GUARD, not a preference. It exists because the refused combinations
    return a plausible number rather than an error, and that number reads as a
    privacy result to anyone who did not check which arm it came from.
    """
    if arm not in ARMS:
        raise ValueError(f"arm must be one of {ARMS}, not {arm!r}")
    if kind not in KINDS:
        raise ValueError(f"kind must be one of {KINDS}, not {kind!r}")
    if arm == "arm2" and kind in _ILL_POSED_IN_ARM2:
        raise ValueError(
            f"attack kind {kind!r} is ILL POSED in arm 2: "
            f"{_ILL_POSED_IN_ARM2[kind]}"
        )


@dataclass(frozen=True)
class AttackSpec:
    """One pre-registered attack. The arm is part of its identity, not a runtime flag.

    ⚠ `positive_control` is required and must be non-empty. An attack with no stated
    condition under which it must SUCCEED cannot have its failure reported as
    privacy (standing constraint 7), and the repository's forecast is that every
    attack here WILL fail at the operating configuration -- so the control is the
    load-bearing field, not the metric.
    """

    name: str
    arm: str
    kind: str
    statistic: str
    chance_level: float
    positive_control: str

    def __post_init__(self) -> None:
        assert_well_posed(self.kind, self.arm)
        if not self.positive_control.strip():
            raise ValueError(
                f"attack {self.name!r} has no positive control. An attacker that "
                f"fails is not evidence the bound holds (standing constraint 7), so "
                f"a stated condition under which it must SUCCEED is mandatory."
            )
        if not 0.0 < self.chance_level < 1.0:
            raise ValueError(f"chance_level must be in (0,1), got {self.chance_level}")


def lrt_statistic(u_p: np.ndarray, u_q: np.ndarray, eta: float,
                  log_z_p: float, log_z_q: float) -> np.ndarray:
    """`Lambda(y)` for every `y`: the Neyman-Pearson optimal statistic.

    `R_D(y)` cancels between the two hypotheses, which is why the baseline never
    appears here and why the statistic needs no knowledge of the prior beyond the
    two normalisers. Positive means "input p is likelier".
    """
    u_p = np.asarray(u_p, dtype=np.float64)
    u_q = np.asarray(u_q, dtype=np.float64)
    if u_p.shape != u_q.shape:
        raise ValueError(f"u_p is {u_p.shape}, u_q is {u_q.shape}")
    if eta < 0:
        raise ValueError(f"eta_tau must be non-negative, got {eta}")
    return eta * (u_p - u_q) - (log_z_p - log_z_q)


def exact_advantage(q_p: np.ndarray, q_q: np.ndarray) -> float:
    """The optimal attacker's accuracy at equal priors, EXACTLY: `(1 + TV)/2`.

    By Neyman-Pearson the best achievable accuracy in the two-hypothesis game is
    `(1 + TV(Q_p, Q_q)) / 2`, so this is not an estimate of an attacker's
    performance -- it is the ceiling no attacker can exceed, and the quantity
    Theorem 1 bounds by `(1 + tau) / 2`.

    ⚠ A MEASURED accuracy above this falsifies the IMPLEMENTATION, not the theorem.
    """
    q_p = np.asarray(q_p, dtype=np.float64)
    q_q = np.asarray(q_q, dtype=np.float64)
    if q_p.shape != q_q.shape:
        raise ValueError(f"supports differ: {q_p.shape} vs {q_q.shape}")
    for nm, q in (("q_p", q_p), ("q_q", q_q)):
        if abs(q.sum() - 1.0) > 1e-9:
            raise ValueError(f"{nm} sums to {q.sum()!r}, not 1")
    return 0.5 * (1.0 + 0.5 * float(np.abs(q_p - q_q).sum()))


def empirical_accuracy(stat_under_p: np.ndarray, stat_under_q: np.ndarray) -> float:
    """Accuracy of the sign rule on draws from each hypothesis, at equal priors.

    Ties are split, which is the right convention at equal priors and matters on a
    small support where `Lambda` is exactly 0 for outputs that both inputs explain
    identically -- at `tau` = 0 that is EVERY output, so this returns exactly 0.5
    there rather than a number that depends on a tie-breaking accident.
    """
    a = np.asarray(stat_under_p, dtype=np.float64)
    b = np.asarray(stat_under_q, dtype=np.float64)
    correct_p = (a > 0).sum() + 0.5 * (a == 0).sum()
    correct_q = (b < 0).sum() + 0.5 * (b == 0).sum()
    return float((correct_p / len(a) + correct_q / len(b)) / 2.0)
