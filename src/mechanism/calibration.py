"""tau and the three constants derived from it. One file, because they must agree.

`eta_tau = arctanh(tau)` is only correct when the utility is bounded in exactly
[0,1]. If a raw utility spans `delta_u != 1` the calibration is
`arctanh(tau)/delta_u`, and using `arctanh(tau)` regardless does not weaken the
guarantee -- it removes it. That is standing constraint 2 and the highest-severity
error class in the project, so the range is an explicit argument here rather than
an assumption made somewhere else.

⚠ CORRECTED 2026-09-30. The first version of this docstring said "the range is a
REQUIRED argument" -- it is an optional keyword defaulting to 1.0 -- and claimed
`tests/test_docs_contract.py` checks these functions against the table in
`docs/terms.md`. IT DOES NOT: that test recomputes the table from `math.atanh` and
never imports this module, and `test_mechanism.py` compares against a third
hardcoded copy of the table. So nothing tied the code to the glossary, which is
how the `epsilon_tau` defect below survived at every `delta_u` != 1.
`test_calibration_is_checked_against_the_table_parsed_from_the_glossary` now parses
the table out of `docs/terms.md` and checks these functions against it.
"""

from __future__ import annotations

import math


def eta_tau(tau: float, *, delta_u: float = 1.0) -> float:
    """The tilt magnitude: `arctanh(tau) / delta_u`.

    `delta_u` is the RANGE of the utility being tilted with. It defaults to 1
    because `utility.py` builds a utility bounded in [0,1] by construction and
    asserts it at runtime -- but it is an argument, not a constant, so that a
    later utility with a different range has to state it here instead of
    silently inheriting a calibration that no longer holds.
    """
    if not 0.0 <= tau < 1.0:
        raise ValueError(f"tau must be in [0, 1), got {tau}")
    if delta_u <= 0:
        raise ValueError(f"delta_u must be positive, got {delta_u}")
    return math.atanh(tau) / delta_u


def epsilon_tau(tau: float, *, delta_u: float = 1.0) -> float:
    """The local DP budget: `2 * arctanh(tau)`. INDEPENDENT of the utility range.

    ⚠ Never `arctanh(tau)` alone -- that halves the stated budget and understates
    the leakage. It comes from the pointwise likelihood-ratio bound: the tilt
    contributes a factor `exp(eta_tau * delta_u)` and the ratio of normalisers
    contributes another, so the ratio spans `exp(2 * eta_tau * delta_u)` =
    `exp(2 * arctanh(tau))`.

    ⚠ FIXED 2026-09-30, AND IT WAS THE ERROR THIS FUNCTION EXISTS TO PREVENT. The
    first version returned `2.0 * eta_tau(tau, delta_u=delta_u)`, i.e.
    `2*arctanh(tau)/delta_u`. That is wrong, and wrong in the forbidden direction:
    `delta_u` CANCELS here, because a utility spanning `delta_u` is tilted by an
    exponent of `eta_tau` whose range is `eta_tau * delta_u = arctanh(tau)`. At
    `delta_u` = 2 the old code returned exactly `arctanh(tau)` -- literally the
    expression `docs/terms.md` forbids in capitals as "a silent understatement of
    the leakage". It was unreachable (nothing passes `delta_u` here) and untested
    at any `delta_u` != 1, so the single scenario the parameter exists for was the
    one that corrupted the headline privacy number. `delta_u` is now validated and
    deliberately unused.
    """
    if delta_u <= 0:
        raise ValueError(f"delta_u must be positive, got {delta_u}")
    if not 0.0 <= tau < 1.0:
        raise ValueError(f"tau must be in [0, 1), got {tau}")
    return 2.0 * math.atanh(tau)


def p_succ_bound(tau: float) -> float:
    """The equal-prior bound on an optimal binary attacker: `(1 + tau)/2`.

    An EQUAL-PRIOR bound. An attacker with an informative prior over cohort
    membership is a different question this does not answer. A measured
    `attacker_accuracy` above this at the same `tau` means the implementation is
    wrong, not the theorem -- and one below it is not evidence of anything
    (standing constraint 7).
    """
    if not 0.0 <= tau < 1.0:
        raise ValueError(f"tau must be in [0, 1), got {tau}")
    return (1.0 + tau) / 2.0
