# terms.md — glossary

Search this before naming anything. Every symbol, quantity and metric is defined
once, with provenance and status. If a name is not here, it is not an approved
name. Symbols follow `paper/Private_Genome_Path_Release.pdf`; the identifier in
the heading is the form used in code, tables and filenames.

**Status values:** `PRIMARY` | `SECONDARY` | `DIAGNOSTIC` | `DEMOTED` | `RETIRED`

**Tags:** [theory] [privacy] [utility] [mechanism] [implementation] [attack] [dataset]

---

## The privacy parameter and its derived constants

### tau  [theory] [privacy]
```
Is:                   The pre-specified privacy parameter, tau in [0, 1). It is the
                      ceiling on total variation distance between the release
                      distributions induced by ANY two admissible input paths.
                      tau = 0 means the release ignores the target entirely.
Computed by:          not computed — CHOSEN, and pre-registered before any attack
                      outcome is inspected. See constraint 5 in CLAUDE.md.
CHANGES MEANING WITH: nothing. tau is definitional. But every quantity that DEPENDS
                      on it (utility_retained, attacker_accuracy) is meaningless
                      without its tau stamped alongside.
Valid vs:             nothing. tau is never compared across runs; it is the axis
                      other things are compared along.
Status:               PRIMARY
Aliases:              "the privacy budget" — WRONG and forbidden. The budget is
                      epsilon_tau. tau is a TV bound, not a budget.
```

### eta_tau  [theory] [mechanism]
```
Is:                   The tilt magnitude, eta_tau = arctanh(tau) = 0.5*log((1+tau)/(1-tau)).
                      It is the exponent multiplier applied to the utility in the
                      release distribution. It limits how strongly the private target
                      may bend the public baseline.
Computed by:          math.atanh(tau)  (see the calibration table below)
CHANGES MEANING WITH: the utility range. If the utility is NOT bounded in [0,1] but
                      spans delta_u, the correct calibration is arctanh(tau)/delta_u.
                      Using arctanh(tau) with an unnormalized utility silently
                      breaks the theorem — this is the highest-severity error class
                      in the project.
Valid vs:             nothing; it is a deterministic function of tau.
Status:               PRIMARY
Aliases:              "eta". Do not write "the temperature" — the sign convention
                      is opposite to the usual one.
```

### epsilon_tau  [theory] [privacy]
```
Is:                   The local differential privacy budget the mechanism satisfies
                      over the admissible input-path domain:
                      epsilon_tau = 2*arctanh(tau) = log((1+tau)/(1-tau)).
                      It is EXACTLY TWICE eta_tau. It comes from the pointwise
                      likelihood-ratio bound in Theorem 1, which is strictly
                      stronger than the TV statement.
Computed by:          2 * math.atanh(tau)
CHANGES MEANING WITH: nothing, but note it is the LDP budget over INPUT PATHS, not
                      over genomes or individuals. It says nothing about what a
                      person's other released statistics leak.
Valid vs:             other epsilon-LDP mechanisms over the same input domain ONLY.
                      Never compared to a central-DP epsilon.
Status:               SECONDARY — reported for readers who want the DP framing.
                      The paper leads with tau because of its direct identification
                      interpretation.
Aliases:              "epsilon". NEVER write epsilon_tau = arctanh(tau) — that
                      halves the stated budget and is a silent understatement of
                      the leakage. tests/test_docs_contract.py checks this relation.
```

### p_succ_bound  [theory] [privacy] [gate]
```
Is:                   The upper bound on an optimal binary attacker's success
                      probability under equal priors: (1 + tau)/2. "Did this release
                      come from the target, or from a cohort member?" cannot be
                      answered better than this by ANY attacker, including one who
                      knows the graph, the cohort, the utility function, tau and the
                      full implementation — everything except the private random draw.
Computed by:          (1 + tau) / 2
CHANGES MEANING WITH: the prior. It is an EQUAL-PRIOR bound. An attacker with an
                      informative prior over cohort membership is a different
                      question this bound does not answer.
Valid vs:             a measured attacker_accuracy at the SAME tau. The measured
                      value must not exceed it; if it does, the implementation is
                      wrong, not the theorem.
Status:               PRIMARY — this is the gate. Mark it * in every table.
Aliases:              "the identification bound"
```

### Calibration table

Recomputed and checked by `tests/test_docs_contract.py`. Values rounded to 6 dp.

| `tau` | `eta_tau` | `epsilon_tau` | `p_succ_bound` |
|---|---|---|---|
| 0.00 | 0.000000 | 0.000000 | 0.5000 |
| 0.01 | 0.010000 | 0.020001 | 0.5050 |
| 0.05 | 0.050042 | 0.100083 | 0.5250 |
| 0.10 | 0.100335 | 0.200671 | 0.5500 |
| 0.25 | 0.255413 | 0.510826 | 0.6250 |
| 0.50 | 0.549306 | 1.098612 | 0.7500 |
| 0.90 | 1.472219 | 2.944439 | 0.9500 |

---

## The mechanism

### baseline_model  [mechanism] [theory]
```
Is:                   R_D, the target-independent probability distribution over the
                      output support Y_D, constructed ONLY from the public graph G
                      and the public cohort D. Realized as the cohort haplotype HMM
                      (rho, A_2:T). Every y in Y_D must have R_D(y) > 0.
Computed by:          PLANNED — src/ranpanmixer/hmm.py:CohortHMM
CHANGES MEANING WITH: the cohort D, the block segmentation, and the state-collapse
                      policy for duplicate haplotypes. A baseline fitted with ANY
                      target-derived quantity is not R_D and voids Theorem 1.
Valid vs:             another baseline over the same Y_D and same block segmentation.
Status:               PRIMARY
Aliases:              "R_D", "the prior", "the cohort model". Never "the null model" —
                      it is not a null hypothesis.
```

### output_support  [mechanism] [theory]
```
Is:                   Y_D, the fixed set of cohort-supported output paths the
                      mechanism may emit. May be complete cohort haplotype paths,
                      synthetic paths from the cohort model, or mosaics of
                      cohort-supported local traversals.
Computed by:          PLANNED — implied by the HMM state space and transition sparsity
CHANGES MEANING WITH: the transition sparsity pattern (which concatenations are
                      graph-valid). MUST be fixed independently of the target:
                      target-specific candidate lists introduce target-dependent
                      zero probabilities and generally invalidate the guarantee.
Valid vs:             nothing — two mechanisms with different Y_D are not comparable.
Status:               PRIMARY
Aliases:              "Y_D", "the support"
```

### u_path  [utility] [mechanism]
```
Is:                   The bounded utility u(p, y) in [0,1]: how well releasing y
                      preserves utility for input path p. Blockwise construction:
                      u(p,y) = sum_t beta_t * phi_t(p_t, y_t).
Computed by:          PLANNED — src/ranpanmixer/utility.py
CHANGES MEANING WITH: the choice of phi_t, the weights beta_t, and above all its
                      RANGE. The theorem assumes range exactly [0,1]; see eta_tau.
Valid vs:             the same utility function at another tau. NEVER across
                      different phi_t or beta_t — that is a different quantity
                      wearing the same name.
Status:               PRIMARY
Aliases:              "u". Never "similarity" — phi_t is the similarity; u is its
                      weighted aggregate.
```

### phi_t  [utility]
```
Is:                   The local per-block similarity phi_t(a, b) in [0,1] between the
                      target's traversal and a candidate traversal at block t. The
                      paper's worked example is a population-weighted k-mer Jaccard:
                      sum_x w(x) min(1[x in K(a)], 1[x in K(b)]) over
                      sum_x w(x) max(...), where w downweights ubiquitous human
                      k-mers and upweights population-informative ones.
Computed by:          PLANNED — src/ranpanmixer/utility.py
CHANGES MEANING WITH: the k-mer length, the weighting function w, and whether K(.)
                      is a set or a multiset. Two phi_t with different k are
                      different metrics; do not average them.
Valid vs:             the same phi_t on the same block segmentation.
Status:               PRIMARY
Aliases:              "block similarity", "local agreement"
```

### beta_t  [utility]
```
Is:                   Non-negative per-block weights with sum_t beta_t = 1. They
                      decide which parts of the path the mechanism tries to preserve.
Computed by:          PLANNED
CHANGES MEANING WITH: whether they are uniform (1/T) or length-weighted. The
                      normalization is what keeps u_path inside [0,1]; weights that
                      do not sum to 1 break the calibration, not just the emphasis.
Valid vs:             the same weighting scheme.
Status:               SECONDARY
Aliases:              none
```

### psi_t  [implementation] [mechanism]
```
Is:                   The local potential psi_t(j) = exp(eta_tau * beta_t * phi_t(p_t, h_t,j)).
                      The ONLY place the private target enters the sampler. The cohort
                      HMM (rho, A) stays fixed.
Computed by:          PLANNED — src/ranpanmixer/sampler.py
CHANGES MEANING WITH: eta_tau, and therefore tau. Potentials cached across tau values
                      are a correctness hazard.
Valid vs:             nothing; an intermediate quantity.
Status:               DIAGNOSTIC
Aliases:              "the tilt potential", "the emission"
```

### release_distribution  [mechanism] [theory]
```
Is:                   Q^tau_p(y) = R_D(y) * exp(eta_tau * u(p,y)) / Z_p, the
                      distribution the sanitized path is drawn from for input path p.
Computed by:          PLANNED — never materialized explicitly; sampled from via ffbs.
CHANGES MEANING WITH: tau, R_D, u_path. All three must be stamped on any figure
                      derived from it.
Valid vs:             another release distribution at matched tau, R_D and u_path.
Status:               PRIMARY
Aliases:              "Q_p", "the tilted distribution"
```

### log_z_p  [implementation] [diagnostic]
```
Is:                   The log normalizing constant, log Z_p = sum_t log c_t, where
                      c_t is the forward-pass per-position normalizer. Theory pins it:
                      0 <= log Z_p <= eta_tau, for every input path p.
Computed by:          PLANNED — src/ranpanmixer/sampler.py, forward pass
CHANGES MEANING WITH: tau. The bound scales with eta_tau, so a log_z_p quoted
                      without its tau cannot be checked.
Valid vs:             its own theoretical bracket. This is the cheapest correctness
                      assertion in the whole implementation: if log Z_p leaves
                      [0, eta_tau], the sampler is wrong. Assert it, do not report it.
Status:               DIAGNOSTIC
Aliases:              "the normalizer", "Z_p" (Z_p is the unlogged form)
```

### ffbs  [implementation] [mechanism]
```
Is:                   Forward-filtering / backward-sampling: the exact sampler for
                      Q^tau_p. Forward pass computes per-position normalized alpha_t;
                      backward pass draws Z_T ~ Categorical(alpha_T) then walks back
                      sampling Z_t | Z_t+1. Algorithm 1 of the proposal.
Computed by:          PLANNED — src/ranpanmixer/sampler.py
CHANGES MEANING WITH: nothing — it is EXACT, not approximate, given the potentials.
Valid vs:             a brute-force enumeration of Q^tau_p on a toy graph. That
                      equivalence is the primary implementation gate.
Status:               PRIMARY
Aliases:              "the sampler", "FFBS". NEVER confused with Viterbi — see the
                      RETIRED section.
```

### stay_switch  [implementation] [mechanism]
```
Is:                   PanMixer's transition model: probability a_t of remaining in
                      the same donor state, b_t of switching, for every other state.
                      Its rank-one-plus-diagonal structure is what reduces the
                      forward recursion from O(TK^2) to O(TK).
Computed by:          PLANNED — src/ranpanmixer/hmm.py
CHANGES MEANING WITH: whether b_t is the per-target or the total switch mass.
                      Getting that wrong rescales every transition by K.
Valid vs:             a dense-transition HMM giving identical forward values. That
                      equivalence is the gate on the O(TK) optimization.
Status:               PRIMARY
Aliases:              "the PanMixer transitions"
```

### T_blocks  [implementation] [dataset]
```
Is:                   T, the number of ordered anchor sites or graph blocks a path
                      is divided into.
Computed by:          PLANNED — set by the block segmentation of G
CHANGES MEANING WITH: the segmentation policy (LD blocks vs graph subpaths vs fixed
                      anchors). Runtime and privacy granularity both scale with it.
Valid vs:             the same segmentation.
Status:               SECONDARY
Aliases:              "T", "blocks", "sites". Pick one per table and stay with it.
```

### K_states  [implementation] [dataset]
```
Is:                   K, the number of cohort haplotype states. At most ~2n for a
                      diploid cohort of n individuals, before duplicate haplotypes
                      are collapsed.
Computed by:          PLANNED
CHANGES MEANING WITH: whether duplicates have been collapsed. A K quoted before
                      collapse and a K quoted after are different numbers; runtime
                      claims must say which.
Valid vs:             the same cohort at the same collapse policy.
Status:               SECONDARY
Aliases:              "K", "donor states"
```

### diploid_utility  [utility] [theory]
```
Is:                   For a diploid target, the haplotype pair is ONE private input:
                      u_diploid(p, y) = (u(p1, y1) + u(p2, y2)) / 2. The average keeps
                      the joint utility in [0,1], so the same tau covers the complete
                      diploid release.
Computed by:          PLANNED
CHANGES MEANING WITH: the aggregator. A SUM instead of a MEAN leaves [0,1] and
                      doubles the effective tilt — the guarantee no longer holds.
Valid vs:             a haploid release only after the aggregator is stated.
Status:               PRIMARY
Aliases:              none
```

---

## Measured quantities

These are estimates from finite samples. The theory quantities above are exact;
these are not, and every one of them needs an `n` and a `tau`.

### utility_retained  [utility] [evaluation]
```
Is:                   u_path(p_g, released_path) achieved by an actual release, for
                      the true target path p_g. Averaged over independent draws.
Computed by:          PLANNED — src/ranpanmixer/eval.py
CHANGES MEANING WITH: tau (monotonically), the utility function, AND the denominator
                      (see the two measurement paths below). Report with n draws.
Valid vs:             the same utility function at another tau, against the floor
                      (tau = 0, an untilted draw from R_D) and the ceiling
                      (the target's own path, u = 1 by construction).
Status:               PRIMARY
Aliases:              "utility". Never quoted bare without its tau.
```

### attacker_accuracy  [attack] [evaluation] [gate]
```
Is:                   The empirical success rate of an IMPLEMENTED binary attacker
                      asked to decide whether a release came from the target or from
                      a cohort member, under equal priors.
Computed by:          PLANNED — src/ranpanmixer/attacks.py
CHANGES MEANING WITH: the attacker. It is a LOWER bound on distinguishability and
                      nothing else. It cannot validate the mechanism; it can only
                      falsify an implementation.
Valid vs:             p_succ_bound at the SAME tau — as a one-sided check that the
                      implementation is not leaking more than the theory allows.
                      Never as evidence the mechanism is "more private than tau".
Status:               PRIMARY — mark * in tables.
Aliases:              "attack accuracy", "identification rate"
```

### tv_empirical  [privacy] [evaluation] [diagnostic]
```
Is:                   An ESTIMATE of TV(Q^tau_p, Q^tau_q) from finite samples.
Computed by:          PLANNED — on toy graphs by exact enumeration; at scale it is
                      not directly estimable and must be reported as such.
CHANGES MEANING WITH: the estimator. Plug-in TV over a large support is severely
                      BIASED UPWARD at finite n — an estimate above tau on a big
                      support is the expected behaviour of the estimator, not a
                      violation of the theorem. Only exact enumeration on a toy
                      support can falsify Theorem 1.
Valid vs:             tau, on a support small enough to enumerate exactly. Nothing else.
Status:               DIAGNOSTIC — never a headline number.
Aliases:              "measured TV". Never write "the TV distance" for the estimate;
                      that name belongs to the exact quantity.
```

---

## THE TWO MEASUREMENT PATHS

These must never be mixed in one comparison. Naming one when you mean the other
inverts conclusions.

- **PATH A — ALL BLOCKS.** Denominator: every block in the path. "How much utility
  does a release retain overall?" Every utility metric MUST be defined here, and
  must be defined for a release that changed nothing.
- **PATH B — SWITCHED BLOCKS ONLY.** Denominator: blocks where the sampled donor
  state differs from the target's own traversal. "Given that the mechanism
  intervened, what did it cost?" Systematically LOWER than Path A, because the
  unchanged blocks contributing u = 1 are excluded.
- **A Path-B number quoted over all blocks is not a more conservative result — it
  is a DIFFERENT and usually wrong quantity.** At small tau most blocks do not
  switch, so the two diverge exactly where the mechanism is most useful.

The same split applies to targets: metrics over ALL targets versus over targets
carrying rare alleles. State which, always.

---

## Retired

`viterbi_release` — RETIRED 2026-08-29, before use. MAP/Viterbi decoding of the
tilted HMM returns the single highest-probability state sequence. It is NOT a draw
from Q^tau_p and carries NO privacy guarantee whatsoever, while producing output
that looks entirely plausible. Recorded here because it is the natural thing to
reach for in an HMM codebase and its failure is silent. Do not re-add.
