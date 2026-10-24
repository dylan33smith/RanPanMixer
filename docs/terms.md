# terms.md — glossary

Search this before naming anything. Every symbol, quantity and metric is defined
once, with provenance and status. If a name is not here, it is not an approved
name. Symbols follow `paper/Private_Genome_Path_Release.pdf`; the identifier in
the heading is the form used in code, tables and filenames.

**Status values:** `PRIMARY` | `SECONDARY` | `DIAGNOSTIC` | `DEMOTED` | `RETIRED`

**Tags:** [theory] [privacy] [utility] [mechanism] [implementation] [attack] [dataset]

---

## The arm — the stamp every number needs first

### arm  [method] [evaluation]
```
Is:                   WHICH OF THE PAPER'S TWO SECTIONS a number belongs to. It is
                      the first stamp on every figure, before tau, the denominator
                      and n, because it decides what is true of the TARGET and
                      therefore what half the glossary means. Two values, and only
                      two:
                      **arm1** — the HEAD-TO-HEAD. The target is one of the 44 HPRC
                      cohort members. We release the obfuscated paths of all 44,
                      DIPLOID, over the 22 autosomes, and run PanMixer's own
                      evaluations so the comparison is direct. Requires per-target
                      leave-one-out at K = 86: BOTH of the target's haplotypes leave
                      the donor panel (see `loo_cohort`). Every leave-one-out,
                      reduced-K and cohort-distortion quantity belongs to this arm
                      alone, and cohort-member protection is in scope HERE ONLY.
                      **arm2** — the EXTERNAL target. The target is outside G, D,
                      every panel and the attack database, so the full cohort is the
                      donor panel (K = 88) and nothing is held out. Nothing in the
                      mechanism is a function of the target but the tilt.
Computed by:          not computed — it is a property of the RUN, set by the
                      configuration and carried on every artifact and every table
                      column. Settled by Dylan 2026-10-08 (the two-arm scope change).
CHANGES MEANING WITH: nothing; it is definitional. But it changes the meaning of
                      `u_ceiling`, `target_fidelity`, `beta_v`'s scope, `loo_cohort`,
                      `K_states`, `af_loss`, `ld_loss`, `p_succ_bound`'s pair and
                      `new_haplotypes`. ⚠ **AN UNSTAMPED NUMBER IS UNINTERPRETABLE,
                      NOT IMPRECISE** — the 2026-09-28 scope correction and the
                      2026-10-08 rewrite were both about exactly this confusion.
Valid vs:             nothing. The two arms are NOT comparable to each other on any
                      ceiling-bearing axis: each number is read against ITS OWN arm's
                      ceiling and floor. Both arms DO run on one consistent GRCh38
                      build (decision F) — the thing PanMixer's own release got wrong
                      — so a site set is shared even where a ceiling is not.
Status:               PRIMARY — every table column names it.
Aliases:              "the head-to-head arm" = arm1; "the external-target arm" =
                      arm2. ⚠ "the comparison arm" is the OLD name for arm1 from when
                      it was optional; entries still carrying `COMPARISON-ARM ONLY`
                      are stale by that much. Never write "the arm" with no value.
```

---

## The privacy parameter and its derived constants

### tau  [theory] [privacy]
```
Is:                   The pre-specified privacy parameter, tau in [0, 1). It is the
                      ceiling on total variation distance between the release
                      distributions induced by ANY two admissible input paths.
                      tau = 0 means the release ignores the target entirely.
Computed by:          not computed — CHOSEN, and pre-registered before any attack
                      outcome is inspected. See constraint 6 in CLAUDE.md (it was
                      constraint 5 before the 2026-10-08 two-arm rewrite, and it now
                      also pre-registers the attacker set and the attack database).
                      ✅ THE tau GRID IS WRITTEN DOWN (2026-10-10), in
                      paper/prereg_2026-10-10.md section 1: SEVEN points, 0, 0.05,
                      0.10, 0.25, 0.50, 0.90 and 0.99, one file covering both arms.
                      Each point's ROLE is part of the pin and must be carried into
                      every table: tau = 0 is the mandatory FLOOR (an untilted draw
                      from baseline_model, zero target information, the control that
                      stops a cohort-level fidelity figure reading as a free lunch);
                      tau = 0.05 to 0.25 is the ONLY range claimed as an operating
                      setting; tau = 0.50 is reported as context, not claimed; tau =
                      0.90 and 0.99 are BOUND-CHECKING ONLY.
                      ✅ AND THE FILE IS IN FORCE: SIGNED 2026-10-10 by Dylan,
                      git_rev b5a220f122ca. The grid is FROZEN under constraint 6 — it
                      does not change because a result is disappointing; a change
                      means a NEW dated file and this one is never edited.
                      ⚠⚠⚠ **`tau` IS NO LONGER THE TOOL'S KNOB (Dylan, 2026-10-24).**
                      It remains fully defined here and in the proof, and every
                      relation below it holds unchanged — but the method, the tool
                      and the results are parameterised by `eta_tau`, because
                      `tau`'s whole usable range moves the release about four cells
                      in 7.6 million and `tanh(eta)` is not representable below 1
                      past `eta_tau` ~ 19. The pre-registered grid is SUPERSEDED
                      (paper/superseded_2026-10-24_prereg_withdrawn.md); the seven
                      points below are now history, not a constraint. `tau` is
                      still the right quantity to DEFINE the guarantee with, and
                      "the `tau` range buys nothing" is a reportable finding.
                      ⚠⚠ THREE FURTHER POINTS EXIST AND ARE **NOT** PRE-REGISTERED:
                      tau = 0.999, 0.9999 and 0.99999, run 2026-10-23 at the PI's
                      request and AFTER the signed grid's outcomes were seen. They
                      are EXPLORATORY, live in logs/B_tau_0p999_0p9999_0p99999/ and
                      are reported APART from the pre-registered result; promoting
                      any of them needs a NEW dated pre-registration file. They are
                      named in SEED_GRID in scripts/arm1_release.py and
                      scripts/arm2_attack.py so each has a FIXED seed index — see
                      docs/bugs.md for why a list position is not one. ✅ What they
                      establish: eta_tau = arctanh(tau) is LOGARITHMIC, so 1 - tau
                      shrinking 1000x multiplies the tilt only 2.3059x, and the
                      measured movement is linear in eta_tau to three decimals at
                      all nine points. tau is a very weak lever, by construction.
                      ⚠ Field 6 was reopened and re-closed on 2026-10-10: the attack
                      design it had pinned was ONE-SIDED, drawing only from Q_p, and
                      so could not measure the equal-prior quantity p_succ_bound
                      bounds at all. Re-closed at the BALANCED 25-draw design, delta
                      unchanged at 0.056 because averaging two arms is more precise
                      per observation.
                      [INCORRECT] - The superseded text said the grid "IS NOT PINNED
                      (2026-10-08)" and that the pre-registration "does not exist
                      yet"; it was unblocked by B-EVL-exact-endpoint on 2026-10-10.
                      ⚠ That endpoint must be E[differ] via the MARGINALS route; the
                      finite-difference derivative route returns a NEGATIVE variance
                      at a genome-wide share and is forbidden there.
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
                      probability under equal priors: (1 + tau)/2. The question it
                      bounds is over INPUT PATHS, which is what Theorem 1 is stated
                      over: "did this release come from input path p, or from input
                      path q?", for ANY two admissible inputs at the SAME release
                      slot. It cannot be answered better than this by ANY attacker,
                      including one who knows the graph, the cohort, the utility
                      function, tau and the full implementation — everything except
                      the private random draw.
                      ⚠ **READ THE PAIR PER ARM (restated 2026-10-08).** In ARM 2, p
                      is the external target and q is any other admissible path. In
                      ARM 1 the target IS one of the 44 cohort members, so the older
                      phrasing of this field — "did this release come from the target,
                      or from a cohort member?" — is VACUOUS there: both branches name
                      the same event. Arm 1's pair is member k's own true diploid path
                      against an admissible DECOY path at the same release slot k. The
                      bound is unchanged either way (it is over pairs of admissible
                      input paths, not over populations); only the sentence that
                      explained it was Arm-2-only, and it is the sentence the paper
                      will quote. See `two_target_indistinguishability`.
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
Computed by:          src/mechanism/hmm.py:CohortHMM (built 2026-09-30). ⚠ Realized as the
                      Li-Stephens chain CONDITIONED on never occupying a donor with no
                      called allele -- see the missing_policy correction below.
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
Computed by:          src/mechanism/hmm.py, implied by the state mask and the switch
                      scalars; enumerated exactly by src/mechanism/exact.py on the toy.
                      ⚠ A Delta_x = 0 step makes Y_D strictly smaller: switching donors
                      there has probability 0, so those output paths are not in Y_D at
                      any tau. Asserted target-independent over the whole tau grid.
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
Computed by:          src/mechanism/utility.py:u_of_output (built 2026-09-30). The [0,1]
                      bound is ASSERTED at runtime, not documented -- _assert_bounded
                      raises rather than returning an out-of-range value.
CHANGES MEANING WITH: the choice of phi_t, the weights beta_t, and above all its
                      RANGE. The theorem assumes range exactly [0,1]; see eta_tau.
Valid vs:             the same utility function at another tau. NEVER across
                      different phi_t or beta_t — that is a different quantity
                      wearing the same name.
Status:               PRIMARY
Aliases:              "u". Never "similarity" — phi_t is the similarity; u is its
                      weighted aggregate.
```


### tv_pinsker  [theory] [privacy] [evaluation]
```
Is:                   The Pinsker UPPER BOUND on the total variation between the release
                      distributions of two inputs, `sqrt(kl_exact / 2)`, clipped at 1.
                      `acc_ceiling` is `(1 + tv_pinsker) / 2`.
                      ⚠⚠ **IT IS A BOUND, NOT A MEASUREMENT.** The KL it is derived from
                      is exact; the total variation is not. Pinsker is an inequality, so
                      the TRUE total variation is SMALLER than this and the true attacker
                      ceiling is closer to 0.5 than `acc_ceiling` prints. The looseness
                      runs in the direction that FLATTERS US, which is why the caveat is
                      printed on every run and must travel with the number.
                      ⚠ Never call it "the measured TV". "measured TV" is an alias of
                      `tv_empirical`, a sampled diagnostic and a different quantity.
Computed by:          scripts/exact_kl.py:140, `min(1.0, sqrt(kl / 2.0))`, from
                      `kl_exact` computed by the state-marginals route.
CHANGES MEANING WITH: everything `kl_exact` changes with -- the arm, `tau`, and above all
                      `beta_v`'s SCOPE AND SHARE. ⚠⚠ The recorded values are chr21 at
                      `beta_share` = 1, the WHOLE budget on one chromosome, because
                      scripts/exact_kl.py calls `beta_v` with no `share` and takes the
                      default of 1.0. That is 66.6x more generous per position than the
                      genome-wide release, so these figures are NOT the genome-wide
                      release's and must never be quoted as such.
Valid vs:             `p_succ_bound` at matched `tau` and matched scope, as a gate.
                      ⚠ Not against a SAMPLED attacker accuracy as though they were the
                      same kind of quantity: this bounds every attacker including
                      unwritten ones, while a measured accuracy is a LOWER bound on
                      distinguishability (constraint 7).
Status:               REPORTED. **Measured 0.000e+00 to 3.109e-03 across `tau` = 0 to
                      0.99** (ARM 2, chr21, `beta_share` = 1, T = 332,521, K = 88,
                      n = 10 targets, max over targets; `logs/B_exact_kl.md`).
                      ⚠⚠ **THE GENOME-WIDE VALUE IS NOT MEASURED.**
Aliases:              "the Pinsker ceiling". Never "the measured TV", never "the TV".
```

### u_ceiling  [utility] [gate]
```
Is:                   The largest u_path a GIVEN target could achieve at any tau:
                      sum of beta_v over the positions where that target has a call.
                      It is the honest denominator for utility_retained, and it is
                      BELOW 1 whenever the target has no call at a weighted position.
Computed by:          src/mechanism/utility.py:u_ceiling (built 2026-09-30)
CHANGES MEANING WITH: the ARM first, then the beta_v SCOPE (deferred decision R1)
                      and the target.
                      ⚠ **THE TWO ARMS HAVE DIFFERENT CEILINGS AND THEY ARE NOT
                      COMPARABLE TO EACH OTHER** (2026-10-08). In ARM 2 the ceiling is
                      set by what the graph can READ of an external target. In ARM 1
                      the target's own assembly is one of the 44 G was built from, so
                      it is readable essentially everywhere and the binding limit is a
                      DIFFERENT one: DONOR-PANEL REACHABILITY under per-target
                      leave-one-out at K = 86. An allele that no REMAINING donor
                      carries cannot be emitted at any tau, including tau -> 1, and
                      the target's own two haplotypes are exactly what left the panel.
                      **Arm 1's ceiling is UNMEASURED** — it needs the LOO panel, so
                      it cannot be inferred from the Arm-2 figure below. Never put the
                      two arms' ceilings in one column.
                      ⚠⚠ **SUPERSEDED BY THE CHAIN REBUILD (2026-10-12). THE LIVE
                      VALUE IS 0.692509**, measured on HG00096, chr21, scope="chain",
                      ARM 2, over the decision-K chain of T = 332,521, where
                      30.749104% of the weight sits on positions the target does not
                      call. The rebuild added 26,635 positions of which only 4,362
                      are called, so the called fraction fell 0.829871 -> 0.776516.
                      ⚠ The figure below is CORRECT FOR ITS BUILD and is kept because
                      many recorded numbers divide by it: measured 2026-09-30,
                      **ARM 2** (external stand-in target), on
                      HG00096, chr21, scope="chain": **0.815455** over 253,846 of
                      305,886 chain positions -- independently reproducing the 0.815
                      that R1 predicted from the weight arithmetic. Restricting the
                      live arrays to the old span reproduces both exactly, which pins
                      the change to the BUILD rather than to the arithmetic. Under
                      scope="chain_readable" it is 1.000 for a target drawn from the
                      same panel `readable` was built from.
                      ⚠ **AND THAT NUMBER IS CHR21-NORMALISED.** It was measured with
                      beta_v summing to 1 over chr21's chain alone. Under decision B
                      beta_v sums to 1 GENOME-WIDE over the 22 autosomes and chr21
                      carries only a SHARE of the budget (`beta_share`), so the
                      genome-wide ceiling is a different number and is UNMEASURED. Do
                      not carry 0.815455 into a genome-wide table.
Valid vs:             the same ARM, the same scope and the same target. ⚠ NEVER quote
                      utility_retained against an implicit ceiling of 1 -- that
                      understates the mechanism by exactly this factor. And never
                      quote it against this number without saying so, which
                      overstates it. Report both, always.
Status:               PRIMARY — mark * in tables.
Aliases:              none. ⚠ NOT target_fidelity's ceiling, which is a different
                      loss in the opposite direction: that one is about target
                      variants with NO record in G — **ARM 2 ONLY** (about 9.9% of
                      non-reference calls, 20.5% of the -log f information). In ARM 1
                      the target's own assembly helped build G, so that particular
                      loss is near zero and Arm 1's binding ceiling is the LOO panel
                      instead. This one is about graph positions where we cannot READ
                      the target. `docs/plan.md` item 1b keeps the three losses apart;
                      conflating them inverts the conclusion, and conflating them
                      ACROSS ARMS inverts it twice.
```

### phi_v  [utility]
```
Is:                   The local per-position similarity in [0,1] between the target's
                      allele and a candidate donor's allele at chain position v.
                      ✅ **UNDER MODEL C (decided 2026-09-28) THIS COLLAPSES TO AN
                      INDICATOR**: phi_v = 1 if the drawn donor carries the target's
                      allele at v, else 0. There is no block to aggregate over, no
                      per-block normaliser, and no aggregation choice to make.
Computed by:          src/mechanism/utility.py:phi_v (built 2026-09-30)
CHANGES MEANING WITH: the missing_policy (an excluded donor has no phi_v at all), and
                      whether the target's own allele is missing at v — in which case
                      phi_v is **ZERO for every donor**, so the tilt POTENTIAL is 1
                      everywhere and the position is a pure prior draw. ⚠ That is NOT
                      a wildcard match (settled; `src/mechanism/utility.py:phi_v`
                      returns zeros at an uncalled target allele): rewarding donors
                      for agreeing with an unknown value would invent information.
Valid vs:             the same phi_v on the same site axis.
Status:               PRIMARY
Aliases:              "local agreement". ⚠ RENAMED from `phi_t` 2026-09-29. The old
                      name came from the per-block design; under model C the index
                      is a variant, not a block. `phi_t` in older entries and in
                      docs/memory.md means this quantity.
⚠ THE k-MER JACCARD IS NOT THIS. The proposal's worked example — a population-weighted
                      k-mer Jaccard, sum_x w(x) min(1[x in K(a)], 1[x in K(b)]) over
                      sum_x w(x) max(...) — is RULED OUT as the in-sampler potential
                      for two independent reasons. Structural: a k-mer straddling two
                      positions depends on both at once and a Jaccard's denominator is
                      a union over the whole span, so it does not decompose into
                      per-position factors and ffbs cannot draw from the tilted
                      distribution exactly — and an approximate draw carries NO bound.
                      Practical: it needs sequence, and the data model is integer
                      allele indices. ✅ DEFERRED 2026-09-28 to a post-build
                      evaluation/reporting tool, where nothing has to factor. Keep it
                      out of the mechanism.
```

### beta_v  [utility]
```
Is:                   Non-negative per-position weights with sum_v beta_v = 1. They
                      decide which parts of the path the mechanism tries to preserve.
                      Under model C: beta_v = w_v / sum_u w_u with w_v = 1/support(v).
                      ⚠⚠⚠ **THE TOOL'S WEIGHTING IS NOW `target_alt` (Dylan,
                      2026-10-24)**: uniform over the positions where THAT haplotype
                      carries a NON-REF allele, normalised by the genome-wide count
                      of them per haplotype. Built by
                      scripts/eta_sweep.py:build_betas and
                      scripts/exact_utility.py:alt_totals, one definition shared
                      with scripts/alt_beta_sweep.py. Worth **29.7x** the endpoint
                      movement of 1/support(v) at matched tau
                      (logs/B_alt_beta_genome_10targets.md). ⚠ It is
                      TARGET-DEPENDENT. Theorem 1 still holds -- beta >= 0, sums to
                      1, u stays in [0,1] -- but it is a DIFFERENT utility function
                      from the production one and constraint 8 forbids them sharing
                      a column. w_v = 1/support(v) is retained as the comparison.
                      ⚠⚠ **CONTIGUOUS-WINDOW AND TOP-n WEIGHTINGS ARE DEAD**
                      (Dylan, 2026-10-24): dropped as a tool direction, absent from
                      the paper. logs/B_beta_concentration.md and
                      logs/B_topn_sweep_10targets.md are ARCHIVED records of what
                      was measured, cited only to explain why they were dropped.
                      ✅ **THE SUM IS GENOME-WIDE (decision B, 2026-10-08).** beta_v
                      sums to 1 over the whole scope tau covers — the 22 AUTOSOMES
                      TOGETHER (decision D: no chrX, no chrY) — and each chromosome
                      takes a SHARE of that one budget, which is `beta_share`. NOT per
                      chromosome: 22 chromosomes each normalised to 1 would be 22
                      independent releases of the same person and would compose, which
                      standing constraint 5 forbids. Constraint 4 is the rule that
                      stops it arising; constraint 5 is the one it would violate.
                      ⚠ w_v = 1/support(v) is INHERITED, not argued for — deferred
                      decision R3 is still open and is a pre-registration item under
                      constraint 6.
Computed by:          src/mechanism/utility.py:beta_v (built 2026-09-30), which takes an
                      explicit `scope` -- see u_ceiling, `beta_scope` and deferred
                      decision R1.
                      ✅ **THE CODE EXPRESSES A SHARE AS OF DECISION J (2026-10-11).**
                      This paragraph used to read "THE CODE CANNOT EXPRESS A SHARE
                      TODAY, AND SAYS SO". beta_v normalises over the chain positions
                      it is handed, multiplies by `share`, and HARD ASSERTS the result
                      AGAINST THE SHARE -- src/mechanism/utility.py raises "sum beta_v
                      = ..., not the declared share ...; the calibration is void." The
                      22 autosome shares are T_c / T_total by variant count and sum to
                      1.000000 over T = 22,134,562, so `tau` is genome-wide.
                      ⚠ A CHR21 FIGURE AT share = 1 IS STILL CHR21-NORMALISED, and
                      every pre-2026-10-12 beta_v, u_ceiling, E_Q[u] and
                      utility_retained figure on record is one. That is legitimate as a
                      SINGLE-CHROMOSOME release, not as a genome-wide number; read the
                      share off the figure before comparing two. Genome-wide E_Q[u] is
                      measured -- see `e_q_u`.
                      ⚠ support(v) needs NO computation — measured 2026-09-28,
                      the pangenome VCF's own INFO/AN field is bit-identical to
                      np.sum(pangenome != -1, axis=(0,2)) across all 340,849 chr21
                      records.
CHANGES MEANING WITH: whether they are uniform (1/T) or support-weighted, and whether
                      length-weighting is ever adopted. The normalization is what keeps
                      u_path inside [0,1]; weights that do not sum to 1 break the
                      calibration, not just the emphasis.
                      ⚠ AND WITH THE ARM, through the scope (2026-10-08). **ARM 1's
                      scope is SETTLED as "chain" on every branch:** `readable` is
                      produced only by the external-target pipeline, and a
                      readable-off-G scope is degenerate in Arm 1 because the chain
                      axis IS G's own deconstructed site list. **ARM 2's scope is
                      STILL OPEN** (R1/R4): "chain" against a named-public-callset
                      "chain_readable". See `beta_scope`.
                      ⚠ AND WITH THE ARM a second way: under per-target leave-one-out
                      each Arm 1 target has its own support(v) and therefore its own
                      beta_v, so the 44 releases carry 44 DIFFERENT utility functions.
                      Pool measured attack outcomes, never utilities (constraint 8).
Valid vs:             the same weighting scheme, the same arm, and the same
                      normalisation axis (genome-wide, or one chromosome).
Status:               SECONDARY
Aliases:              ⚠ RENAMED from `beta_t` 2026-09-29, same reason as `phi_v`.
⚠ THE BUDGET IS INVARIANT IN T. Since u = sum_v beta_v phi_v with sum beta = 1, u lands
                      in [0,1] however the path is partitioned, and the total tilt
                      between a perfectly-matching and a non-matching path is exp(eta_tau)
                      regardless of T. Partitioning finer does NOT spread the budget
                      thinner: each position moves less and there are proportionally
                      more of them. Measured at tau = 0.5 (eta = 0.5493): per-position
                      tilt 1.0000055 at T = 100,757 (blocks) against 1.0000016 at
                      T = 340,849 (variants), total exp(eta) = 1.7321 either way.
                      Model C is in fact BETTER on utility at matched tau, because a
                      block design must take a whole block from one donor while C can
                      collect partial credit position by position.
⚠ THE REAL CONSEQUENCE OF FINE T IS NUMERICAL. The per-position tilt under C is only
                      **15.1x** the float32 resolution of 1.19e-07, against 45.7x under
                      a block partition. float64 for the stored alphas is load-bearing,
                      not a preference — float32 would consume a meaningful fraction of
                      the target's entire influence.
                      ⚠ CORRECTED 2026-10-10: this read **13.5x**, computed over the
                      340,849-record site axis instead of the chain's T = 305,886. The
                      chain is the denominator, so it is 0.549306/305,886 = 1.80e-06 and
                      15.1x. ⚠ Under decision B (genome-wide `beta_v`) the tilt falls
                      BELOW float32 eps outright — see `docs/bugs.md`.
⚠ THE THIN-BUDGET PROBLEM IS A PROPERTY OF tau, NOT OF T. At tau = 0.5 the best-matching
                      path is only 1.73x more likely than the worst; that is what
                      tau = 0.5 MEANS. exp(arctanh(tau)) = 1.1055 / 1.2910 / 1.7321 /
                      4.3589 at tau = 0.1 / 0.25 / 0.5 / 0.9 (derived). If the achieved
                      utility is too low, the lever is tau — a weaker privacy claim —
                      not a different segmentation.
⚠ REFINED 2026-09-30, MEASURED: tau IS A VERY WEAK LEVER, AND THE LEVER IS sum beta^2.
                      The two paragraphs above are correct and are about the wrong
                      quantity for predicting utility. The tilt can only move mass within
                      the range of u that R_D actually REACHES, and that range is not
                      [0,1]: measured on chr21, HG00096, scope "chain", ARM 2, with
                      beta CHR21-NORMALISED (not the decision-B genome-wide axis),
                      Var_Q(u) = 4.046e-06, so sd(u) = 2.011e-03 around E_Q[u] = 0.6951.
                      These are `e_q_u` and `var_q_u`; both must be re-measured on the
                      genome-wide beta before any endpoint is pre-registered.
                      The achievable gain is about eta_tau * Var = 2.222e-06 at
                      tau = 0.5 — a signal-to-noise ratio against a single draw of about
                      1:900, needing ~8.3e+05 draws to resolve. MEASURED CONSEQUENCE: the
                      released path MOVES between tau = 0 and tau = 0.9 — 8 of 10
                      disjoint-stream draws differ, by a mean of 417 state positions and
                      54 allele positions per haplotype — and the movement is
                      UTILITY-NEUTRAL. ⚠ An earlier version of this entry said the path
                      was "BIT-IDENTICAL ... 0 of 611,772 positions differ"; that was
                      measured on one seed under common random numbers and is FALSE. See
                      the correction in `docs/memory.md` 2026-09-30.
                      Part of the reason is exact: sum_v eta_tau * beta_v = eta_tau, since
                      beta sums to 1, so the expected number of INDEPENDENT flips is of
                      order eta_tau (0.549 at tau = 0.5, 1.47 at tau = 0.9) and CANNOT
                      grow with T — but one flip re-seeds the whole downstream path
                      through the backward recursion, so that bounds the flips, not the
                      number of differing positions.
⚠ AND REWEIGHTING beta_v CANNOT FIX IT (measured 2026-09-30, corrects this entry's own
                      earlier claim that "THE LEVER IS sum beta^2"). The real
                      sum_v beta_v^2 is 4.556723e-06; the 3.27e-06 once quoted here is
                      1/T = mean(beta_v), which equals sum beta^2 only for uniform
                      weights. Concentrating beta 3,000-fold raises sum beta^2 2,200x and
                      Var_Q(u) only 16.5x — the gain stays 55x below the per-draw sd(u) —
                      while E_Q[u] collapses from 0.695 to 0.049. 94.8% of Var_R(u) is
                      between-position COVARIANCE from the mosaic's segment structure,
                      which no beta_v touches. So R3 (rarity weighting instead of call
                      rate) is worth settling on its own merits and is a PRE-REGISTRATION
                      item under constraint 6, but it is NOT the answer to R5. What
                      remains is the tension between exact sampling — phi_v must decompose
                      per position for ffbs — and a utility with real variance. See
                      `docs/plan.md` R5 and `docs/memory.md` 2026-09-30.
```

### beta_scope  [utility] [method]
```
Is:                   WHICH POSITION SET beta_v's normaliser runs over, and therefore
                      what "sum_v beta_v = 1" is a statement about. The axis is
                      GENOME-WIDE over the 22 autosomes (decision B); this names the
                      position set WITHIN that axis. Three values:
                      **chain** — every chain position. Honest, the current default,
                      and SETTLED for ARM 1 on every branch.
                      **chain_readable** — chain positions the target's callset can
                      also speak to. Legitimate only if `readable` is a function of
                      something fixed and public (G itself, or a named public callset
                      pinned in advance), never of the target's own VCF — that is
                      deferred decision R4, and it is ARM 2's open question (R1/R4).
                      **window** — ⚠⚠ RETIRED 2026-10-24, contiguous windows are
                      dropped as a tool direction and are not in the paper; the
                      value is documented below only so an old artifact naming it
                      remains readable. Was: a contiguous-window scope for decision A's
                      SECONDARY utility arm. ⚠ NOT IMPLEMENTED: BETA_SCOPES at
                      src/mechanism/utility.py:51 is ("chain", "chain_readable") only,
                      so this value is a glossary placeholder until the code admits it.
Computed by:          src/mechanism/utility.py:beta_v's `scope` argument, validated
                      against BETA_SCOPES (src/mechanism/utility.py:51,136-137).
                      `Utility.scope` carries it onto every figure.
CHANGES MEANING WITH: the arm, and whether `readable` is target-derived. A
                      target-derived scope is a standing-constraint-1 violation, not a
                      utility choice: it makes the normaliser a function of the
                      private input. ⚠ It also moves `u_ceiling`, so no ceiling is
                      quotable without it.
Valid vs:             the same scope only. Two runs under different scopes are
                      different quantities wearing one name.
Status:               PRIMARY — a pre-registration item under constraint 6 (R1/R4),
                      settled for Arm 1 and OPEN for Arm 2.
Aliases:              "the beta scope", `scope` in code. Never "the window" bare —
                      the secondary window's width, count, and whether it is ever
                      DRAWN are all still open.
```

### beta_share  [utility] [method]

⚠ **THE ALLOCATION RULE IS DECISION J (2026-10-10): `w_c = T_c / T_total`,
by VARIANT COUNT.** Mean `beta_v` per position is then `w_c / T_c = 1 / T_total` on
every chromosome, so this rule is what "normalise `beta_v` over the whole genome in
one go" MEANS — the per-chromosome share is bookkeeping, not a second weighting.
The rejected alternative, allocation by genetic length, would make a position's
weight depend on which chromosome it sits on.
⚠ `T_c` IS DEFINED ON THE FULL 88-HAPLOTYPE COHORT, fixed before any target is
selected. Under Arm 1's leave-one-out, up to 5 chr21 positions of 305,886 (0.0016%)
drop to `support == 0` and would otherwise make `w_c` weakly target-dependent —
a standing-constraint-1 problem avoided by construction rather than by smallness.
⚠ The shares therefore cannot be stated as 22 numbers until preprocessing has run,
because `T_c` is an output of it. The pre-registration pins the RULE; a rule is
re-derivable and falsifiable where 22 bare numbers are not. UNMEASURED: `T_total`.
```
Is:                   The fraction of the ONE genome-wide beta_v budget that a single
                      autosome carries: sum of beta_v over that chromosome's positions,
                      with the 22 shares summing to 1 (decision B). It is what replaces
                      the per-chromosome "sums to 1" the code enforces today, and it is
                      why a chr21-normalised figure is not a genome-wide figure: on the
                      genome-wide axis chr21's weight is its SHARE, not 1.
Computed by:          NOT YET IMPLEMENTED. It needs a genome-wide normaliser over the
                      22 autosomes; src/mechanism/utility.py:beta_v normalises over the
                      positions it is handed and asserts the sum is 1
                      (src/mechanism/utility.py:171), so it cannot emit a share.
CHANGES MEANING WITH: the weighting scheme, through w_v = 1/support(v): the share is
                      NOT the chromosome's position count over the total, because
                      support varies along the genome. Also with `chain_span` under
                      decision C, since declared-rate regions add positions that carry
                      weight. ⚠ Any share computed on a subset of the 22 autosomes is
                      not this quantity.
Valid vs:             the other 21 shares from the SAME genome-wide normaliser, and
                      their sum against 1.
Status:               PRIMARY — the endpoint under decision A is the genome-wide
                      utility, and it cannot be read without these.
Aliases:              "the per-chromosome budget share". ⚠ UNMEASURED for every
                      autosome, chr21 included. Never estimate it from T alone.
```

### psi_t  [implementation] [mechanism]
```
Is:                   The local potential psi_t(j) = exp(eta_tau * beta_t * phi_t(p_t, h_t,j)).
                      The ONLY place the private target enters the sampler. The cohort
                      HMM (rho, A) stays fixed.
Computed by:          src/mechanism/utility.py:tilt_potentials (built 2026-09-30), which
                      also multiplies in the allowed-state mask so that applying one and
                      forgetting the other is not expressible.
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
Computed by:          never materialized explicitly; sampled from via ffbs. Materialized
                      ONLY on the toy, by src/mechanism/exact.py:release_distribution,
                      which is where the bound is checked.
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
Computed by:          src/mechanism/ffbs.py:log_z_p (built 2026-09-30) = the tilted
                      forward pass's sum_t log c_t MINUS baseline_log_norm.
                      ⚠ THE SUBTRACTION IS NOT OPTIONAL — see baseline_log_norm.
CHANGES MEANING WITH: tau. The bound scales with eta_tau, so a log_z_p quoted
                      without its tau cannot be checked.
Valid vs:             its own theoretical bracket. This is the cheapest correctness
                      assertion in the whole implementation: if log Z_p leaves
                      [0, eta_tau], the sampler is wrong. Assert it, do not report it.
Status:               DIAGNOSTIC
Aliases:              "the normalizer", "Z_p" (Z_p is the unlogged form)
```


### baseline_log_norm  [implementation] [diagnostic]
```
Is:                   The log conditioning constant of R_D: sum_t log c_t from an
                      UNTILTED forward pass, whose only potentials are the
                      allowed-state mask. It exists because missing_policy =
                      RENORMALISE enters as a 0/1 potential rather than as a
                      renormalised transition row, so R_D is the Li-Stephens chain
                      CONDITIONED on never occupying a donor with no called allele,
                      and that conditioning has its own normaliser.
Computed by:          src/mechanism/ffbs.py:baseline_log_norm (built 2026-09-30)
CHANGES MEANING WITH: the cohort and the chain span — and NOTHING ELSE. It is
                      target-independent by construction (its potentials are the mask
                      alone), which is what lets it be computed once per cohort and
                      reused across every target and every tau. Measured 2026-09-30 on
                      our chr21 build, K = 88, T = 305,886: **-584.8135**, in 1.9 s.
Valid vs:             itself, as a cached constant. Asserted equal to the tau = 0
                      pass through the target-facing code path.
Status:               DIAGNOSTIC — but load-bearing. ⚠ log_z_p IS THE TILTED PASS
                      MINUS THIS. Forgetting the subtraction puts log_z_p at about
                      -584 against a bracket of [0, 0.549], and the tempting "fix" of
                      widening the bracket would discard the cheapest correctness
                      check in the project.
Aliases:              none. Never "Z_base" in a table; it is not a reportable number.
```

### delta_u  [theory] [calibration] [gate]
```
Is:                   the GLOBAL RANGE of the utility,
                      sup over all admissible p, q, y of |u(p, y) - u(q, y)|. For a
                      non-negative indicator utility it is exactly `sum_v beta_v`.
                      ⚠⚠ **IT IS THE ONLY PROPERTY OF THE UTILITY THAT THE PRIVACY
                      BOUND DEPENDS ON.** Established 2026-10-12: `beta_v` appears
                      nowhere in Theorem 1's proof, only in the example construction
                      and in the factorisation, which is about COMPUTABILITY. So the
                      DISTRIBUTION of `beta_v` across positions is free -- concentrate
                      the whole budget on 10 positions or spread it over 22,134,562 and
                      `tau`, `eta_tau`, `epsilon_tau` and `p_succ_bound` are IDENTICAL,
                      provided `delta_u` is unchanged. This is what makes budget
                      reallocation privacy-neutral.
Computed by:          nothing computes it as a named value -- it is 1 BY CONSTRUCTION,
                      enforced as `sum beta_v` == the declared share inside
                      src/mechanism/utility.py:beta_v, and as `sum beta_v` <= 1 inside
                      src/mechanism/utility.py:tilt_potentials (added 2026-10-12, the
                      only check on the path every tilted release must take).
CHANGES MEANING WITH: the utility FUNCTION, and nothing else. Not the target, not the
                      cohort, not `tau`. ⚠ If a raw utility has `delta_u` != 1 the
                      calibration is `eta_tau = arctanh(tau) / delta_u` -- standing
                      constraint 2. An unnormalised utility does not weaken the bound,
                      it REMOVES it.
Valid vs:             nothing; it is a property of the utility definition, not a
                      measurement. Quote it only to state which calibration is in force.
Status:               STRUCTURAL. 1 by construction on every release on record.
                      ⚠⚠ **DO NOT SHRINK IT USING THE DATA.** `u_ceiling` is the range
                      achievable for a GIVEN TARGET and is a readout of that target's
                      missingness, so calibrating `eta_tau` to it is local-sensitivity
                      calibration, which Nissim-Raskhodnikova-Smith show is not private
                      at any `epsilon_tau`. And the SAFE form, dividing by
                      sup over p of `u_ceiling`(p), gains precisely nothing, because
                      that supremum is exactly 1 -- attained by any admissible path
                      called at every chain position. Ruled on 2026-10-12.
                      ⚠ Measured consequence of getting it wrong: at a stated
                      `tau` = 0.5, a `delta_u` of 1.979595 means a TRUE total variation
                      of 0.795929 and a true accuracy of 0.897964.
Aliases:              "the sensitivity", "the range of u". Never "the ceiling", which
                      is `u_ceiling` and is a different, target-dependent quantity.
```

### e_differ  [utility] [evaluation] [gate]
```
Is:                   `E[# released cells differing from the target]`, EXACT, over the
                      target's CALLED cells only. The raw COUNT whose rate is
                      `target_fidelity`: target_fidelity = e_differ / called_cells.
                      ✅ **THE PRE-REGISTERED PRIMARY ENDPOINT** is this quantity via
                      the MARGINALS route (prereg_2026-10-10.md section 3), reported
                      against its own arm's ceiling and never against 1.
                      ⚠ It is a DIFFER count, so LOWER is more faithful -- the same
                      direction as `target_fidelity`, whose entry records the 2026-10-08
                      inversion that this direction corrects.
Computed by:          src/mechanism/ffbs.py:expected_mismatches from the state marginals,
                      summed over chromosomes and over both haplotypes. Since 2026-10-12
                      it rides along inside `ffbs.release` via `exact_out`, computed from
                      the SAME `psi` and `forward` the draw already holds, which is
                      pinned inert with respect to the draw by tests/test_exact_out.py.
                      ⚠⚠ **THE FINITE-DIFFERENCE ROUTE IS FORBIDDEN** at genome-wide
                      scope: it returns a NEGATIVE `Var_Q(u)` there. The marginals route
                      needs no step size and is the only admissible one.
CHANGES MEANING WITH: the ARM (K = 86 under leave-one-out, K = 88 in arm 2), `tau`, the
                      utility function, and the CALLED-CELL DENOMINATOR, which is per
                      target. Report the count and its denominator together; a count
                      without its denominator is uninterpretable, not merely imprecise.
Valid vs:             the same target at another `tau`, same arm. ⚠ Across targets only
                      as a distribution, never pooled: leave-one-out gives each arm 1
                      target its own chain and its own denominator.
Status:               PRIMARY. **Measured ARM 1, 2026-10-12: median 7,630,132.36 cells
                      of 43,958,219 called at `tau` = 0, moving -4.3056 cells across
                      `tau` = 0 to 0.99** (n = 44 held-out individuals, K = 86, 22
                      autosomes, exact). 44/44 targets move in the faithful direction.
                      ⚠⚠ **THE MOVEMENT IS THE CLAIM, NOT THE LEVEL.** The level is
                      reached at `tau` = 0 by an untilted draw carrying ZERO target
                      information, so it belongs to human allele sharing and not to the
                      mechanism. Four cells in seven and a half million is what the whole
                      privacy range buys, and it reproduces arm 2's +8.345e-08 `e_q_u`
                      movement on a different arm.
Aliases:              "E[differ]" in the pre-registration. Never "mismatches" bare,
                      which hides whether a count or a rate is meant.
```

### loo_id  [dataset] [provenance] [gate]
```
Is:                   The identifier of the individual HELD OUT of the donor panel for a
                      release -- the machine-checked half of `loo_cohort`. In ARM 1 it
                      is the target's own id, because the target IS a cohort member and
                      the guarantee holds only when both its haplotypes leave the panel.
                      In ARM 2 it is literally absent: the target is external to `G`,
                      `D`, every panel and the attack database, so nothing is held out.
Computed by:          stamped on every record by src/mechanism/ffbs.py:release and
                      enforced by src/releases.py:_check_arm, which REFUSES an arm 1
                      release with no `loo_id` and an arm 2 release that carries one.
                      ⚠ Both directions are refused because each is a different false
                      claim: the first would be a release drawn against a panel
                      containing its own target, the second an arm 1 release mislabelled.
CHANGES MEANING WITH: nothing -- it is an identifier. ⚠ But it is part of the CONFIG
                      STAMP that decides whether two results are comparable at all
                      (`docs/plan.md` Reporting contract): graph+cohort id, arm, `tau`,
                      utility id, `loo_id`.
Valid vs:             itself. Two releases with different `loo_id` were drawn against
                      DIFFERENT donor panels and therefore different `R_D`, so their
                      utilities are not comparable (constraint 8).
Status:               REQUIRED throughout ARM 1 and absent throughout ARM 2. Verified
                      present and equal to the target on all 308 arm 1 releases,
                      2026-10-12.
                      ⚠ It records the hold-out, which is what makes arm 1's guarantee
                      conditional rather than void. It does NOT record the other
                      condition: `G` was built from all HPRC assemblies INCLUDING the
                      target's, so arm 1's guarantee is conditional on the graph too.
                      Nothing on the record captures that; it must be said in prose.
Aliases:              none. Never "held out" bare in a table -- the field name is this.
```

### e_match_nonref  [utility] [evaluation]
```
Is:                   The EXPECTED number of the target's NON-REF positions at which
                      the released path carries the target's allele. A count, not a
                      weighted average: sum over positions with target allele index
                      > 0 of P(release matches there) under the posterior.
                      ✅ **THE ONLY UTILITY NUMBER THAT MAY BE COMPARED ACROSS
                      WEIGHTINGS**, which is why it exists (2026-10-24). e_q_u is
                      beta-weighted and therefore weighting-SPECIFIC, so comparing
                      two weightings' e_q_u is the constraint-7 error and gave an
                      artifactual 5.44x; e_differ counts every called cell and is
                      dominated by the ~19.6M REF positions a cohort mosaic matches
                      by chance, so it FAVOURS whichever weighting recovers common
                      alleles. This one is independent of beta entirely.
Computed by:          scripts/eta_sweep.py:one_eta, from the same posterior state
                      marginals as e_q_u and e_differ -- no extra forward pass, zero
                      Monte-Carlo error. Summed over the 22 autosomes by
                      scripts/eta_sweep.py:run_genome.
CHANGES MEANING WITH: the TARGET (its non-ref set is its own) and the cohort (which
                      alleles are callable at all). NOT with beta -- that is the
                      point. ⚠ Report the MOVEMENT from eta_tau = 0, never the level:
                      the untilted mosaic already matches 2,357,272 of 4,994,881
                      non-ref positions (47.19%, HG00096, genome-wide) carrying zero
                      target information.
Valid vs:             the same quantity under another weighting at MATCHED
                      acc_ceiling. ⚠ Never at matched eta_tau -- two weightings at
                      one eta_tau are at different privacy levels.
Status:               PRIMARY for choosing a weighting.
Aliases:              none. ⚠ NOT "informativeness": it counts non-ref alleles
                      equally, and an allele 80 of 88 donors carry is nearly as
                      uninformative as a reference one. It is a better comparison
                      metric than e_q_u or e_differ, not a biological utility
                      measurement.
```

### target_alt  [utility] [mechanism]
```
Is:                   A beta_v weighting: UNIFORM over the positions where THAT
                      haplotype carries a non-ref allele, zero elsewhere,
                      normalised by the GENOME-WIDE count of them per haplotype.
                      ✅ **THE TOOL'S WEIGHTING (Dylan, 2026-10-24.)** Recovers
                      1.633x the non-ref alleles of the production weighting at
                      matched acc_ceiling, altering 0.613x the cells -- 2.67x more
                      efficient per altered cell -- with 98.1% of what it recovers
                      being non-ref against production's 36.8%.
Computed by:          scripts/eta_sweep.py:build_betas with
                      scripts/exact_utility.py:alt_totals for the normaliser; the
                      same expression as scripts/alt_beta_sweep.py, one definition.
                      HG00096: 2,493,299 and 2,501,582 positions per haplotype, so
                      beta_v = 4.0108e-07 there.
CHANGES MEANING WITH: the TARGET. ⚠ It is TARGET-DEPENDENT, unlike w_v = 1/support(v).
                      Theorem 1 still holds -- beta >= 0, sums to 1, u stays in
                      [0,1] -- but it is a DIFFERENT utility function from the
                      production one and constraint 7 forbids them sharing a column.
                      Per-haplotype: the two haplotypes' non-ref sets differ, so the
                      two chains get DIFFERENT beta arrays.
Valid vs:             production and target_rarity, via e_match_nonref at matched
                      acc_ceiling only.
Status:               PRIMARY -- the tool's default.
Aliases:              "the non-ref weighting". ⚠ Never "the informative weighting":
                      it does not distinguish a singleton from an allele 80 of 88
                      donors carry.
```

### target_rarity  [utility] [evaluation]
```
Is:                   A beta_v weighting: w_v = -log f_v, with f_v the COHORT
                      frequency of the target's own allele (carriers / called over
                      the donor panel). Zero where the target has no call, zero
                      where no donor carries the allele (unmatchable at any eta_tau,
                      so budget there is wasted), and zero where every donor carries
                      it (-log 1 = 0). Bounded by log(K) = 4.4773 at K = 88, which
                      is what keeps u inside [0,1] for standing constraint 2.
Computed by:          scripts/exact_utility.py:rarity_weights, normalised by
                      rarity_totals over all 22 autosomes;
                      scripts/eta_sweep.py:build_betas.
CHANGES MEANING WITH: the TARGET and the COHORT (f_v is a cohort frequency, so the
                      same target against a different panel gives different weights).
Valid vs:             production and target_alt, via e_match_nonref at matched
                      acc_ceiling only.
Status:               ⚠ TRIED AND REJECTED 2026-10-24, retained as the comparison.
                      It recovers only 1.060x production and **0.649x target_alt**.
                      Better TARGETED than a ref/non-ref cut and still worse, because
                      it is less CONCENTRATED: effective support 48,700 positions
                      against target_alt's 33,303 on chr21, var_q_u 4.7268e-07
                      against 9.3736e-07, and the exchange rate goes as
                      sqrt(var_q_u). **Concentration beats targeting here.**
Aliases:              "the rarity weighting", "-log f weighting".
```

### e_q_u  [utility] [theory] [evaluation]
```
Is:                   E_Q[u], the EXACT expected utility of a release under Q^tau_p —
                      the first derivative of the log normaliser,
                      E_Q[u] = d log Z_p / d eta. No sampling and no Monte-Carlo error.
                      ✅ **THE PRIMARY PRE-REGISTERED ENDPOINT (decision A, 2026-10-08)
                      IS THIS QUANTITY, GENOME-WIDE**, reported exactly. The sampled
                      average (`utility_retained`) is its DIAGNOSTIC counterpart, not
                      the endpoint. ⚠ The contiguous-window version that this line
                      once called a labelled SECONDARY arm is RETIRED (2026-10-24);
                      the live comparison is production vs `target_alt` weighting.
Computed by:          scripts/results_table.py:var_and_mean — three forward passes of
                      src/mechanism/ffbs.py:log_z_p at eta - h, eta, eta + h, central
                      difference, with the step PINNED at H_STEP = 0.05
                      (scripts/results_table.py:39). Cross-checked against direct
                      enumeration on the toy. The marginals route —
                      src/mechanism/ffbs.py:state_marginals, exact for any quantity
                      that is a per-position sum — is the alternative and is what
                      src/mechanism/ffbs.py:expected_mismatches already uses.
                      ✅ **AND THE GENOME-WIDE ROUTINE NOW EXISTS (2026-10-12):**
                      scripts/exact_utility.py assembles it across the 22 autosomes by
                      the marginals route under a genome-wide beta, one forward +
                      state_marginals pass per (chromosome, haplotype, tau), 44 passes
                      in ~30 min. This paragraph used to read "NOT YET A GENOME-WIDE
                      ROUTINE". It carries NO finite-difference step: the marginals
                      route needs no h.
CHANGES MEANING WITH: the ESTIMATOR ROUTE and, for the finite-difference route, the
                      step h — a derivative taken at a different h is a different
                      number, so h must be stamped. And with `beta_v`'s scope and
                      normalisation axis: the recorded value **E_Q[u] = 0.6951**
                      (HG00096, chr21, scope "chain", ARM 2, tau = 0.5) was taken with
                      beta CHR21-NORMALISED, so it is NOT the genome-wide endpoint
                      decision A pre-registers. ✅ **THE GENOME-WIDE VALUE IS NOW
                      MEASURED (2026-10-12): E_Q[u] = 0.708967646 at tau = 0 and
                      0.708967729 at tau = 0.99** — HG00096, ARM 2, 22 autosomes,
                      T = 22,134,562, shares by decision J summing to 1.000000, scope
                      "chain", n = 1 target, exact marginals route, zero Monte-Carlo
                      error. **The MOVEMENT is +8.345e-08 and that is the figure**; the
                      LEVEL is reached at tau = 0 too, where the release carries zero
                      target information, so the level is not attributable to the
                      mechanism.
                      ⚠ AND WITH THE ARM, through the ceiling it must be read against:
                      quote it beside `u_ceiling` for its OWN arm, never against 1.
Valid vs:             the same route at matched tau, beta scope and normalisation axis;
                      and against its own arm's `u_ceiling` and the tau = 0 floor.
                      On the toy, against direct enumeration — that equivalence is what
                      licenses the derivative route at scale.
Status:               PRIMARY — the pre-registered endpoint under decision A.
Aliases:              "E_Q[u]", "the exact expected utility". ⚠ Never "the utility"
                      bare, and never interchangeable with `utility_retained`, which is
                      a finite-sample average of the same quantity.
```


### e_q_u_over_ceiling  [utility] [evaluation] [derived]
```
Is:                   `e_q_u` / `u_ceiling` — the EXACT expected utility as a fraction
                      of the best any release could achieve for that target. A pure
                      ratio of two exact quantities, so it carries no sampling error.
                      ⚠ **IT IS NOT `utility_retained`**, which is the SAMPLED u_path
                      average and is a LEVEL, not a ratio. This term exists because
                      scripts/exact_utility.py emitted the ratio under that name on
                      2026-10-12 and the label was wrong; the saved JSON key was
                      relabelled, values bit-identical.
Computed by:          scripts/exact_utility.py — a derived row, from the same two
                      exact quantities on the rows above it (built 2026-10-12).
CHANGES MEANING WITH: everything `e_q_u` and `u_ceiling` each change with — the arm,
                      tau, beta_v's scope AND its share, and the target. Both the
                      numerator and the denominator are arm-indexed, so a ratio built
                      across arms is meaningless.
Valid vs:             the SAME target's ratio at another tau, same arm, same share.
                      ⚠ Never across targets: `u_ceiling` differs per target, so the
                      ratio's denominator moves and the comparison is not like-for-like.
Status:               DERIVED and SECONDARY. The primary endpoint is `e_q_u` itself,
                      and the figure that carries a claim is its MOVEMENT in tau, not
                      any level. **Measured 0.858617102 at tau = 0 and 0.858617203 at
                      tau = 0.99** — HG00096, ARM 2, 22 autosomes, share 1.000000,
                      scope "chain", n = 1 target, exact.
                      ⚠⚠ **DO NOT QUOTE THIS AS "WE RETAIN 86% OF UTILITY".** 0.8586 is
                      reached at tau = 0, where the release is a pure cohort mosaic
                      carrying ZERO target information — any two humans share most of
                      their alleles. Quoting the level is exactly the free-lunch error
                      `target_fidelity` was defined to prevent.
Aliases:              none. Never "utility retained", which is a different quantity.
```
### var_q_u  [utility] [theory] [diagnostic]
```
Is:                   Var_Q(u), the EXACT variance of the utility under Q^tau_p — the
                      second derivative of the log normaliser,
                      Var_Q(u) = d2 log Z_p / d eta2. It is what says whether a
                      measurement of `e_q_u` is resolvable at a feasible number of
                      draws, which is why it is reported beside it rather than hidden.
Computed by:          ✅ **scripts/exact_variance.py:haplotype_variance, EXACTLY and in
                      O(TK)** (2026-10-22). One centred forward sweep carrying the first
                      AND second conditional moments of the running prefix sum over
                      FFBS's backward chain; no sampling, no finite difference, no pair
                      marginals. Gated against brute-force enumeration on the toy
                      (7.7e-15), against 30 real chr21 draws (SD ratio 1.0334), and by
                      its own centred first moment, which must be 0 and comes back at
                      2.4e-13 or below.
                      ⚠⚠ **IT HAS A FLOOR NEAR 1e-33 AND IT ERRS IN THE FLATTERING
                      DIRECTION — measured 2026-10-23.** The centring subtracts
                      quantities of order 1, so once the posterior concentrates onto
                      one path the result is cancellation noise: on the toy it returns
                      a CONSTANT 7.703720e-34 where enumeration gives 3.52e-127 and
                      4.42e-276, overstated by 10^93 and 10^242. Var is the SLOPE of
                      the response curve, so an overstated floor overstates how much
                      the next unit of `eta` would buy. **NEVER quote var_q_u below
                      ~1e-28.** `e_q_u` is unaffected (no centring, no cancellation)
                      and agrees to 4.76e-16 at every `eta`. See docs/bugs.md.
                      ⚠ THE SUPERSEDED ROUTE: scripts/results_table.py:var_and_mean,
                      a second central difference with H_STEP = 0.05. It is the route the
                      pre-registration FORBIDS at genome scope — it returns a NEGATIVE
                      variance there — and it is what left this quantity unmeasured
                      genome-wide until the sweep above existed.
                      ⚠ **[CORRECTION - 2026-10-22]** `docs/plan.md` and
                      `scripts/exact_utility.py` both say settling the 43.6x "needs
                      Var_Q(u) from PAIR marginals". Wrong about the METHOD, not the
                      quantity: all-pairs is O(T^2 K^2) and unreachable at
                      T = 22,134,562, while `u` is an ADDITIVE functional of the path, so
                      one O(TK) sweep is exact. The stale wording is still in that
                      script's docstring and in paper/draft_2026-10-13.md.
CHANGES MEANING WITH: the BUILD, the arm, the ploidy and `beta_v`'s scope and share, and
                      every figure needs all five. Diploid is the release's own quantity,
                      Var(u_dip) = (Var(u1) + Var(u2))/4 over two independent chains at
                      eta_tau/2.
                      **CURRENT BUILD, exact** — ARM 2, HG00096, 22 autosomes at
                      decision-J shares, diploid: **3.153126468e-08**
                      (`logs/B_exact_variance_genome.md`). Same target, chr21 at
                      share 1: **1.374101750e-06** at tau = 0 and 1.374107483e-06 at
                      tau = 0.99 (`logs/B_exact_variance_chr21_share1.md`). ARM 1, HG00438 held out
                      at K = 86, chr21 at its own share 0.0150248984, diploid:
                      **4.652085850e-10** (`logs/B_exact_variance_arm1_chr21.md`).
                      ⚠ The recorded **4.046e-06** / sd(u) = 2.011e-03 (HG00096, chr21,
                      scope "chain", ARM 2, tau = 0.5) is PRE-REBUILD and HAPLOID, by the
                      finite-difference route. The matched configuration on the current
                      chain gives 2.790926e-06, and `e_q_u` moved 0.6951 -> 0.5898 over
                      the same rebuild, so the gap is the BUILD and not the estimator.
                      Neither figure is wrong; they describe different chains.
                      ⚠ 94.8% of Var(u) is between-position COVARIANCE from the mosaic's
                      segment structure, which no reweighting of beta_v touches — see
                      `cov_share`, and `beta_v`'s last block before concluding that a
                      utility can be given variance by reweighting. ✅ That share is now
                      reproduced INDEPENDENTLY: 0.948906 by the exact route on the current
                      build against 0.948 by finite differences on the old one.
Valid vs:             itself across tau at matched beta scope, arm and build; against the
                      per-draw sd it implies; and against the `e_q_u` MOVEMENT, which it
                      predicts — `d E_Q[u] / d eta = Var_Q(u)`, so the movement across the
                      grid is the variance integrated. ✅ MEASURED: eta_tau * Var(tau = 0)
                      = 3.636770e-06 against a measured movement of 3.636777e-06 on chr21
                      at share 1 (relative 1e-05), and 8.345230e-08 against a logged
                      8.345e-08 genome-wide (relative 3e-05). BOTH scopes.
                      ✅ **THE 43.6x BETWEEN THE TWO SCOPES IS A PRODUCT OF THREE
                      MEASURED FACTORS**: sum beta_v^2 128.6181 / 2.2679 (the
                      beta_v^2-weighted mean p(1-p) rises, 0.007675 -> 0.017405) =
                      56.7130 diagonal part / 1.3014 (the covariance multiplier rises,
                      18.27 -> 24.08) = **43.5790**, reproducing the measured ratio to
                      1.6e-16. ⚠⚠ The standing correlation hypothesis is only the 1.30x;
                      the larger 2.27x is WHICH positions carry beta_v weight at each
                      scope, which no correlation argument touches.
                      ⚠ It is the denominator of the project's own feasibility argument:
                      quoting a sampled utility difference without it is how an
                      unresolvable measurement gets reported as a result.
Status:               DIAGNOSTIC — reported beside `e_q_u`, never alone. ⚠ It is NOT a
                      pre-registered endpoint; it explains and predicts the primary
                      endpoint's movement and is reported as explanation, not as a result
                      of its own.
Aliases:              "Var_Q(u)", "sd(u)" for its square root. Never "the noise".
```

### var_q_u_diag  [utility] [theory] [diagnostic]
```
Is:                   what `var_q_u` WOULD BE if the per-position allele matches were
                      independent: sum_v beta_v^2 * p_v * (1 - p_v), with p_v the
                      posterior probability that the release carries the target's allele
                      at v. It is the quantity the project's `sum(beta_v^2)` argument
                      implicitly assumes, written out so the assumption is measured
                      rather than reasoned about.
Computed by:          scripts/exact_variance.py:haplotype_variance, from the SAME
                      posterior marginals as `var_q_u` — so the comparison between them
                      isolates the independence assumption and nothing else.
CHANGES MEANING WITH: exactly what `var_q_u` does. Never quote one without the other.
Valid vs:             `var_q_u` at matched everything, and only that. ⚠ It is NOT an
                      estimate of the variance and must never be reported as one: on
                      chr21 at share 1 it is **7.130104363e-08** against the true
                      1.374101750e-06, low by a factor of 19.3. It is the null it exists
                      to reject.
Status:               DIAGNOSTIC — the denominator of `cov_share`.
Aliases:              "the independence prediction". Never "the diagonal variance"
                      alone, which reads as a property of the data rather than of an
                      assumption about it.
```

### cov_share  [utility] [theory] [diagnostic]
```
Is:                   `1 - var_q_u_diag / var_q_u`: the share of the utility's variance
                      that comes from COVARIANCE between positions rather than from the
                      positions individually. It is why the tilt has more to push on than
                      a per-position account predicts, and it is the mechanism behind the
                      43.6x that two first-order predictions both missed.
Computed by:          scripts/exact_variance.py, as the ratio of the two above.
CHANGES MEANING WITH: the arm and the beta scope, but FAR LESS than the variance itself,
                      which is what makes it the quotable quantity: measured **0.948111**
                      (ARM 2, chr21, share 1, diploid, current build), **0.948906** (the
                      same by the pre-rebuild route's matched configuration) and
                      **0.957022** (ARM 1, HG00438, chr21 at its own share) and
                      **0.960128** (ARM 2, GENOME-WIDE, 22 autosomes). Per chromosome it
                      spans 0.9334 (chr22) to 0.9684 (chr4). Two arms, two builds, two
                      estimators, every scope: **0.93-0.97**.
                      ⚠ It is a LEVEL, and the two-scope 43.6x is NOT this number — see
                      `var_q_u`'s three-factor decomposition, where the covariance
                      multiplier contributes only 1.30x of 2.95x.
Valid vs:             itself across arms, scopes and builds — the comparison this quantity
                      is FOR, since it is a ratio and the scale cancels.
Status:               DIAGNOSTIC. ⚠ Not pre-registered; it is the explanation of a
                      pre-registered endpoint's behaviour.
Aliases:              "the covariance share". ⚠ Never "the correlation", which is a
                      different quantity, and never phrase it as a share of the UTILITY —
                      it is a share of the utility's VARIANCE.
```

### var_q_differ  [utility] [evaluation] [diagnostic]
```
Is:                   Var_Q(`e_differ`) — the exact posterior variance of the NUMBER of
                      released cells differing from the target, summed over both
                      haplotypes. Its square root, `sd_q_differ`, is the draw-to-draw
                      standard deviation of a single release, so it is what says whether
                      a difference between two single-draw table columns means anything.
Computed by:          scripts/exact_variance.py:haplotype_variance with unit weights and
                      the differs-from-target indicator — the same sweep as `var_q_u`, so
                      one gate covers both. ⚠ Diploid adds OUTRIGHT here
                      (Var(h1) + Var(h2)), not at a quarter: the cell count is a SUM over
                      haplotypes while `u_dip` is their MEAN. Using the wrong one reports
                      an SD 2x off for the row it is meant to describe.
CHANGES MEANING WITH: the arm, the target, the chromosome scope and the denominator — and
                      it does NOT change meaningfully with `tau`, which is the point.
                      Measured: **6.467224125e+07** genome-wide (ARM 2, HG00096, 22
                      autosomes), 7.739900306e+05 on chr21 at share 1, and
                      1.127535146e+06 in ARM 1 (HG00438 held out, chr21 at its own
                      share), all flat across tau = 0 to 0.99.
Valid vs:             the observed scatter of single-draw columns at matched arm and
                      scope. ✅ MEASURED: the head-to-head's seven `tau` columns give an
                      empirical SD of 1,312 cells for `cells changed` against the exact
                      1,061.9 — ratio 1.24, where a 7-draw SD estimate carries 1 sigma of
                      +-29%. So that row's spread IS draw noise.
                      ⚠ It does NOT cover `utility_loss_panmixer`, whose weights are
                      PanMixer's own cost vector and are not loaded by this script. The
                      `cells changed` and `target_fidelity` rows are covered; the cost row
                      is not.
Status:               DIAGNOSTIC — report it beside any single-draw column, because a
                      column quoted without it invites a reader to compare noise.
Aliases:              none. Its square root has its own entry, `sd_q_differ`, because it
                      carries a different unit. Never "the error bar", which implies an
                      estimate; this is exact.
```

### sd_q_differ  [utility] [evaluation] [diagnostic]
```
Is:                   sqrt(`var_q_differ`), in CELLS — the draw-to-draw standard
                      deviation of one release's differing-cell count. It exists as its
                      own name because it is the form every table needs: the unit of the
                      variance is cells SQUARED, which no reported column is in.
Computed by:          scripts/exact_variance.py, as the square root of `var_q_differ`.
CHANGES MEANING WITH: everything `var_q_differ` does. Measured **8,041.9 cells** (ARM 2,
                      HG00096, GENOME-WIDE), **879.8 cells** (same target, chr21 at
                      share 1) and **1,061.9 cells** (ARM 1, HG00438 held out, chr21 at
                      its own share), all flat across tau.
Valid vs:             the observed spread of single-draw columns at matched arm and
                      scope, and the difference between any two such columns — a gap
                      smaller than about 2 of these is not a finding. ⚠ Do NOT compare it
                      across arms: ARM 1's panel, chain and utility function all differ.
Status:               DIAGNOSTIC — quote it beside any single-draw column.
Aliases:              none. Never "the error bar" or "the noise".
```

### ratio_spread  [privacy] [evaluation] [gate]
```
Is:                   max_y (Q_p(y)/Q_q(y)) divided by min_y of the same ratio, over
                      the output support, for a pair of admissible inputs p, q.
                      Theory pins it: the ratio is
                      exp(eta_tau (u_p(y) - u_q(y))) * Z_q/Z_p, and the Z factor is
                      common to every y, so it CANCELS OUT of the spread -- leaving
                      exp(eta_tau * range_y(u_p - u_q)), which is at most
                      exp(2 eta_tau) = exp(epsilon_tau) because each utility is in
                      [0,1].
Computed by:          src/mechanism/exact.py:ratio_spread, on the toy support only
CHANGES MEANING WITH: the PAIR. The bound is over all pairs, so one favourable pair
                      proves nothing; report the worst found. Nothing else -- there is
                      no estimator here and no support-size bias, which is exactly
                      what distinguishes it from tv_empirical.
Valid vs:             exp(epsilon_tau), at matched tau, on a stated support.
Status:               PRIMARY — mark * in tables. ⚠ IT IS THE ONLY GATE WITH POWER
                      against a miscalibrated tilt. Measured 2026-09-30 on the toy: a
                      constructed input pair SATURATES it exactly at every tau tested,
                      while the worst exact TV reaches only 27.4-43.1% of tau. So a
                      tilt applied with epsilon_tau in place of eta_tau -- twice the
                      intended exponent, a one-token error, the highest-severity error
                      class in the project -- PASSES the TV gate at every tau and
                      passes the log_z_p bracket below tau = 0.5, and fails only this
                      one. See `docs/memory.md` 2026-09-30.
Aliases:              "the pointwise bound". Never "the likelihood ratio" bare -- that
                      is the per-y quantity; this is its range over the support.
```

### ffbs  [implementation] [mechanism]
```
Is:                   Forward-filtering / backward-sampling: the exact sampler for
                      Q^tau_p. Forward pass computes per-position normalized alpha_t;
                      backward pass draws Z_T ~ Categorical(alpha_T) then walks back
                      sampling Z_t | Z_t+1. Algorithm 1 of the proposal.
Computed by:          src/mechanism/ffbs.py (built 2026-09-30). ✅ VERIFIED against exact
                      enumeration 2026-09-30: output distribution and state-path JOINT
                      both within Monte-Carlo error at K=4, T=7.
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
Computed by:          src/mechanism/hmm.py:switch_mass_from_cm (built 2026-09-30).
                      `switch_mass` is s_t, the coefficient of J in A = (1-s)I + (s/K)J,
                      i.e. the TOTAL mass; the per-state value is s_t/K and the chance of
                      leaving a state is s_t(K-1)/K. The decomposition is cross-checked
                      against the published closed form at construction.
CHANGES MEANING WITH: whether b_t is the per-target or the total switch mass.
                      Getting that wrong rescales every transition by K.
Valid vs:             a dense-transition HMM giving identical forward values. That
                      equivalence is the gate on the O(TK) optimization.
Status:               PRIMARY
Aliases:              "the PanMixer transitions"
```

### T_blocks  [implementation] [dataset]
```
✅ SETTLED 2026-09-28 (was open question 8): **one chain position per VARIANT** —
   model C. Model B (one position per LD block) is the retained fallback.
Is:                   T, the number of ordered positions in the hidden Markov chain —
                      the steps the sampler takes along the chromosome, at each of
                      which it picks a donor haplotype and emits that donor's allele.
Computed by:          src/cohort/genmap.py:chain_mask (built 2026-09-29) plus
                      src/mechanism/hmm.py:from_arrays, which applies the mask.
                      One position per VCF record inside `chain_span`.
                      chr21: 305,886 positions (340,849 records less the 34,371 below
                      the coordinate cut).
CHANGES MEANING WITH: the model, AND the `chain_span` policy. ⚠ Every T on record was
                      measured under the SUPERSEDED coordinate cut. Under decision C the
                      map-silent positions are kept and sampled under `declared_rate`,
                      so T GROWS; the post-decision-C T is UNMEASURED on every autosome,
                      chr21 included. The figures below are stamped to the old cut.
                      C = 340,849 raw / 305,886 after the cut. ✅ MODEL B IS
                      BUILT (2026-10-07) and on OUR blocks gives **59,472** chain units,
                      not the 100,757 quoted from PanMixer's block set — our blocks come
                      from the 30x panel with the 10 test targets removed, SNPs only, and
                      orphans stay singletons. Runtime and numerical headroom both scale
                      with T — see beta_v for why finer T costs float32 headroom but NOT
                      privacy budget.
                      ⚠ UNDER MODEL B A UNIT IS ATOMIC: one donor emits its whole span, so
                      `allowed` is an AND across the unit and the transition is the GAP
                      between units (last position of one to first of the next), which is
                      what discards the within-block recombination. Measured: expected
                      switches 3,394.1 against C's 7,021.0, and the longest verbatim donor
                      run is 1,212 variants — departure (b), which `tau` does not bound.
                      4 units have no donor callable across them and are dropped, costing
                      1,222 chain positions — and the reason is exponential: the atomic
                      requirement needs ONE donor called at EVERY position, which at 4%
                      missingness is about 0.96^n for a unit of n variants, so only the
                      very largest units fail (534, 257, 246 and 185 variants).
Valid vs:             the same model and the same chain_span.
Status:               SECONDARY
Aliases:              "T", "sites". ⚠ The name `T_blocks` is now a misnomer under C,
                      where a position is a variant and no blocks exist in the
                      mechanism. Retained so existing references resolve; write "T".
```

### n_sites  [implementation] [dataset]
```
Is:                   The length of the SITE AXIS — the number of VCF records for a
                      chromosome, and therefore the first dimension every preprocessing
                      artifact is indexed by. chr21: **340,849** (OUR build; PanMixer's is 340,824).
Computed by:          src/cohort/arrays.py (built 2026-09-29) — the record count of the
                      deconstructed per-chromosome VCF.
CHANGES MEANING WITH: the VCF. A re-sort, a re-filter or a different `vg` version changes
                      which row is which SILENTLY, because every join in the pipeline is
                      positional. ⚠ Distinct from T: `n_sites` is the axis, T is the
                      number of CHAIN positions, which is `n_sites` restricted to
                      `chain_span` (chr21: 305,886 of 340,849 under the SUPERSEDED cut;
                      under decision C the restriction is support(v) > 0 alone, so T
                      rises toward `n_sites` by an amount that is UNMEASURED).
                      ⚠ And with the BUILD: decision F rebuilds both arms' starting data
                      on one consistent GRCh38, and the 22 autosomes each have their own
                      `n_sites` — only chr21's is measured.
Valid vs:             the same VCF, verified by the `sites.tsv` digest.
Status:               SECONDARY
Aliases:              "the site axis", "rows". Never "T".
```

### allele_lengths  [implementation] [dataset]
```
Is:                   The length in base pairs of every declared allele at every record —
                      ragged, `num_alleles[i]` entries for row i. chr21: 767,376 declared
                      alleles over 340,849 records, max 90 at one site.
Computed by:          src/cohort/arrays.py (built 2026-09-29) — our VCF-to-matrix step,
                      the only moment REF/ALT are in hand.
CHANGES MEANING WITH: nothing, but its ABSENCE changes several things. Carried because
                      the integer allele matrix knows an allele's INDEX and not one base
                      of its SEQUENCE, hence not its length — so without this array no
                      downstream quantity can see how much sequence an edit rewrote.
                      PanMixer's per-record utility has exactly this blindness: 3.13 Mb
                      was rewritten on one measured run and no metric in its suite could
                      see it. Needed for any length-weighted utility, and cheap:
                      **2.9 MB as int32 against 118.3 MB for the sequences themselves**,
                      whose distribution is extreme — 10 records hold 58.1% of chr21's
                      total ALT bytes.
Valid vs:             the same site axis.
Status:               SECONDARY — not consumed by the v1 mechanism; carried because it is
                      cheap NOW and unrecoverable later without re-reading the VCF.
Aliases:              none
```

### chain_span  [dataset] [method]
```
Is:                   The coordinate interval of a chromosome the chain actually runs
                      over, and — since decision C — the cM SOURCE of every position
                      inside it.
                      ✅ **THE COORDINATE CUT IS SUPERSEDED FOR THE MAP-SILENT REGIONS
                      (decision C, Dylan, 2026-10-08).** Regions where the genetic map
                      gives NO cM AT ALL are no longer EXCLUDED — they are SAMPLED
                      under a DECLARED CONSTANT RECOMBINATION RATE (`declared_rate`).
                      So the chain runs over the whole assembled autosome, the map's
                      marker span becomes an INPUT PROPERTY that selects which
                      positions get a MAP-INTERPOLATED cM and which get a declared-rate
                      cM, and nothing is dropped for want of a cM. This RESOLVES
                      deferred decision R2 (the provisional cut). Positions may still
                      leave the chain on named, NON-MAP grounds: support(v) == 0, where
                      no donor has a called allele, so the sampler has nothing to emit
                      and RENORMALISE has an empty pool.
                      ⚠ **SCOPE, NARROW AND DELIBERATE: decision C REACHES THE
                      MAP-SILENT REGIONS ONLY.** The INTERIOR zero-cM map PLATEAUS
                      inside the mapped span are NOT overridden — see the block below,
                      which still stands. That question is DEFERRED until a per-autosome
                      gap measurement exists.
                      ⚠ **THE RATE'S VALUE IS NOT YET CHOSEN** and whether it is global
                      or per-chromosome is open, so the bp and position counts the two
                      cM sources cover are UNMEASURED on every autosome, chr21
                      included (`map_coverage`). Do not quote a post-decision-C span
                      until they are measured.
                      THE SUPERSEDED chr21 CUT, retained because every number on record
                      was measured under it: chr21 (decided 2026-09-28, closed at both
                      ends 2026-09-29) ran **12,968,320 < POS <= 46,680,243** — the
                      interval over which the genetic map carries real markers —
                      keeping **305,886** of 340,849 variants (89.74%), with everything
                      outside excluded from the chain and from the release.
Computed by:          configs/inputs.chr21.yaml `chain_span`, applied by
                      src/cohort/genmap.py:chain_mask (built 2026-09-29). Stored as a
                      MASK over the full axis, never as a row filter.
                      ⚠ **BOTH STILL ENCODE THE SUPERSEDED CUT, and neither can express
                      a declared-rate position.** The config block carries `start`,
                      `end`, `keeps: 305886`, `dropped_below`, `dropped_above`,
                      `provisional: true` and two `unresolved` entries;
                      src/cohort/genmap.py:chain_mask RAISES on any position inside the
                      span with no interpolated cM. Not yet rewritten — the rate's value
                      is unchosen, so there is nothing to encode yet.
CHANGES MEANING WITH: the chromosome and the genetic map. It is a property of the
                      ASSEMBLY and the MAP, never of the target — which is what makes it
                      safe: it is fixed before the target is opened, so R_D stays a
                      function of D alone and Theorem 1 is untouched. ⚠ NEVER let an
                      exclusion become target-dependent (e.g. "drop positions where the
                      target is missing"); that would leak.
Valid vs:             the same span AND the same cM-source policy. Any utility or
                      fidelity number must name both, since the denominator changes
                      with the span and the transition model changes with the source.
Status:               PRIMARY — a release that does not state its span is not quotable,
                      and after decision C a span stated without its declared rate and
                      its per-position cM source is not a stated span.
Aliases:              "the coordinate cut" — ⚠ now a HISTORICAL name: under decision C
                      the map-silent regions are no longer cut. Do not use it for the
                      post-decision-C span.
```

**WHAT THE p ARM IS, and why chr21's WAS excluded** (the exclusion is superseded by
decision C; the p arm is now sampled under `declared_rate`, and everything below is the
recorded reason it needs a DECLARED rate rather than interpolation). Every chromosome has a pinch point,
the **centromere**, dividing it into a short arm — **p**, from French *petit* — and a long
arm, **q** (simply the next letter). On most chromosomes both arms carry genes. But five
human chromosomes — **13, 14, 15, 21 and 22** — are **acrocentric**: the centromere sits
very near one end, so the p arm is small and consists almost entirely of repetitive DNA plus
long tandem arrays of **rDNA** genes (the ones encoding ribosomal RNA). Those arrays are
near-identical *between* the five acrocentric chromosomes, which is why assemblers cannot
place them uniquely and why GRCh38 fills them with N. chr21's p arm is the 0–~12 Mb region.

**The live defect the cut fixed — and which decision C fixes only HALF of.**
⚠ **CORRECTED 2026-10-09.** This block briefly read that a declared rate "gives these
positions a NON-ZERO `Delta_x`, so `P(switch) > 0` and the verbatim-emission failure
below is answered by sampling rather than by exclusion." **That over-states decision C's
reach**, because the failure below has two halves and the rate reaches only one:
- **Half one — MAP-SILENT, and the rate DOES apply.** `np.interp` clamps by default, so
  all **26,044** variants below the map's start (10,326,676) receive the map's first
  value. ⚠ **MEASURED ON OUR BUILD 2026-10-10**, replacing the 26,018 that stood here —
  that figure is PanMixer's vg-1.36 340,824-record build. On our 340,849-record build it
  is 26,044, of which 26,043 have `support > 0`. This is the quantity the declared rate
  actually applies to on chr21.
- **Half two — A MEASURED ZERO INSIDE THE MAPPED SPAN, and the rate does NOT apply.** The
  map's own first interval (10,326,676 → 12,968,320, spanning the centromere) carries
  **dcM = 0.000000**. Under decision C as scoped this is an interior zero-cM plateau and
  is deliberately untouched.

⚠⚠ **THE CONSEQUENCE, AND IT IS A LIVE RESIDUAL — NOW WITH A SIZE.** The superseded cut
ran to 12,968,320, which sits **above** the map's first marker at 10,326,676 — so 2.64 Mb
of what the cut removed is INSIDE the mapped span. Decision C re-admits those positions to
the chain, and because their `dcM` is a measured zero they will be emitted as a **verbatim
donor run**. So decision C does not simply restore the p arm: it restores part of it under
the declared rate and part of it under `P(switch) = 0`.

**MEASURED 2026-10-10 on our 340,849-record build**, resolving two figures in these docs
that did not reconcile:

| class | positions | `support > 0` | declared rate? |
|---|---|---|---|
| below the map start (< 10,326,676) — map-silent | **26,044** | 26,043 | **applies** |
| in [10,326,676, 12,968,320] — in-span zero-cM plateau | **8,327** | 8,327 | **does NOT** |
| sum | **34,371** | — | — |
| above the map end (> 46,680,243) — map-silent | **592** | — | **applies** |

The sum matches the recorded below-cut count exactly, and the plateau is genuinely flat:
**one** distinct cM value (0.584144) among all map markers in that interval. ⚠ The 26,018
that circulated in these docs is PanMixer's 340,824-record build, and 26,018 + 8,327 =
34,345 is why a second total was also on record. Ours is 26,044 + 8,327 = 34,371.

**So the residual is 8,327 chr21 positions — 2.7% of a 305,886-position chain — all of
them usable (`support > 0`) and all at `P(switch) = 0`.** Whether that is acceptable is
open question 2 and is Dylan's call; the point is that it is a known quantity, not an
unknown. The same split on the other 21 autosomes IS unmeasured and belongs to
`B-DAT-map-coverage`. In Arm 1 the donors are the protected population, so a verbatim run
there is a privacy cost, not merely a utility one.
Together that is **34,371 positions excluded by the cut** (26,044 NaN below the map plus 8,327 on the map's zero-cM centromeric plateau, all at cM 0.584144 exactly), so `Delta_x = 0`,
`P(switch) = 0`, and the sampler emits **one donor's real haplotype verbatim across
7.25 Mb** — departure (b) at roughly 200x the scale of the worst single block. This is a
privacy failure, not a utility one, and it is silent.

**Five independent reasons, each sufficient** (all measured or cited 2026-09-28). ⚠ Under
decision C these are no longer grounds for EXCLUSION — they are the recorded justification
for a DECLARED rate instead of an interpolated one, and the record of what the released
p-arm calls rest on. **Reason 1 inverts into a live hazard:** a constant cM/Mb applied to
fabricated N-gap distance MANUFACTURES genetic distance, so the declared rate must NOT be
accumulated across the fabricated N-gaps — see `declared_rate`.

| leg | evidence |
|---|---|
| the sequence is not real | GRCh38 chr21 5,010,000–10,814,560 is placeholder model sequence: 25 contigs separated by fabricated N-gaps, 22 of them exactly 50 kb. Base-pair distance there is fictitious |
| unmappable | 1000G's own strict accessibility mask rules **26,014 of 26,018** (99.98%) of those variants inaccessible to short reads |
| our own data is thin | cohort matrix is **36.76%** missing below the cut against **6.57%** above |
| no measured recombination | deCODE's pedigree map (4.5M observed crossovers) assigns **0.0 cM/Mb to every 1 Mb bin from 0–13 Mb** |
| none observed directly | a four-generation CEPH pedigree: *"Not a single allelic recombination was observed on the p-arm"* — 107 transmissions across ~38 Mbp of acrocentric short arms |

✅ **REVISITED AND SETTLED IN PART (decision C, 2026-10-08).** The option taken is none of
the ones below: the region is kept and given a DECLARED CONSTANT RATE, which is why the
rejections below still matter — they are the reasons an interpolated map was refused. What
remains OPEN is the rate's numeric VALUE and whether it is global or per-chromosome. The
original revisit list, recorded verbatim:
- **Maps that nominally cover the region exist.** pyrho hg38 (Spence & Song 2019, Zenodo
  11437540) covers chr21 from **5,088,754 bp** and gives every variant a distinct cM; the
  Eagle redistribution of our own HapMap map carries **66 markers** and 0.6166 cM inside the
  cut region. Both were rejected because they spend ~10 cM across intervals defined by two
  markers — 23–26% of chr21's genetic length attributed to a region with zero observed
  crossovers. That is interpolation across a void, not measurement.
- **The graph DOES carry the p arm.** The same HPRC graph's CHM13-backbone VCF has
  **185,414 records** in `chm13#chr21:1-10,000,000` against 22,746 for GRCh38, which cannot
  express anything below 5,010,001. Recovering it would mean a CHM13-backbone build with its
  own segmentation — it recovers sequence, not genetic distance — and missingness there is
  **43.9%** against 1.5% on the q arm, which collides hard with RENORMALISE.
- **T2T maps cannot simply be lifted back.** 774,011 bp of the liftable p arm lands on CHM13
  chr13/14/15/22, and the part reaching hs1 chr21 spreads from 2.48 to 44.74 Mb, so any cM
  assignment built that way is non-monotonic in GRCh38 coordinates.
- **What the cut costs is UNMEASURED.** The share of a target's `-log f` information
  discarded by excluding the p arm has not been computed. Compute it and report it beside
  the `target_fidelity` ceiling.

✅ **THE TELOMERIC END IS NOW CUT TOO (2026-09-29).** The **592 variants above the map
end** (map ends 46,680,243; positions run to 46,699,788) formed a **593-position
zero-distance run** — structurally the identical defect to the p arm, 58x shorter. They are
excluded, on the same principle: extrapolating there would invent genetic distance in a
subtelomeric region, which is exactly what we declined to do at the other end. Cost: 592 of
306,479 = **0.19%**.

⚠ **AND THE TELOMERIC END IS MAP-SILENT, SO DECISION C REACHES IT TOO.** Those 592
positions are above the map's last marker — no cM at all — so under decision C they are
kept and sampled under `declared_rate` rather than excluded, on the same footing as the
p arm. The refusal to EXTRAPOLATE the map stands; a declared rate is not an extrapolation
of the map.

⚠ **THE no-NaN INVARIANT IS RESTATED, NOT KEPT.** The old invariant — *every chain position
carries a genuine interpolated cM* — is **false by design** under decision C, and so is the
`assert no NaN` reading of the gate that rested on it. The replacement invariant: *every
chain position carries a FINITE cM from exactly ONE of the two DECLARED sources — map
interpolation or the declared constant rate — the SOURCE is recorded per position, and
cumulative cM is non-decreasing along the chromosome.* ⚠ NOT YET IMPLEMENTED:
src/cohort/genmap.py:chain_mask still RAISES on a non-finite cM inside the span
(src/cohort/genmap.py:59-88), which is the old invariant, and no per-position source array
exists. The gate must be replaced rather than widened — widening it to accept NaN would
discard the check entirely.

⚠ **ONE ZERO-DISTANCE CASE REMAINS, and it is correct to leave.** **618 interior map
intervals have dcM = 0**, covering 2.71 Mb, with ~1,650 kept variants inside them. These are
REAL map plateaus, and `P(switch) = 0` between two variants with no genetic separation is
the CORRECT model, not a defect.

✅ **AND DECISION C DOES NOT TOUCH THEM (2026-10-08).** The declared constant rate applies
ONLY to the MAP-SILENT regions — no cM at all. Interior plateaus sit INSIDE the mapped span,
where the map does speak and says zero, so the repo's position above STANDS and is not
stale. It is enforced as a Theorem-1 SUPPORT property by
tests/test_toy_enumeration.py:154-185 (`test_a_forced_stay_step_really_forbids_switching`),
whose docstring cites this very block: a tilt that pushed mass onto a path requiring a
switch across a zero-distance step would be emitting a path `R_D` assigns zero probability,
which is outside the support the theorem is stated over.

⚠ **WHAT IS DEFERRED, AND IT IS NOT NOTHING.** Whether interior plateaus are eventually
overridden is OPEN, and it is deferred until a PER-AUTOSOME gap measurement exists — the
chr21 figures above (618 intervals, 2.71 Mb) are the only ones measured, and the genome-wide
extent over the 22 autosomes is **UNMEASURED** (`map_coverage`). Measure it before arguing
the question either way, and do not carry a chr21 figure into a genome-wide claim. The
reason it needs settling at all: where `P(switch) = 0` over a long run, the sampler emits
one donor's real haplotype verbatim, and in **ARM 1** that donor is a cohort member whose
protection is in scope — the same departure (b) the p-arm cut was made to avoid. In **ARM 2**
the donors are the public cohort and this project makes no claim about protecting them.

✅ **WHAT CLOSING BOTH ENDS SIMPLIFIES.** Two planned artifacts fall away. `genetic_pos_valid`
is unnecessary — inside the span every position has a real cM by construction. And explicit
chain SEGMENTS are unnecessary: a long variant-free stretch inside the span produces a large
`Delta_x`, hence a switch probability approaching uniform, which is already the correct
behaviour. The chain does not need to be broken by hand anywhere inside `chain_span`.

⚠ **INTERPOLATE WITH `left=nan, right=nan`**, then assert no NaN survives the cut. That one
change converts a silent 34,371-position failure into a crash. ⚠ Under decision C the FIRST
half of that rule stays — never clamp, so the map-silent positions are identifiable rather
than silently given the map's first value — and the SECOND half changes: a NaN from the
interpolation is now the SIGNAL that a position takes the declared rate, so the assertion
becomes "every chain position ends with a finite cM and a recorded source", per the restated
invariant above. Clamping remains forbidden: it is what hid 34,371 positions in the first
place.


### support  [dataset] [implementation]
```
Is:                   The number of COHORT HAPLOTYPES with a CALLED allele at position
                      v, written support(v). Out of `K_states` = 88 (the 44 diploid
                      cohort members). It is a property of D alone — never of the
                      target — which is what makes every use of it below safe under
                      standing constraint 1.
                      It does THREE separate jobs, and conflating them is the drift
                      this entry exists to stop:
                      (1) CHAIN MEMBERSHIP. support(v) == 0 is the ONLY ground on which
                          a position leaves the chain after decision C. With no called
                          allele anywhere in D the sampler has nothing to emit and
                          RENORMALISE has an empty pool, so the position is not
                          sampleable rather than merely uninformative.
                      (2) THE WEIGHT, via `phi_v` = 1/support(v) under
                          `w_v = inv_support`, so a position only a few donors speak to
                          counts for more than a position everyone shares.
                      (3) THE DENOMINATOR of the per-position emission, since the
                          baseline R_D draws a donor from the called ones.
Computed by:          src/cohort/vcf.py at parse time, stored as data/cohort/chrN/
                      support.npy (int, one row per site on the site axis). ALSO read
                      independently from the VCF's INFO/AN field into an_info.npy.
                      ⚠⚠ **TWO INDEPENDENT ROUTES TO ONE QUANTITY, SO THEY CAN
                      DISAGREE** — one counts non-missing genotypes, one reads a field
                      the variant caller wrote. scripts/verify_cohorts.py checks all
                      three against each other (support.npy, the recomputed count, and
                      an_info.npy) on every built chromosome, and it is the
                      highest-power check in that script because a row misalignment
                      shows up here and nowhere else.
CHANGES MEANING WITH: the cohort D and the missingness policy (`missing_policy`). ⚠ It
                      changes under ARM 1 LEAVE-ONE-OUT: both of the target's
                      haplotypes leave the panel, so support(v) drops by up to 2 and
                      positions the target alone carried fall to 0 and LEAVE the chain.
                      So support(v) is per-target in arm 1 and cohort-wide in arm 2 —
                      never quote one for the other (standing constraint 8).
Valid vs:             the same cohort, the same arm and the same leave-one-out status.
Status:               PRIMARY — an input to the weight, the chain and the emission.
Aliases:              "AN", "allele number", "called count" — ⚠ `an_info` is the VCF
                      FIELD, a corroborating second route, not a synonym for the
                      quantity.
```

### chain_rule  [dataset] [method]
```
Is:                   WHICH policy decided each position's cM, recorded as one string
                      in every chromosome's manifest. Exactly two values exist:
                      `declared_rate_fill` — the post-decision-C rule: the chain runs
                          over the whole assembled autosome, positions inside the map's
                          marker span take an INTERPOLATED cM, positions outside it take
                          a cM accumulated at `declared_rate` over gap-free sequence,
                          and nothing is dropped for want of a cM.
                      `coordinate_cut` — the SUPERSEDED rule: positions outside the
                          marker span were EXCLUDED from the chain entirely
                          (`chain_span`).
                      ⚠ It exists because the two rules produce trees of the SAME SHAPE
                      with DIFFERENT CONTENT. Under the cut a map-silent position is
                      absent; under the fill it is present with a declared-rate cM. Both
                      yield a complete, self-consistent chromosome, so without this
                      field the only way to tell which one you are holding is to notice
                      that `cm_source.npy` is missing.
Computed by:          src/cli.py, written into data/cohort/chrN/manifest.json. It
                      REPLACED a required `chain_span` key — a manifest reader that
                      demands the span cannot express a declared-rate build.
CHANGES MEANING WITH: nothing; it is the policy label itself. But EVERY cM-dependent
                      quantity changes with it: `stay_switch`, `expected_switches`,
                      `T_blocks`, and every denominator that counts chain positions.
Valid vs:             the SAME rule. ⚠⚠ **MIXING RULES ACROSS CHROMOSOMES IS A
                      STANDING-CONSTRAINT-8 VIOLATION AT GENOME SCALE**, and it is the
                      live risk, not a hypothetical: `beta_v` sums to 1 over all 22
                      autosomes (standing constraint 4), so ONE chromosome built under
                      the other rule contaminates the genome-wide budget. chr21's tree
                      was pre-decision-C and carried no `chain_rule` at all.
                      scripts/verify_cohorts.py fails any chromosome not built under
                      `declared_rate_fill`.
Status:               PRIMARY — a genome-wide number is uninterpretable without it.
Aliases:              none. ⚠ Do not call it "the span"; `chain_span` is one of its two
                      values' input, not the field.
```

### cm_source  [dataset] [implementation]
```
Is:                   PER POSITION, which of the two cM sources `chain_rule`
                      = `declared_rate_fill` actually used: 0 = interpolated from the
                      genetic map's markers, 1 = accumulated at `declared_rate` over
                      gap-free sequence. The `chain_rule` label says which POLICY built
                      the chromosome; this says what the policy DID at each position,
                      which is the resolution any cM-dependent claim needs.
                      ⚠ Its existence is what makes a declared-rate build auditable. A
                      position's cM is a number either way, and a number carries no
                      record of how it was obtained.
Computed by:          src/cohort/genmap.py:fill_map_silent, which returns
                      (cm_filled, source) as a PAIR so the two cannot be written
                      separately and drift. Stored as data/cohort/chrN/cm_source.npy.
CHANGES MEANING WITH: the genetic map's marker span and the declared rate. ⚠ A
                      position's source is a property of the MAP and the ASSEMBLY, never
                      of the target, so it is fixed before any target is opened.
Valid vs:             the same map and the same rate.
Status:               PRIMARY for provenance, DIAGNOSTIC as a quantity — the SHARE of
                      chain positions on the declared rate is what makes a chromosome's
                      map coverage legible (`map_coverage`), and it varies by an order
                      of magnitude across autosomes.
                      ⚠ numpy.interp cannot leave an INTERIOR gap, so a declared-rate
                      position inside the mapped span means the fill was applied where
                      it should not have been. scripts/verify_cohorts.py checks exactly
                      that, and it is a check that CAN fail rather than a restatement of
                      the implementation.
Aliases:              none. ⚠ Not "the gap mask" — assembly gaps are SKIPPED when
                      accumulating declared distance (`declared_rate`) and are a
                      different thing.
```

### deconstruct  [implementation] [dataset]
```
Is:                   The FIRST and most expensive pipeline step: `vg deconstruct -P
                      grch38 -a` turns the HPRC PGGB graph into a VCF of variant
                      records on the GRCh38 backbone, one column per assembly. It is
                      what converts a pangenome graph into the site axis everything
                      else is indexed by, so `positions.npy`, `haplotypes.npy` and
                      every denominator in the project descend from it.
                      ⚠ Records may SHARE a start coordinate, because the graph's
                      bubbles NEST: a parent bubble and a variant inside it can begin
                      at the same reference position. The `lv` field is the nesting
                      level and 10-18% of records on every autosome carry `lv` > 0, so
                      a tied position is normal output (`nested_variant`). chr9 has one
                      such tie at 37,752,316 and it is correct.
Computed by:          src/cohort/vcf.py:deconstruct, via src/cli.py step
                      `deconstruct`, writing data/temp/chrN/deconstruct.vcf plus its
                      provenance sidecar.
                      ⚠ COSTS ONE TO TWO HOURS OF vg PER CHROMOSOME and about 10-15 GB
                      of RSS, which is what sizes every batch in
                      scripts/build_all_cohorts.py (22.3 MB of RSS per MB of gzipped
                      GFA, measured). It is RESUMABLE: a complete VCF on disk is
                      reused, and a complete VCF whose sidecar is missing gets one from
                      vcf.attach_deconstruct_prov without re-running vg.
CHANGES MEANING WITH: **the vg VERSION, which is load-bearing** — it decides ALT allele
                      ordering and therefore every allele index downstream, so the
                      binary's sha256 is verified against the config before it runs and
                      a mismatch REFUSES rather than warns. Also with the reference
                      prefix: `-P grch38` is what makes GRCh38 the backbone, and a
                      different prefix would produce a different coordinate system
                      under the same file name.
Valid vs:             the same graph, the same vg binary and the same reference prefix.
                      ⚠ Never compare records across two deconstruct runs without
                      checking both sidecars' `gfa` and `vg` digests.
Status:               PRIMARY — the root of the whole data lineage.
Aliases:              "the VCF" is ambiguous; this is the 45-sample
                      `deconstruct.vcf` BEFORE chm13 is dropped. The 44-sample product
                      is `cohort44.vcf.gz`, from a separate step.
```

### in_chain  [dataset] [implementation]
```
Is:                   The CHAIN MASK: a bool over the full site axis saying which
                      positions the sampler actually runs over. `T` is its sum, and
                      **every denominator in this project is a count of chain
                      positions**, so a change here moves every rate, share and
                      fidelity figure at once.
                      ⚠ EXACTLY TWO EXCLUSIONS, and both are target-INDEPENDENT:
                      (1) support(v) == 0 — no donor has a called allele, so the
                          sampler has nothing to emit and RENORMALISE has an empty pool.
                      (2) an `untraversable_block` (decision K) — a maximal
                          zero-genetic-distance run that no single donor spans, which
                          the copying model has no path through.
                      ⚠ A THIRD exclusion is a defect. The coordinate cut
                      (`chain_span`) was the third and is SUPERSEDED by decision C.
Computed by:          src/cohort/genmap.py:chain_mask then :drop_untraversable, both
                      inside genmap.build, stored as data/cohort/chrN/in_chain.npy as a
                      MASK over the full axis — never as a row filter, so the site axis
                      stays the join key.
                      ⚠ scripts/verify_cohorts.py RECOMPUTES both exclusions from the
                      data rather than reading the manifest's own account of them, which
                      is what lets it catch a chain that was never rebuilt.
CHANGES MEANING WITH: the cohort D, the genetic map, AND THE ARM. ⚠ Under ARM 1
                      LEAVE-ONE-OUT both of the target's haplotypes leave the panel, so
                      support(v) drops by up to 2, positions the target alone carried
                      fall to 0 and LEAVE the chain, and an intersection that was
                      non-empty can empty. Arm 1's chain is therefore PER TARGET. Never
                      reuse arm 2's mask for an arm 1 release, and never compare a
                      denominator across arms (standing constraint 8).
Valid vs:             the same cohort, map, arm and leave-one-out status.
Status:               PRIMARY. MEASURED 2026-10-10, genome-wide: 22,307,009 records ->
                      **22,134,562** chain positions, excluding 305 for support == 0
                      (0.0014%) and 172,142 in 20 untraversable blocks (0.772%).
Aliases:              "the chain", "T" (its sum). ⚠ `T_blocks` is the MODEL-B unit count
                      and a different quantity; do not use "T" for both in one table.
```

### untraversable_block  [dataset] [method] [implementation]
```
Is:                   A maximal run of chain positions with ZERO genetic distance
                      between neighbours that NO single cohort donor spans. Decision K
                      (Dylan, 2026-10-10) EXCLUDES the whole block from the chain.
                      Why no path gets through: the transition is
                      A_t = (1 - s_t) I + (s_t / K) J and s_t is exactly 0 when
                      Delta_x is 0, so A_t is exactly the identity and no switch is
                      possible anywhere inside the run. The donors able to traverse it
                      are therefore the INTERSECTION of `allowed` over EVERY position in
                      the run. Empty intersection -> the forward pass dead-ends with
                      zero total mass.
                      ⚠⚠ IT IS NOT A PAIRWISE PROPERTY. Every adjacent pair in a run can
                      share a donor while no donor spans the whole run, so a pairwise
                      check passes and the release still dies. Three donors are the
                      minimum needed to construct that case.
                      ⚠ The model is not wrong here: P(switch) = 0 on a zero-cM plateau
                      is CORRECT (see `chain_span`). It is the DATA that makes the state
                      space disconnected.
Computed by:          src/cohort/genmap.py:untraversable_blocks, applied by
                      genmap.build in decision-C mode only, recorded per chromosome in
                      the manifest as `untraversable_blocks` with each block's bp span.
                      ⚠ scripts/verify_cohorts.py recomputes them INDEPENDENTLY of the
                      code that wrote the mask, which is what gives that check power in
                      both directions -- an extra exclusion AND a block that should have
                      been dropped and was not.
                      ⚠ NO ITERATION: cm is non-decreasing and each block is maximal, so
                      a dropped block's neighbours are strictly separated and cannot form
                      a new zero-distance step. genmap.build asserts it anyway.
CHANGES MEANING WITH: the genetic map and the cohort D. ⚠ NEVER with the target -- it
                      reads only cm and `allowed`, which is what keeps it inside standing
                      constraint 1; a target-dependent exclusion would void Theorem 1.
                      Asserted from the function signatures in
                      tests/test_untraversable.py.
                      ⚠ IT WOULD CHANGE UNDER ARM 1 LEAVE-ONE-OUT, because removing a
                      target's two haplotypes can empty an intersection that was
                      non-empty. Arm 1's chain is therefore per-target, like its
                      support(v) and its beta_v -- do not reuse arm 2's mask there.
Valid vs:             the same cohort, the same map, and the same arm.
Status:               PRIMARY — it is one of the chain's TWO exclusions and every
                      denominator in the project counts chain positions.
                      MEASURED 2026-10-10: of 52,870 zero-distance runs genome-wide only
                      **20** are untraversable, **172,142 positions, 0.772% of the
                      chain**; 14 of the 20 overlap a centromere. chr9 loses 7.15% of its
                      chain, chr20 3.58%, chr21 2.44%; chr10, chr14 and chr22 lose none.
                      ⚠ THE WHOLE BLOCK GOES, not the minimum needed. The best-covering
                      donor spans 72-100% of each block, so keeping its sub-run would
                      retain about 96% of these positions -- rejected, because a single
                      donor is FORCED there and those positions would enter every
                      per-position denominator while carrying no variability.
Aliases:              "the centromeric plateaus" — ⚠ APPROXIMATE and not a synonym: 14
                      of 20 overlap a centromere, 6 do not (chr1 twice, chr6, chr13,
                      chr15). "dead-end block" is acceptable in prose.
```

### declared_rate  [dataset] [method]
```
Is:                   The DECLARED CONSTANT RECOMBINATION RATE, in cM/Mb, used to give
                      a genetic distance to the MAP-SILENT regions of an autosome —
                      the regions where the genetic map gives NO cM AT ALL. Introduced
                      by decision C (Dylan, 2026-10-08), which replaced the coordinate
                      cut: those regions are now SAMPLED rather than excluded, so the
                      chain needs a `Delta_x` there, and a declared constant is the one
                      choice that does not pretend to be a measurement.
                      ✅ **VALUE SETTLED 2026-10-10 (Dylan): ONE GLOBAL
                      1.2707 cM/Mb** for all 22 autosomes. It is the autosomal mean
                      MAPPED rate — 3,543.8 cM over 2,788.8 Mb of mapped span,
                      measured by `scripts/map_coverage.py`.
                      ⚠ **IT IS A DECLARED CONVENTION, NOT AN ESTIMATE OF THE REGION
                      IT IS APPLIED TO.** It is used where the map is silent precisely
                      because no map speaks there, and the biology points the other
                      way: deCODE's pedigree map assigns 0.0 cM/Mb to every 1 Mb bin
                      from 0-13 Mb of chr21. Report it as a declaration.
                      ⚠ **WHY GLOBAL AND NOT PER-CHROMOSOME, which reversed the earlier
                      recommendation.** 96% of the map-silent sequence is the five
                      acrocentric p arms, and those five have the genome's HIGHEST
                      mapped rates (chr22 2.046, chr15 1.727, chr21 1.711 cM/Mb, all
                      above the mean). An acrocentric's rate is measured entirely on
                      its q arm, because the p arm IS the region being filled — so a
                      per-chromosome rate would extrapolate an unusually hot q-arm rate
                      onto heterochromatin and rDNA, and would do it on exactly the
                      chromosomes where the rule matters most. Spread 1.84x.
Computed by:          DECLARED, in `configs/inputs.chr21.yaml` under `declared_rate`,
                      and pre-registered with tau, the utility function, the endpoints,
                      the attacker set and the attack database under constraint 6.
                      ⚠ A test asserts the declared number still EQUALS the measured
                      autosomal mean, so the config and the maps cannot drift apart —
                      a declared constant that no longer matches its own stated
                      derivation looks measured and is not.
                      ⚠ STILL NOT IMPLEMENTED: nothing in `src/cohort/genmap.py` can
                      assign a cM to a map-silent position — it interpolates with
                      left=nan/right=nan and raises if a NaN survives, which is the
                      invariant decision C makes false by design. That fill is
                      preprocessing work, not part of the declaration.
CHANGES MEANING WITH: the reach it is given and the distance it is applied to.
                      ⚠ IT APPLIES TO THE MAP-SILENT REGIONS ONLY. The INTERIOR
                      zero-cM map plateaus inside the mapped span are NOT overridden —
                      see `chain_span`. Applying it there is a different decision,
                      still open.
                      ⚠ AND IT MUST NOT BE ACCUMULATED ACROSS FABRICATED SEQUENCE. The
                      largest map-silent region on chr21 is placeholder GRCh38 model
                      sequence (5,010,000-10,814,560: 25 contigs separated by fabricated
                      N-gaps, 22 of them exactly 50 kb), where base-pair distance is
                      fictitious. A constant cM/Mb times a fabricated bp span
                      MANUFACTURES genetic distance. The declared policy must say what
                      happens at an N-gap before any rate is applied.
                      ⚠⚠ **STILL UNDECLARED AND BLOCKED as of 2026-10-10**, and the
                      settled 1.2707 does NOT settle it. It needs the GRCh38 assembly's
                      gap coordinates. ✅ ACQUIRED 2026-10-10 as
                      `data/raw/ucsc_gap/gap.txt.gz`, so this is no longer BLOCKED —
                      only undeclared. ⚠ CORRECTION: this read "no FASTA, AGP or
                      assembly report is on disk", which was false; PanMixer's tree
                      holds a full GRCh38 FASTA. The true statement was "none staged in
                      OUR tree", and reading inside that tree is forbidden.
                      ✅ MEASURED: **79.8% of the map-silent sequence is gap** (68.8 of
                      86.2 Mb), so accumulating the rate across it asserts 109.5 cM of
                      which **87.4 cM is on placeholder bp**. **No variant of ours is
                      inside any gap** (chr21: 0 of 340,849), so skipping gap spans
                      changes NO position's chain membership and drops the map-silent
                      distance to 22.1 cM.
                      ⚠ Until it is decided the implementation can only ACCUMULATE
                      across the gaps, because it cannot locate them — a consequence of
                      the missing input, not a choice, and it must be stated in any
                      writeup. On chr21 that is 1.1 Mb of fabricated gap inside the
                      10.33 Mb the rate applies to; the extent on the other 21
                      autosomes is unmeasured.
                      ⚠ AND WITH THE MONOTONICITY IT IMPLIES: the declared-rate cM has
                      to be stitched to the mapped cM so cumulative cM stays
                      non-decreasing across the boundary in GRCh38 coordinates.
Valid vs:             itself, as a declared constant, and against the measured maps
                      that were REJECTED for the same region (pyrho hg38 from
                      5,088,754 bp; the Eagle redistribution's 66 markers and 0.6166 cM
                      inside the old cut) — rejected because they spend ~10 cM across
                      intervals defined by two markers, 23-26% of chr21's genetic length
                      in a region with ZERO observed crossovers. A declared constant is
                      honest about being a declaration; those were interpolation across
                      a void wearing the clothes of measurement.
Status:               PRIMARY — it supersedes deferred decision R2 and it is on the
                      critical path: `chain_span`, T, `beta_v`, `expected_switches` and
                      every denominator move with it.
Aliases:              "the constant rate", "the declared rate". ⚠ NEVER "the
                      genome-wide average rate" (~1.2 cM/Mb is a POPULATION-GENETIC
                      fact about mapped sequence, not this project's declaration) and
                      never "the interpolated rate" — the whole point is that it is
                      declared, not interpolated.
```

### map_coverage  [dataset] [method]
```
Is:                   How much of an autosome the genetic map actually speaks to, and
                      therefore how much of it decision C hands to `declared_rate`.
                      Per autosome, three quantities that must be reported together:
                      (a) the MAPPED span's bp and the count of VCF records inside it;
                      (b) the MAP-SILENT bp and record count — below the first marker,
                      above the last, and any gap in the middle with no marker at all,
                      which is what decision C's declared rate covers; and (c) the
                      INTERIOR ZERO-cM bp and record count — inside the mapped span,
                      where the map speaks and says zero — which decision C does NOT
                      cover and which is reported separately so the two are never added.
Computed by:          ✅ `scripts/map_coverage.py` (2026-10-10), output at
                      `logs/B_map_coverage.md`. MEASURED in base pairs for all 22
                      autosomes: map-silent **3.00%** (86.2 Mb, 96% of it the five
                      acrocentric p arms), interior zero-cM **3.27%** (93.9 Mb,
                      chr9 alone 25.61 Mb), autosomal mean mapped rate **1.2707
                      cM/Mb** with a 1.84x per-chromosome spread.
                      ⚠ POSITION counts are chr21-only, because they need each
                      chromosome's own deconstructed VCF: 26,044 below the map
                      start, 8,327 on the in-span plateau, 592 above the map end.
                      Base pairs are NOT positions and the rate acts per position,
                      so the bp table must not be read as a substitute.
                      ⚠ The N-GAP class is still BLOCKED: it needs the GRCh38
                      assembly's gap coordinates and no FASTA, AGP or assembly
                      report is on disk.
                      ✅ **chr21 IS measured, 2026-10-10, on our 340,849-record build:**
                      **26,044** records below the map's start (10,326,676), of which
                      26,043 have `support > 0`; **8,327** in [10,326,676, 12,968,320],
                      the in-span zero-cM plateau, all with `support > 0`; the two sum to
                      the recorded **34,371** dropped by the superseded cut; and **592**
                      above the map's end (46,680,243). ⚠ The **26,018** previously
                      quoted here is PanMixer's 340,824-record build — pairing it with
                      8,327 gives 34,345, which is why two different totals were on
                      record. Ours is 26,044. Plus 618 interior zero-cM
                      intervals covering 2.71 Mb with ~1,650 kept variants.
CHANGES MEANING WITH: the MAP (a different map has different markers, so this is a
                      property of the map-and-assembly pair, never of the target), the
                      autosome, and above all WHICH OF THE THREE QUANTITIES is being
                      quoted. ⚠ Conflating (b) with (c) is the specific error decision
                      C's scope narrowing exists to prevent: they get different
                      treatment, so a single "zero-distance bp" number is not usable.
Valid vs:             the same map on the same build, per autosome. A genome-wide total
                      is quotable only once all 22 are measured — the 22 are NOT
                      interchangeable (chr21 is acrocentric; most are not).
Status:               PRIMARY — the deferred interior-plateau question is explicitly
                      deferred UNTIL this exists, and `declared_rate`'s reach cannot be
                      stated without it.
Aliases:              "map silence", "the no-cM regions". ⚠ UNMEASURED genome-wide. Any
                      genome-wide Mb or percentage for this quantity is unverified: do
                      not quote one, and do not scale chr21's figures to 22 autosomes.
```

### expected_switches  [implementation] [gate]
```
Is:                   E[number of donor changes] along the chain under the PRIOR,
                      sum_t s_t (K-1)/K. A property of the cohort, the genetic map
                      and the chain architecture -- never of the target.
Computed by:          src/mechanism/hmm.py:CohortHMM.expected_switches (built
                      2026-09-30)
CHANGES MEANING WITH: the chain ARCHITECTURE above all. Measured on chr21 at K = 88,
                      ARM 2, under the SUPERSEDED coordinate cut (T = 305,886):
                      **7,021.0** under model C and **3,394.1** under model B, against
                      7,155 and 3,532 predicted independently on 2026-09-28.
                      ⚠ BOTH MOVE UNDER DECISION C, in both of their inputs: the chain
                      gains the map-silent positions, and those positions carry a
                      `declared_rate` cM instead of a zero one, so they contribute
                      switch mass where they previously contributed none and were then
                      excluded. The post-decision-C value is UNMEASURED and its
                      independent prediction has to be re-derived with the rate, which
                      is unchosen. ⚠ AND WITH K: at ARM 1's K = 86 the (K-1)/K factor
                      and the per-state mass both change, so an Arm-1 figure is not an
                      Arm-2 figure. Also with
                      the genetic map, which is what makes it a canary -- PanMixer's
                      shipped sampler gives about 0 because its distance is a VCF row
                      index.
Valid vs:             its own independent prediction, and the analytic segment length
                      n/(Ne*r) = 0.0068 cM. It is the REPLACEMENT for mis-specified
                      exit gate 3, which demanded switch mass of order 0.1-1 between
                      adjacent positions and would have rejected a correct model.
Status:               PRIMARY — mark * in tables. ⚠ It is the only quantity in the
                      project with an independently predicted value, and that is what
                      caught the model-B transition bug (`docs/bugs.md` 2026-10-07):
                      nothing else disagreed, no test failed, and the wrong number
                      looked entirely plausible. Check it whenever the chain changes.
Aliases:              "E[switches]". ⚠ NOT `donor_switches`, which is the REALISED
                      count in one release; this is the expectation before any draw.
```

### donor_switches  [implementation] [evaluation]
```
Is:                   The number of positions at which a RELEASED path changes donor
                      haplotype, counted on the sampled state path. The measured
                      counterpart of E[donor switches] under the prior, and the
                      quantity the model-C architecture decision turned on.
Computed by:          src/mechanism/ffbs.py:Release.record (built 2026-09-30);
                      the prior expectation is CohortHMM.expected_switches.
CHANGES MEANING WITH: the PLOIDY — a diploid release reports the sum over both
                      haplotypes, so halve it before comparing against a per-haplotype
                      figure. And with the chain architecture: measured 2026-09-28,
                      104,379 under model A, 3,532 under B, 7,155 under C, against an
                      analytic segment length of n/(Ne*r) = 0.0068 cM (~4.5 kb).
                      ⚠ Also with the genetic map, which is what makes it a useful
                      canary: PanMixer's shipped sampler gives ~0 because its distance
                      is a VCF row index.
Valid vs:             expected_switches under the same cohort, the same arm and the
                      same chain_span, and against the analytic segment length.
                      ⚠ **EVERY FIGURE BELOW IS ARM 2, K = 88, UNDER THE SUPERSEDED
                      COORDINATE CUT (T = 305,886)** — stamped 2026-10-10 to match
                      `expected_switches`. All of them MOVE under decision C (more
                      positions enter the chain, and map-silent positions gain a
                      declared-rate cM where they previously had none) and under
                      decision E (Arm 1's K = 86 changes the (K-1)/K factor). The
                      post-decision-C values are UNMEASURED.
                      Measured 2026-09-30 on our
                      chr21 build at K = 88: prior expectation **7,021.0** per
                      haplotype, against the **7,155** the model-C decision predicted
                      independently; realised **14,169** for a diploid release
                      (7,084 per haplotype). ⚠ An earlier version said "IDENTICAL at
                      every tau from 0 to 0.9"; the per-draw values DO move —
                      14,147/14,211/14,147 at tau = 0.9 against 14,148/14,212/14,148 below
                      it — and only the rounded MEAN (14,169 either way) was identical.
Status:               SECONDARY — structure, not a gate. But a value near 0 or near T
                      means the transition model is broken, which is exactly how
                      PanMixer's defect would have shown up had anyone looked.
Aliases:              "switches", "recombinations". ⚠ Never "crossovers" — those are a
                      biological event; this is a copying-model bookkeeping count.
```

### K_states  [implementation] [dataset]
```
Is:                   K, the number of cohort haplotype states. At most ~2n for a
                      diploid cohort of n individuals, before duplicate haplotypes
                      are collapsed.
Computed by:          src/mechanism/hmm.py:CohortHMM.K = haplotypes.shape[0]
                      (built 2026-09-30).
CHANGES MEANING WITH: ⚠ **THE ARM, FIRST OF ALL (decision E, 2026-10-08).**
                      **ARM 2: K = 88** — the target is external, so the full cohort is
                      the donor panel and nothing is held out. **ARM 1: K = 86** — the
                      target is a cohort member and per-target leave-one-out removes
                      BOTH of its haplotypes. A K quoted without its arm is ambiguous,
                      and the two are not interchangeable: alpha storage, runtime, the
                      uniform initial distribution and every per-state switch
                      probability s_t/K are all functions of K. The K = 86 figures are
                      ~2.3% smaller than the K = 88 ones.
                      ⚠ Arm 1's K = 86 is applied by MASKING two haplotype rows in
                      memory at release time, not by building a second cohort tree, so
                      the K = 88 artifacts on disk remain the source and the release
                      record must carry the masked haplotype ids.
                      Also with whether duplicates have been collapsed. A K quoted before
                      collapse and a K quoted after are different numbers.
                      ⚠ **COLLAPSE IS NOT IMPLEMENTED (gap found by audit, 2026-09-30).**
                      No collapse, no policy, no declaration, and no assertion that
                      duplicates are absent. **EVERY K this project quotes is therefore
                      PRE-COLLAPSE**, including the runtime and alpha-storage figures.
                      Two donors identical across all chain positions would be two
                      states the mechanism treats as distinct, putting twice the prior
                      mass on one haplotype.
                      ✅ **BUT IT DOES NOT OCCUR ON THIS COHORT — MEASURED 2026-09-30:**
                      all 88 haplotypes are DISTINCT over the 305,886 chain positions,
                      and the closest pair differs at **32,741** of them (10.7%). So
                      every K = 88 quoted is correct; what was missing is that nothing
                      said so and nothing checked. The gap is a missing ASSERTION, not a
                      wrong number. ⚠ CORRECTED 2026-10-10: this field previously read
                      "Whether it OCCURS on this cohort is unmeasured. Measure it before
                      quoting K again" — stale, since the measurement was taken the same
                      day the gap was filed and is recorded in `docs/plan.md`.
                      ⚠ Removing two haplotypes cannot CREATE a duplicate, so the
                      measurement carries to Arm 1's K = 86 — but a different cohort or a
                      shorter chain could produce one, so the assertion is still worth
                      adding.
Valid vs:             the same cohort at the same collapse policy AND the same arm.
Status:               SECONDARY
Aliases:              "K", "donor states"
```

### diploid_utility  [utility] [theory]
```
Is:                   For a diploid target, the haplotype pair is ONE private input:
                      u_diploid(p, y) = (u(p1, y1) + u(p2, y2)) / 2. The average keeps
                      the joint utility in [0,1], so the same tau covers the complete
                      diploid release.
Computed by:          src/mechanism/utility.py:diploid_utility (built 2026-09-30).
                      ✅ THE MEAN FACTORISES EXACTLY, which is what makes it samplable:
                      exp(eta (u1+u2)/2) = exp((eta/2) u1) exp((eta/2) u2), so the two
                      haplotype chains are INDEPENDENT at HALF the exponent and log Z_p is
                      the sum of the two, landing in [0, eta/2] + [0, eta/2] = [0, eta_tau].
                      One diploid release, one tau, one bracket. ffbs.release handles
                      ploidy 1 and 2 through the same entry point.
CHANGES MEANING WITH: the aggregator. A SUM instead of a MEAN leaves [0,1] and
                      doubles the effective tilt — the guarantee no longer holds.
                      ⚠ Running each haplotype at the FULL eta_tau is that SUM in
                      disguise: measured on the toy at tau = 0.5, the two log Z_p values
                      then sum past eta_tau.
Valid vs:             a haploid release only after the aggregator is stated.
Status:               PRIMARY
Aliases:              none
```

---

## PanMixer terms (upstream)

Defined here because we consume them. A PanMixer name must never be used for one
of our quantities, and vice versa. Evidence is `external/PanMixer/` at commit
`c182c38`; paper references are to `archive_docs/Blindenbach2026_PanMixer.pdf`.

### nested_variant  [dataset] [implementation]
```
Is:                   A variant that sits INSIDE another variant's allele — a bubble
                      within a bubble in the graph. `vg deconstruct -a` emits every
                      level of the snarl tree, tagging each record with LV (level,
                      0 = top) and PS (parent bubble id).
Computed by:          upstream, by vg deconstruct; visible in the VCF INFO fields.
                      Measured on our chr21: 90.8% LV=0, 8.2% LV=1, 0.9% LV=2,
                      0.1% LV=3, 10 records at LV=4.
CHANGES MEANING WITH: whether the LV/PS tags survive. `VCFtoNP` DROPS them, so the
                      numpy matrix has no record of nesting at all — every row looks
                      independent when 9.2% of them are not.
Valid vs:             nothing; it is a structural property, not a measurement.
Status:               SECONDARY
Aliases:              "child variant". ⚠ A nested variant and its parent are BOTH
                      rows in the matrix, so any per-row sum counts the region twice.
```

### hypervariable_site  [dataset] [privacy]
```
Is:                   A site so multi-allelic that few or no haplotypes share an
                      allele. The extreme on our chr21 is grch38#chr21:14569980,
                      where all 88 haplotypes carry 88 distinct alleles. Typically a
                      tandem repeat or segmental duplication where each assembly
                      lands on a different copy number.
Computed by:          `num_alleles.npy`; 7 sites on chr21 have >=88 alleles, 155
                      have >=44, 618 have >=20.
CHANGES MEANING WITH: whether unused ALT alleles were trimmed. `view -s ^chm13`
                      removes the sample but not its alleles, so the declared allele
                      count can exceed the number actually carried (90 declared,
                      88 carried at the site above).
Valid vs:             the same site set on the same cohort.
Status:               DIAGNOSTIC — but a PRIMARY concern for the mechanism: these are
                      maximally identifying AND scored at zero by `eps_pmi`.
Aliases:              none. Do not call these "SVs" — size is not what defines them.
```

### missing_policy  [implementation] [method]
```
Is:                   What the mechanism does when a donor haplotype carries -1 (no
                      called allele) at a chain position. ✅ **DECIDED 2026-09-28:
                      RENORMALISE** - drop that donor from the state distribution at
                      that position. The rejected options are kept below because a
                      run under one is not comparable to a run under another.
                      Alternatives considered:
                      WILDCARD - treat missing as matching whatever the target has,
                      which makes poorly-assembled donors universally attractive;
                      MISMATCH - treat it as a difference, which penalises assembly
                      gaps as though they were genuine variation;
                      RESTRICT - exclude ill-behaved positions from the chain
                      altogether, which is what PanMixer's anchor rule does in
                      effect. Worth keeping as a baseline precisely because it is
                      the published behaviour, so a v1 run can be compared against
                      it (Dylan, 2026-09-25).
Computed by:          PLANNED
CHANGES MEANING WITH: which positions are in the chain. It barely matters for the
                      anchors PanMixer's HMM steps over (0.26% mean missingness,
                      538 over 25%) and matters a lot for v1, which steps over every
                      variant (4.26%, 30,800 over 25%). Worst are the non-anchor
                      records inside multi-variant blocks, the positions v1 adds
                      inside blocks PanMixer already chains (11.31%, 7,555 over 25%).
                      Scopes measured in `docs/memory.md` 2026-09-26; never quote
                      one of these means without naming its set. Nested variants are the driver: a donor whose path
                      does not traverse the parent bubble has no allele at the child,
                      and the matrix cannot distinguish that from missing data.
                      ⚠ RENORMALISE acts on the STATE SPACE, not the transition. Our
                      baseline is a pure prior over donor paths driven by cM distance
                      and knows nothing about alleles; the target enters only through
                      the tilt. A donor missing at v has no allele, so phi_v is
                      undefined for it and it is excluded from the state set at v.
                      The path therefore CANNOT sit on that donor at v.
                      ⚠ CONSEQUENCE TO MEASURE: leaving a donor and returning costs
                      TWO switch events (~0.0026^2 at typical 12 bp spacing), so ffbs
                      will prefer to switch away and STAY away. Missingness acts as a
                      switch TRIGGER, not a one-position blip — and it is spatially
                      clustered (28,052 runs on strand 0, mean length 20, max 23,478),
                      so expect systematic donor switching in badly-assembled regions.
                      Measure this once the sampler runs; it is not a reason to change
                      the policy now.
                      ⚠ GUARD REQUIRED: measured on OUR chr21 build, 1 site of 340,849
                      has all 88 donors missing and 41 sites leave <= 1 donor. The
                      renormalised denominator can be zero. Declare the fallback
                      (uniform over all donors), count how often it fires, and report it.
Valid vs:             another run under the SAME policy. Never compare across policies.
Status:               PRIMARY — DECIDED. An undeclared default here would be a silent
                      modelling choice; RENORMALISE is the declared one.
Aliases:              none
```

⚠ **`-1` conflates FOUR conditions on chr21, and they do not want the same treatment.**
Measured (see `docs/memory.md`, 2026-09-25 and 2026-09-28). A fifth cause — haploid
genotypes — was expected and **does not occur on this data**: all 237,594 non-piped GT
fields have exactly one distinct value, `.`, so there are ZERO true haploid genotypes and
the converter's `len(genotype) == 1` branch only ever fires on an explicit missing call.
The row is retained struck-through because it is a real hazard for any OTHER input VCF.

| cause | how it arises | how common |
|---|---|---|
| **not applicable** | the child bubble sits inside a parent allele this haplotype does not carry — the DNA does not exist on that chromosome | 66.1% of nested missingness |
| **inherited** | the parent call was itself missing | 12.1% of nested missingness |
| **assembly gap** | the donor's assembly does not cover the region | 349 runs of >=100 consecutive positions, covering 448,283 entries |
| **conflict** | the sample had multiple graph paths that disagreed, so vg wrote `.` | 1,855 records, mean 29.4 of 88 haplotypes |
| ~~**haploid genotype**~~ | ~~the VCF gave one allele, not two; the converter writes `-1` into strand 1~~ | **DOES NOT OCCUR on chr21** — all 237,594 non-piped GT fields are `.` (measured 2026-09-28) |

"Not applicable" argues for RENORMALISE — that donor genuinely has no allele there.
An assembly gap argues against it — the donor almost certainly HAS an allele and we
merely do not know it, so dropping them penalises poorly-assembled donors for reasons
unrelated to genetics. **A single policy is being asked to cover all four.** RENORMALISE was chosen knowing this:
it is CORRECT for "not applicable", which is the largest class (66.1% of nested missingness),
and it does penalise poorly-assembled donors for a defect that is not theirs. That cost was
accepted deliberately (Dylan, 2026-09-28) on the grounds that released pangenomes should be
well assembled; if they are not, it is something to note in the writeup rather than model
around at this stage.

⚠ **The causes ARE distinguishable, but the information is destroyed before the
mechanism sees it.** `LV` and `PS` separate "not applicable" and "inherited";
run-length separates assembly gaps; the `CONFLICT` INFO tag names conflicts;
haploid-vs-diploid is visible in the GT itself. But
`external/PanMixer/tools/common/VCFtoNP.py` keeps ONLY the position and the GT
subfield (`fields[1]` and `fields[9+i].split(':')[0]`) — every INFO field is dropped.
So by the time anything reaches the matrix, all five look identical.
**Preserving the cause is a cheap preprocessing addition — an auxiliary reason array
alongside the allele matrix — and is far cheaper now than reconstructing it later.**


✅ **HOW RENORMALISE IS ACTUALLY IMPLEMENTED (2026-09-30), and why the obvious reading is
wrong.** "Drop that donor from the state distribution at that position" invites masking
row `j` of the transition to the allowed set and renormalising it. **That is not well
defined.** Across a step with `Delta_x = 0` the transition is the identity, so a donor
available at `t-1` and unavailable at `t` has a masked row summing to **exactly zero** --
there is no chain to renormalise and the division is by zero. Measured on our chr21 build:
**1,058 zero-distance steps inside `chain_span`**, against 4.26% mean missingness, so this
is reachable on the real data rather than a contrived case.

So the mask enters as a **0/1 POTENTIAL** and `R_D` is the Li-Stephens chain CONDITIONED on
never occupying an unavailable donor. The forward pass normalises it away. ⚠ The
consequence is `baseline_log_norm`: the conditioning constant is part of `sum_t log c_t`,
so `log_z_p` is that sum MINUS the baseline's own, and a dead-end chain now raises rather
than dividing by zero.

### target_encoding  [dataset] [method]
```
Is:                   The form a target genome arrives in, and therefore what
                      `Map(g, G)` has to do to turn it into a path over G's sites.
Computed by:          n/a -- a property of the input we must build around.
CHANGES MEANING WITH: the encoding. The three that matter:
                      **VCF** (what we have): SPARSE. Lists the target's variants
                      and is SILENT everywhere else. Turning it into a dense path
                      requires INFERRING reference wherever it is silent -- and
                      "silent" conflates "this person is reference here" with "this
                      site was never in my callset". That conflation is the whole
                      of `absent_policy`, and getting it wrong at multi-allelic
                      records is what produced 6,309 fabricated reference calls per
                      target on 2026-09-29.
                      **gVCF**: DENSE. Records reference-confident BLOCKS as well as
                      variants, so "this person is reference from X to Y, with this
                      confidence" is stated rather than inferred. ⚠ **A gVCF would
                      dissolve `absent_policy` entirely** -- the ambiguity is an
                      artifact of the sparse encoding, not of the biology.
                      **Reads / assembly** (BAM, CRAM, FASTA): needs calling or
                      alignment first. Aligning straight to the graph with
                      `vg giraffe` would give a path directly and skip the join
                      altogether -- but HPRC marks short-read mapping "untested" for
                      the PGGB graph, so that route is not free either.
Valid vs:             the same encoding. A fidelity number from a VCF target and one
                      from a gVCF target are not comparable, because the second knows
                      something the first had to guess.
Status:               PRIMARY -- it determines what the private input even is.
Aliases:              none
```

⚠ **A SEPARATE REPRESENTATION MISMATCH, independent of sparseness.** Even a VCF that covers
every site can disagree with the graph about HOW a site is written. `bcftools norm -m-any` --
an extremely common normalisation, and what the 1000G panel ships -- splits a multi-allelic
record into one biallelic record per ALT. The graph keeps them together. So `C -> T,G` in the
graph is two records in the panel and a strict `(POS, REF, full-ALT)` join matches neither.
This is **not** a panel quirk: any target VCF normalised that way has it. `split_alt_join`
handles it, restricted to records the target's VCF fully covers, because only then is
"carries none of our alts" a conclusive statement rather than a gap.

### reason_array  [implementation] [dataset]
```
Is:                   An auxiliary array row- and slot-aligned to the allele matrix,
                      holding a small integer code per entry saying WHY that entry is
                      `-1`. Shape `(88, n_sites)` int8, code 0 = called.
Computed by:          src/cohort/arrays.py (built 2026-09-29) — our VCF-to-matrix step,
                      the only moment the INFO column is in hand.
CHANGES MEANING WITH: the code set and the precedence rule (below). Two reason arrays
                      built under different precedence are NOT comparable.
Valid vs:             the same code set, same precedence, same site axis.
Status:               PRIMARY — `missing_policy` is only expressible if the cause
                      survives conversion, and it does not survive by default.
Aliases:              "the REASON array", "cause codes"
```

✅ **DESIGN DECIDED 2026-09-29 — three per-cell codes plus a SEPARATE RUN TABLE.**
The four causes are not the same kind of object, and collapsing them into one per-cell
code would bake an unvalidated threshold into the stored data.

- **Per-cell codes**, for the causes determinable from the record itself:
  `0` called · `1` not applicable (nested child whose parent allele this haplotype does
  not carry — from `LV`/`PS`) · `2` inherited (the parent record's own call was missing)
  · `3` conflict (the `CONFLICT` tag names this sample) · `4` uncategorised no-call.
- **A separate run table**, `(sample, strand, start_row, end_row)`, listing every run of
  consecutive `-1`. Assembly gaps are NOT given a per-cell code.

**Why the split — measured on chr21, 2026-09-29, over all 1,278,653 `-1` entries:**

| test applied per cell | cells | share |
|---|---|---|
| in a run of >= 100 consecutive `-1` | 1,056,726 | 82.64% |
| at a `CONFLICT` record | 54,626 | 4.27% |
| at a nested record (`LV` > 0) | 624,261 | 48.82% |

Those exceed 100% because they overlap heavily: **long-run AND nested = 470,679 cells,
37% of all missingness**; long-run and conflict 39,095; conflict and nested 25,596; all
three 20,001; and **58,409 (4.57%) match none of the three**, which is why code 4 exists.

⚠ **Reason 1 — assembly gap is a property of a RUN, not of a cell.** The other three are
readable off the record. Classifying assembly gaps per cell means choosing a length
threshold, and the threshold is load-bearing: **>= 10 -> 92.06%**, >= 50 -> 85.85%,
>= 100 -> 82.64%, >= 500 -> 75.38%, **>= 1000 -> 68.67%** of all `-1` (measured). A
23-point swing on an arbitrary number. `docs/memory.md` is explicit that the >= 100
figure was *the length at which runs were counted, not a validated decision rule*.
Keeping runs as intervals lets any threshold be applied, changed or swept downstream
without re-reading the VCF — which is the same argument that justified the array itself.

⚠ **Reason 2 — precedence would otherwise be silent.** 470,679 cells are both inside a
long run and at a nested record, and the two call for OPPOSITE handling: not-applicable
means the DNA genuinely is not on that chromosome and RENORMALISE is correct, while an
assembly gap means the donor almost certainly carries an allele we cannot see. One code
per cell forces a winner; the split lets both facts be recorded and the conflict resolved
by a policy that can be stated and changed.

⚠ **PRECEDENCE AMONG THE PER-CELL CODES** must still be declared, because conflict and
nested also co-occur (25,596 cells). Pinned order: **3 conflict > 1 not applicable >
2 inherited > 4 uncategorised.** Conflict wins because the `CONFLICT` tag is an explicit
statement by the caller that it could not choose, which is strictly more informative than
the structural inference.

⚠ **"Inherited" costs a second pass.** It needs the PARENT's call — `PS` -> parent row ->
that sample's genotype there — so it cannot be decided while streaming the child's line.
Cheap, but it shapes the converter: build the matrix and the `PS` index first, then fill
code 2 in a second pass.

### anchor_snp  [implementation] [dataset]
```
Is:                   A pangenome VCF record that also appears in the PanGenie
                      "bi_all" callset, matched exactly on (POS, REF, whole-ALT-string).
                      PanMixer's within-block HMM runs over these and only these.
Computed by:          starting_data/scripts/src/get_mappings.py:45-49 ->
                      pangenome_to_thousand_g_alignments.pickle, keyed by pangenome ROW INDEX
CHANGES MEANING WITH: the callset used for the match. NOT the paper's "top level
                      SNPs V_SNPs", and not necessarily a SNP — the only claim that
                      these are SNPs is a source comment (hmm.py:32).
Valid vs:             nothing; it is a site set, not a measurement.
Status:               SECONDARY
Aliases:              "top-level SNP" — FORBIDDEN. The paper's V_SNPs and the code's
                      anchors are different sets; conflating them was one of this
                      project's first misreadings.
```

### eps_pmi  [privacy] [evaluation]
```
Is:                   PanMixer's per-block privacy score, eps_j = -log p(h_bj): the
                      self-information of the TARGET's ORIGINAL block under the cohort
                      model. Higher = rarer = more identifying.
Computed by:          external/PanMixer/tools/panmixer/obfuscate.py:230-251
CHANGES MEANING WITH: (a) whether the scoring panel includes the target — as shipped
                      it DOES, capping eps_j near log(2N) = 4.4773; (b) which branch
                      the block took, and the thresholds DISAGREE: scoring uses the
                      HMM at >=1 anchor, sampling at >=2, so 452 chr21 blocks are
                      scored one way and sampled the other; (c) MOST IMPORTANTLY, on
                      the HMM branch the forward algorithm scores ONLY the anchor
                      variants — non-anchor variants contribute nothing at all. In
                      32% of blocks the ignored non-anchor rarity EXCEEDS eps_j
                      itself. So eps_j is the self-information of the target's ANCHOR
                      alleles, not of its block.
Valid vs:             other eps_pmi values from the identical panel. **NEVER against
                      tau** — see the note below.
Status:               DIAGNOSTIC — an axis we can plot on, never an input to our mechanism.
Aliases:              "PMI", "privacy risk". Never write "privacy budget".
```

⚠ **`eps_pmi` and `tau` are formally incomparable.** `eps_pmi` is an average-case
self-information defined only for PanMixer's construction; `tau` is a worst-case
pairwise TV bound. There is no conversion. A head-to-head may align the two arms
only on MEASURED empirical attack success, never on these internal parameters.
Using `eps_pmi` anywhere inside `baseline_model` or `u_path` is FATAL — it is a
function of the private input's own path.

### eps_private  [privacy] [evaluation]
```
Is:                   PanMixer's operating threshold: "the LARGEST privacy risk value
                      at which all obfuscated target individuals could no longer be
                      linked" (published, page 4 lines 509-512). Published value 0.002.
Computed by:          swept empirically against ONE linkage attacker; not derived.
CHANGES MEANING WITH: the attack database (a larger or more genetically similar panel
                      demands a lower threshold — the paper says so explicitly), the
                      graph, and the target set. It is NOT a transferable constant.
Valid vs:             another eps_pmi on the identical panel and graph.
Status:               DIAGNOSTIC — a reproduction target, never an input to our mechanism.
Aliases:              "the privacy threshold". ⚠ The PREPRINT gives 0.001 and defines it
                      as the SMALLEST such value; both the number and the direction were
                      corrected in the published version. Cite 0.002 and "largest".
```

### utility_loss_panmixer  [utility] [evaluation]
```
Is:                   PanMixer's per-block cost of an edit, eta_j = sum over ALTERED
                      variants of 1/support(v), where support(v) counts non-missing
                      haplotype entries at v. Supplementary Remark 1 proves it equals
                      the induced L1 allele-frequency change.
Computed by:          external/PanMixer/tools/panmixer/obfuscate.py:207-211,279-290
CHANGES MEANING WITH: whether support is computed on the cohort INCLUDING the target
                      (as shipped it is); whether -1 is masked before comparing; and
                      a legacy 2/support variant in get_support.py:25.
                      ⚠ AND WITH THE ARM, via that first clause. In ARM 1 the target is
                      a cohort member and `loo_cohort` removes BOTH its haplotypes, so
                      support(v) is over K = 86 on our side while PanMixer's shipped
                      support is over all 88. Report which support each column used; a
                      head-to-head that silently mixes the two is comparing
                      denominators, not mechanisms.
Valid vs:             itself across capacities on one cohort, at a stated support. It
                      is a COST, unbounded above — it is NOT a similarity and must not
                      be used as phi_v without the normalization recorded under phi_v.
Status:               ARM 1: PRIMARY — **the ONLY metric adopted from PanMixer's
                      utility suite** (decision H, 2026-10-08); `af_loss` and `ld_loss`
                      are not adopted.
                      ✅ **AND IT IS SCORED AT K = 88 ON BOTH SIDES (decision M,
                      2026-10-12).** The head-to-head runs PanMixer's evaluator
                      UNMODIFIED, so this metric uses its shipped support INCLUDING
                      the target, for our released arrays and theirs alike. That puts
                      the same denominator on both sides, which is what the warning
                      above asks for. ⚠ So this column is NOT a leave-one-out
                      quantity, while our RELEASE is (K = 86, decision E) and every
                      quantity we compute ourselves stays at K = 86. Never put a
                      K = 88-scored column beside a K = 86 one. ARM 2: DIAGNOSTIC — an axis we can plot on, not
                      an endpoint.
Aliases:              "utility loss". Ours is `utility_retained`; do not mix the names.
```

### af_loss  [utility] [evaluation]
```
Is:                   PanMixer's allele-frequency utility loss U_AF: the L1 change in
                      cohort allele frequencies induced by a release, reported overall
                      and per MAF stratum.
Computed by:          external/PanMixer/tools/downstream/utility_in/af_loss.py:77-183
CHANGES MEANING WITH: the MAF strata, which are read from pangenome_mask.npy — the
                      cross-build match set. Under the GRCh37/GRCh38 mismatch every
                      stratum is 0/0; `docs/data.md` section 6 marks everything derived
                      from that mask DEPRECATED, and decision F requires it rebuilt on
                      OUR consistent GRCh38 build, genome-wide, before the quantity
                      could be computed at all.
                      Also with the row count: comparing a 44-row cohort AF to a 45-row
                      one shifts every frequency by ~1/45. ⚠ That caveat was a property
                      of the EXTERNAL-target reading — pretending the target is a 45th
                      member. In ARM 1 the release REPLACES one of the 44 columns, so
                      the row count is 44 on both sides and the shift does not arise;
                      the identical-site-set and identical-row-count requirement still
                      has to be ASSERTED rather than assumed.
Valid vs:             the other arm at the identical site set and row count.
Status:               PANMIXER-ONLY — **NOT ADOPTED as an endpoint (decision H,
                      2026-10-08), and that is a SCOPE CHOICE, not a definitional
                      impossibility.**
                      ⚠⚠ **MEASURED ANYWAY ON 2026-10-15, EXPLORATORY, AT DYLAN'S
                      REQUEST — AFTER OUTCOMES WERE INSPECTED.** ARM 1, chr21, 44
                      targets x 7 `tau` (`exp_5`) plus PanMixer's own 8 capacities
                      (`exp_6`), their `af_loss.py` run unmodified. Ours **0.003668 ..
                      0.003698** across the grid; theirs **0.000220 .. 0.002734**.
                      **THEY BEAT US AT EVERY CAPACITY** — ⚠ and `af_loss` is close to
                      their optimiser's OWN objective, since `utility_loss_panmixer`
                      equals the induced L1 AF change (their Supplementary Remark 1),
                      so this scores them on what they minimise and us on what we do
                      not optimise.
                      ⚠⚠ **IT FAILS THE DIRECTION CHECK THE PRE-REGISTERED ENDPOINTS
                      PASS: 20/44 targets down, 24/44 up** from `tau` 0 to 0.99, median
                      +6.478e-06. `target_fidelity` and `e_q_u` move 44/44 together.
                      So its movement is NOISE, not a small signal, and no `tau` trend
                      may be read from it. First metric in the project to fail that
                      check.
                      ⚠ Reporting it as an ENDPOINT needs a NEW dated pre-registration,
                      not an addendum. Provenance `logs/B_arm1_af_ld.md`.
                      ⚠ **THE OLD REASON WAS WRONG AND IS CORRECTED HERE.** This entry
                      used to say "we never modify the cohort, so the quantity is not
                      defined for our release object". That is true of ARM 2 only. In
                      ARM 1 the released object IS the 44 cohort members' obfuscated
                      paths, so a release does change the released cohort's allele
                      frequencies and U_AF is FULLY DEFINED for us there. The reason we
                      do not report it is that **Arm 1 reports
                      `utility_loss_panmixer` ONLY** from PanMixer's utility suite
                      (decision H); af_loss and ld_loss are out of scope for the
                      head-to-head, not out of reach.
                      ARM 2: n/a BY DEFINITION — the release never edits the cohort.
                      Our own utility axes remain `target_fidelity` (always against
                      its OWN ARM's ceiling), `utility_retained` and, as the
                      pre-registered endpoint, `e_q_u`.
Aliases:              "U_AF"
```

### donor_coherence  [utility] [evaluation]
```
Is:                   The mean, over random windows of W consecutive positions, of the
                      largest share of a window explained by any SINGLE donor
                      haplotype. It asks whether a release is a MOSAIC of real
                      haplotypes or a patched copy of the target.
                      ⚠ IT IS NOT A PRIVACY METRIC AND NOT AN ENDPOINT. High coherence
                      means the release looks like a real haplotype; it says nothing
                      about WHOSE. A `tau` = 0 release is coherent and carries ZERO
                      target information. It exists to explain the `ld_loss` gap
                      between the two tools, nothing more.
Computed by:          scripts/arm1_donor_coherence.py. W = 200, 2,000 windows, seed 0,
                      both haplotype columns, windows needing >= 50% called positions.
CHANGES MEANING WITH: W above all -- a shorter window is explained better by chance and
                      a longer one worse, so W must be stated with the number.
                      And with THE DONOR POOL: the target MUST be excluded (43 others
                      x 2 = 86 haplotypes). Leaving it in lets a patched copy of the
                      target match itself and score ~100%, measuring nothing.
Valid vs:             the TARGET'S OWN haplotypes at the same W and pool, never against
                      100%. ⚠ A real human genome scores only **98.3934%** at W = 200
                      against 86 others; that shortfall IS human haplotype sharing, and
                      scoring a release against 100% makes every release look worse
                      than it is.
Status:               REPORTED, ARM 1, as an EXPLANATION of `ld_loss` rather than a
                      result of its own. Measured 2026-10-15, chr21, HG00438:
                      **ours 94.83%..95.07% across the `tau` grid; theirs 97.88% at
                      `capacity` 0.01 falling to 92.17% at saturation**, against the
                      real genome's 98.39%.
                      ✅ **AT SATURATION OURS IS THE MORE HAPLOTYPE-LIKE OBJECT** --
                      about half their departure from a real genome -- which tracks
                      `ld_loss` (ours 0.00542, theirs 0.0401). Ours is a path through
                      real donor segments; theirs is the target's sequence with ~16% of
                      cells flipped, with nothing constraining whether the result is a
                      combination that occurs in nature.
                      ⚠ It also answers "is 16% altered just randomising?" -- no, for
                      either tool. A mosaic of real haplotypes naturally differs from
                      the target at ~16-17% of called cells; our `tau` = 0 release
                      alters 17.4% carrying zero target information. 16% is the FLOOR
                      set by haplotype sharing, not evidence of randomisation.
                      Provenance `logs/B_arm1_donor_coherence.md`.
Aliases:              none. ⚠ NOT a measure of similarity to the TARGET -- that is
                      `target_fidelity`, and coherence is deliberately blind to it.
```

### obfuscated_accuracy  [utility] [evaluation]
```
Is:                   PanMixer's own measure of how much of the TARGET a release still
                      carries, computed by their Beagle evaluator: the share of
                      xsol-selected called cells where the RELEASE equals the TRUTH.
                      Higher means LESS obfuscated. It is the complement of
                      `target_fidelity` restricted to the same cell mask, and it is
                      THEIR field name, not ours.
Computed by:          external/PanMixer/tools/beagle/accuracy_stats.py:compute_accuracy,
                      reported per variant subset alongside `beagle_accuracy`,
                      `major_allele_accuracy` and `reconstruction_accuracy`.
CHANGES MEANING WITH: the xsol mask, which decides WHICH CELLS are scored at all --
                      ours is an all-ones substitution, PanMixer's is their optimiser's
                      real block selection. ⚠ Never compare it across two different
                      masks; see `docs/bugs.md`.
                      And with the variant subset: the figures quoted are `snp_only`.
Valid vs:             the same tool across its own knob, and across tools ONLY at a
                      matched mask and matched subset.
Status:               REPORTED, ARM 1, as CONTEXT for the Beagle result rather than as
                      an endpoint of ours. **Measured 2026-10-15 on PanMixer's own 8
                      capacity releases: 0.998846 at `capacity` 0.01 falling to
                      0.854683 at saturation.** ⚠⚠ **IT IS WHAT EXPLAINS THEIR BEAGLE
                      CROSSOVER**: at `capacity` 0.01 their release is 99.88% identical
                      to the truth and alters 597 SNP cells, so there is nearly nothing
                      for the attack to undo and a baseline that ignores the data is
                      easy to beat. Provenance `logs/B_arm1_pm_beagle.md`.
Aliases:              none. ⚠ NOT `target_fidelity`, which is OURS, is a DIFFER rate
                      (lower is better) and is computed over a different cell set.
```

### ld_loss  [utility] [evaluation]
```
Is:                   PanMixer's linkage-disequilibrium utility loss U_LD: the change
                      in pairwise LD among nearby sites induced by a release.
Computed by:          external/PanMixer/tools/downstream/utility_in/ld_loss.py:16-123
CHANGES MEANING WITH: the window — the shipped KB_WINDOW_SIZE = 5
                      (external/PanMixer/tools/downstream/utility_in/ld_loss.py:8) is
                      compared against raw bp at ld_loss.py:107, so the window is
                      5 BASE PAIRS, not the paper's 5 kb; and the site set
                      (pangenome_mask.npy again, DEPRECATED in `docs/data.md`
                      section 6), which under the build mismatch makes it return
                      (0.0, 0) rather than raising.
Valid vs:             the other arm at the identical window and site set. If it is ever
                      run, run it BOTH as-shipped (5 bp, for parity with the published
                      figure) and corrected to 5,000 bp, as two separate rows — never
                      one row whose window is unstated.
Status:               PANMIXER-ONLY — **NOT ADOPTED as an endpoint (decision H,
                      2026-10-08).**
                      ⚠⚠ **MEASURED ANYWAY ON 2026-10-15, EXPLORATORY, AFTER OUTCOMES
                      WERE INSPECTED.** ARM 1, chr21, 44 targets x 7 `tau` (`exp_5`)
                      plus PanMixer's own 8 capacities (`exp_6`), their `ld_loss.py`
                      run unmodified. Ours **0.005228 .. 0.005418** per pair, flat
                      across `tau`; theirs **0.000042** at `capacity` 0.01 then
                      **0.0380 .. 0.0401** from 0.05 up. **WE BEAT THEM FROM
                      `capacity` 0.05 UP BY ~7.4x**; at 0.01 they are ~130x better,
                      where they alter 7,687 cells against our ~107,797.
                      ⚠⚠ **THE WINDOW IS 5 POSITION UNITS, NOT 5 KB.**
                      `KB_WINDOW_SIZE = 5` (`ld_loss.py:8`) is compared against a RAW
                      position difference (`ld_loss.py:107`). Median SNP gap on chr21
                      is 83 bp, so it reaches **2.92% of SNP positions and ~7,450
                      pairs**; a true 5 kb window would be **10,588,565 pairs**, ~1,400x
                      more. We run it AS SHIPPED — that is the faithful comparison —
                      and the window must be stated wherever the number appears,
                      because the name implies something 1,400x larger than what is
                      measured.
                      ⚠ Our 44-target spread MIXES between-person variation with
                      single-draw noise (one release per target per `tau`,
                      constraint 5), so the grid scatter is NOT a trend.
                      Provenance `logs/B_arm1_af_ld.md`.
```

### capacity  [mechanism] [evaluation]
```
Is:                   PanMixer's knapsack budget: the per-chromosome fraction of the
                      maximum achievable utility loss the optimizer may spend.
Computed by:          external/PanMixer/tools/panmixer/obfuscate.py:453-455
CHANGES MEANING WITH: the chromosome (it is applied per chromosome, 22 times
                      independently), and the cohort that sets the maximum.
Valid vs:             another capacity on the same chromosome and cohort. **Not
                      comparable to tau**, which is genome-wide and a different quantity.
Status:               DIAGNOSTIC
Aliases:              none
```

### new_haplotypes  [implementation] [dataset]
```
Is:                   PanMixer's integration artifact: the released allele vector for one
                      target, shape (n_sites, 2) int16, -1 = missing, row-aligned to
                      that chromosome's pangenome VCF record order. Every PanMixer
                      attack and utility evaluator consumes exactly this file.
Computed by:          PanMixer writes it from tools/panmixer/stacker.py. ✅ **ON OUR
                      SIDE, scripts/emit_new_haplotypes.py since 2026-10-13.** This read
                      "PLANNED and UNIMPLEMENTED — nothing in src/ emits this format".
                      It maps our release from OUR per-target chain onto THEIR site axis
                      by POSITION, and ships a sites.tsv beside every array.
                      ⚠⚠ ROWS WE DO NOT RELEASE ARE -1, NOT THE TARGET'S TRUE ALLELE.
                      Their stacker leaves the original in place outside the blocks it
                      swaps, so filling ours the same way would hand a linkage attack
                      ~9,400 alleles our mechanism never released. Measured for HG00438:
                      331,416 released rows of their 340,824, 9,408 at -1.
                      ⚠ Their evaluators mask on both sides being called, so those rows
                      DROP OUT rather than scoring as matches — the conservative
                      direction, and it must be stated with any number from this file.
                      ✅ THE ALIGNMENT CONTROL, which the vacuous shape assertion below
                      cannot provide: map the target's own original onto their axis and
                      compare to THEIR original. Correct alignment gives **98.6475%**
                      jointly-called agreement; shifts of 1, 2, 10, 100, 1,000 and 10,000
                      rows all give 71.9-77.0%, every one caught by the 0.90 threshold.
                      The script REFUSES to write below it.
                      ⚠ It is no longer optional: under decision E Arm 1 runs PanMixer's
                      own evaluations on our release, and this file is what every one of
                      their evaluators consumes. It is ARM 1 PLUMBING.
CHANGES MEANING WITH: the site axis. It is positional — a row-order mismatch silently
                      scores the wrong variants rather than erroring. ⚠ The shape
                      assertion in PanMixer's gap-score evaluators is VACUOUS (both sides
                      equal site_mask.sum() by construction), so nothing upstream catches
                      this. Ship a `sites.tsv` beside any emitted copy and verify it.
                      ⚠ AND WITH THE BUILD. Decision F rebuilds PanMixer's starting data
                      on OUR consistent GRCh38 build and runs BOTH arms on it, so a
                      new_haplotypes written against their shipped axis and one written
                      against ours are not the same object. Name the build.
                      ⚠ AND WITH THE ARM, in extent: Arm 1 emits all 44 targets, DIPLOID,
                      over the 22 AUTOSOMES (decision D) — 44 files per chromosome, one
                      per target, each written under that target's own K = 86 LOO panel.
Valid vs:             another new_haplotypes over the identical site axis, build and
                      cohort.
Status:               ARM 1: REQUIRED — it is the interface the head-to-head runs
                      through, and emitting it is the cheapest way to run PanMixer's
                      ATTACK evaluators on our release unmodified. ARM 2: OPTIONAL — not
                      a requirement of the mechanism.
                      ⚠ **CORRECTED 2026-10-08:** this entry used to read
                      `COMPARISON-ARM ONLY ... NOT a requirement`, and to say that
                      emitting it "does not make PanMixer's cohort-level utility metrics
                      applicable to us — `af_loss` and `ld_loss` would run and return
                      meaningless numbers". In ARM 1 those metrics ARE defined on this
                      object (the release edits a cohort column); we decline them as a
                      SCOPE CHOICE under decision H, reporting
                      `utility_loss_panmixer` only. See `af_loss` and `ld_loss`.
Aliases:              "the release", "the obfuscated haplotypes"
```

### loo_cohort  [dataset] [method]
```
Is:                   The PER-TARGET leave-one-out cohort D = HPRC \ {target}, with
                      every derived artifact recomputed without the target. ✅ **THE
                      CORE CONSTRUCTION OF ARM 1 (decision E, 2026-10-08):** BOTH of
                      the target's haplotypes leave the donor panel, so **K = 86**, and
                      it is applied PER TARGET — 44 times, once per cohort member.
                      What must be recomputed without the target, under model C:
                      the donor panel (K = 86), support(v) and therefore w_v = 1/support
                      and `beta_v`, the allowed-state mask and `output_support` Y_D, the
                      baseline `R_D` with its `baseline_log_norm`, `u_ceiling`, and
                      `expected_switches`. ⚠ NOT allele frequencies and NOT a block
                      assignment: model C has neither — there is one chain position per
                      VARIANT and no block object in the mechanism.
                      ✅ **HOW IT IS APPLIED (settled, 2026-10-08):** by MASKING the
                      target's two haplotype ROWS in memory at release time, not by
                      materialising 44 cohort trees. The artifacts are identical and it
                      is ~160-380 GB cheaper.
                      ✅ The attack database is built ONCE, with the target removed AT
                      USE rather than 44 databases built.
Computed by:          src/cohort/loo.py:loo_chain (built 2026-10-12). This read
                      "PLANNED — the row mask is not implemented". It masks both rows in
                      memory, recomputes support(v) over the remaining donors, and
                      RE-APPLIES decision K per target. Surveyed over all 44 x 22 = 968
                      chains by scripts/loo_chain_survey.py.
                      ⚠ Validated by reconstruction, not by inspection: run at K = 88
                      the same composition reproduces in_chain.npy BIT-IDENTICALLY and
                      support.npy bit-identically, so the only difference on the
                      leave-one-out path is the two masked rows.
                      ⚠⚠ **MEASURED, AND THE UNTRAVERSABLE SET IS NOT MONOTONE IN THE
                      DONOR PANEL.** T spans 22,129,531 to 22,136,532 over the 44
                      targets against arm 2's 22,134,562, and `untraversable_block`
                      counts span 19 to 21 against arm 2's 20 — so leave-one-out
                      sometimes leaves FEWER dead ends, not more. On chr16 the 2,026-
                      position block arm 2 excludes becomes TRAVERSABLE once HG02630 is
                      held out, putting 1,981 positions in that target's chain that are
                      NOT in arm 2's. Mechanism: dropping positions SHORTENS the
                      surviving zero-distance runs, and a shorter run needs a donor to
                      span fewer positions. **So arm 1's chain must never be built by
                      subsetting arm 2's.**
CHANGES MEANING WITH: which artifacts were actually recomputed. Recomputing only the
                      donor panel (what PanMixer does) is NOT a loo_cohort — and under
                      decision E, every item in the list above must be recomputed per
                      target or standing constraint 1 is violated.
                      ⚠ AND WITH THE ARM. It is ARM 1's construction and ARM 2 has no
                      instance: there the target is outside `G`, `D`, every panel and
                      the attack database, so nothing is held out and `K = 88`.
                      ⚠ AND IT GIVES EACH ARM-1 TARGET ITS OWN SETTINGS — its own D, its
                      own `beta_v`, its own Y_D, its own ceiling. So the 44 releases
                      carry 44 utility functions: pool measured attack outcomes, never
                      utilities (constraint 8), and see `joint_release_gap`.
Valid vs:             the same held-out individual under the same recompute list. ⚠ Not
                      against ARM 2: a K = 86 LOO release and a K = 88 external release
                      are different mechanisms on different panels.
Status:               ARM 1: PRIMARY — the construction the whole arm rests on, and what
                      makes standing constraint 1 hold there at all. ARM 2: n/a — there
                      is nothing to hold out.
                      ⚠ **CORRECTED 2026-10-08.** This entry used to read
                      `COMPARISON-ARM ONLY ... Not a concept our mechanism needs`, with
                      "`K = 88`, unconditionally" and "externality must therefore be
                      SIMULATED". All three are now false for half the paper. Arm 1 does
                      NOT simulate externality — the target genuinely IS a cohort member,
                      and that is the point of the arm, since what is obfuscated is a
                      released donor path. The LOO survives for a different reason:
                      keeping Y_D, the support counts and `beta_v` out of the protected
                      individual's own data, which is standing constraint 1.
Aliases:              "LOO". ⚠ Even a full loo_cohort is an APPROXIMATION of the
                      external-target setting: the PGGB graph TOPOLOGY was built from all
                      HPRC assemblies, so G itself saw the stand-in. ⚠ That is now an
                      ACCEPTED LIMITATION OF THE PAPER, stated in Arm 1 rather than a
                      caveat about an optional arm: Arm 1's guarantee is CONDITIONAL ON
                      THE GRAPH. State it; never claim it away.
```

### target_fidelity  [utility] [evaluation]
```
Is:                   **The share of the target's CALLED alleles that the release
                      ALTERS — a DIFFER rate, where LOWER is more faithful.** It is the
                      axis that transfers between the two arms, as opposed to
                      cohort-level utility, which a pure prior draw satisfies perfectly
                      while carrying zero target information.
                      ⚠ **DIRECTION CORRECTED 2026-10-08, and it was an inversion.**
                      This field used to read "how much of the TRUE target the release
                      retains", which is the COMPLEMENT of what every script and every
                      reported number in this repo computes: `docs/plan.md` calls it
                      "the share of the target's called alleles the release alters" and
                      quotes ours at 15.03%; `docs/data.md` section 7 reports our
                      frontier as 14.8077%-14.8070% against PanMixer's 1.20%-16.41% on
                      that same convention; `docs/memory.md` 2026-10-03 tabulates it as
                      "cells differing from the target". Under the old wording 15.03%
                      would have meant the release retains 15% of the target, when it
                      in fact retains 84.97%. The DIFFER rate is kept as the definition
                      because it is what is measured; the identifier is unchanged so
                      existing table rows still resolve. ⚠ The NAME still reads like a
                      retention, which is a rename decision for Dylan, not a silent fix.
Computed by:          external/PanMixer/tools/downstream/privacy/MIA_privacy.py:80-83,140
                      already computes match fractions of a release against the target's
                      ORIGINAL haplotypes, mechanism-agnostically. ⚠ Those are MATCH
                      fractions — take the complement to get this quantity, and say
                      which you printed.
CHANGES MEANING WITH: the ARM (see the two ceilings below), the site set it is measured
                      over, and whether a half-missing genotype is treated as hom-ref
                      (MIA_privacy.py:49-51 does).
                      ⚠ AND WITH WHAT COUNTS AS THE TARGET'S OWN PATH IN ARM 1
                      (decision I, 2026-10-08). Arm 1 scores "self" from the target's
                      OWN ASSEMBLY-DERIVED path, which is defined for all 44 cohort
                      members — not from its 1000G panel row, which exists for only 39
                      of the 44. PanMixer's panel-row semantics run as a separate,
                      LABELLED PARITY COLUMN, never pooled with ours.
⚠ CEILING (ARM 2):    In the EXTERNAL-target setting this quantity has a maximum BELOW 1
                      that does not depend on `tau`. The release is a path through `G`, so
                      a target variant with no record in `G` cannot be retained at any
                      `tau`, including `tau` -> 1. Measured 2026-09-27 on six 1000G samples
                      outside the HPRC 44: about **9.9%** of a target's chr21 non-reference
                      calls have NO record in the graph at all, and those carry **20.5%** of
                      the target's total `-log f` information (2.07x over-representation).
                      A further 9.9% have a record at that POS under a different REF/ALT
                      spelling; those ARE measured, and should be reported separately since
                      `phi_t` charges the mechanism for a join artifact. Ceiling on
                      retainable information: about **79.5%**. **Never report
                      target_fidelity without its ceiling**, or a perfect result is
                      indistinguishable from a mediocre one. Note also that `u(p, y)` is
                      defined over chain positions, which come from VCF rows, so u can
                      reach 1.0 while a fifth of the target's variation was dropped before
                      the mechanism ran. u is fidelity-on-what-the-graph-can-see.
                      ⚠ EVERYTHING IN THIS BLOCK IS ARM-2 REASONING and says so in its
                      first clause. Carrying the 79.5% into ARM 1 would understate the
                      mechanism by about 20 points on the one axis the two arms are
                      aligned on.
⚠ CEILING (ARM 1):    A DIFFERENT ceiling, from a different cause, and **UNMEASURED.**
                      In Arm 1 the target's own assembly is one of the 44 the graph was
                      built from, so essentially every target variant HAS a record in
                      `G` and the graph-representability limit above does not bind. What
                      binds instead is DONOR-PANEL REACHABILITY under `loo_cohort`: at
                      K = 86, with both of the target's own haplotypes removed, an
                      allele no remaining donor carries cannot be emitted at any `tau`.
                      It must be computed ON THE LOO PANEL, per target, and it CANNOT be
                      derived from the Arm-2 figure — see `u_ceiling`, whose own warning
                      is about exactly this conflation, and `private_allele_survival`,
                      which is the extreme case of it.
Valid vs:             the other arm at matched empirical privacy, and only against a
                      stated ceiling computed on the same graph and **in the same arm**.
                      ⚠ The rule is now ARM-INDEXED: never report target_fidelity
                      without ITS OWN ARM's ceiling, and never put the two arms'
                      ceilings in one column.
Status:               PRIMARY — without this axis our tau = 0 point looks like a free lunch.
Aliases:              none. ⚠ **"retained" is FORBIDDEN for this quantity** — it is a
                      DIFFER rate, and `utility_retained` is a different quantity in the
                      opposite direction. "fidelity loss" and "the differ rate" are the
                      safe informal names.
```

---

## Measured quantities

These are estimates from finite samples. The theory quantities above are exact;
these are not, and every one of them needs an `n` and a `tau`.

### utility_retained  [utility] [evaluation]
```
Is:                   u_path(p_g, released_path) achieved by an actual release, for
                      the true target path p_g. Averaged over independent draws. It is
                      the FINITE-SAMPLE counterpart of `e_q_u`.
                      ⚠ **DEMOTED BY DECISION A (2026-10-08).** The primary
                      pre-registered endpoint is the genome-wide EXACT expectation
                      `e_q_u`, reported from derivatives of log Z_p or from state
                      marginals with NO sampling. This sampled average is a diagnostic:
                      per decision G, 3 seeds on a single configuration, as a
                      diagnostic only. The paper's numbers do not rest on it.
Computed by:          src/mechanism/ffbs.py:Release.u_achieved per draw; averaged and
                      tabulated by scripts/run_mechanism_chr21.py (built 2026-09-30).
                      ⚠ Never quotable without u_ceiling beside it.
CHANGES MEANING WITH: tau (monotonically), the utility function, AND the denominator
                      (see the two measurement paths below). Report with n draws.
Valid vs:             the same utility function at another tau, against the floor
                      (tau = 0, an untilted draw from R_D) and the ceiling.
                      ⚠ CORRECTED 2026-09-30: the ceiling is NOT "the target's own path,
                      u = 1 by construction" — that holds only for a target with a call at
                      every weighted position. For a real external target it is u_ceiling,
                      measured 0.815455 on HG00096 chr21 at beta scope "chain". The old
                      wording is what produced the 0.6940-against-1 misreading recorded in
                      docs/bugs.md; see u_ceiling.
                      ⚠ AND THE CEILING IS NOW ARM-INDEXED: Arm 2's is `u_ceiling`
                      measured on graph readability, Arm 1's is set by the K = 86 LOO
                      panel and is UNMEASURED. Read each against its own arm's.
Status:               DIAGNOSTIC — demoted from PRIMARY by decision A, 2026-10-08. The
                      PRIMARY utility endpoint is `e_q_u`, exact and genome-wide.
Aliases:              "utility". Never quoted bare without its tau, its denominator, its
                      arm or its u_ceiling. ⚠ Not `target_fidelity`, which is a DIFFER
                      rate in the opposite direction.
```

### delta  [attack] [evaluation] [gate]
```
Is:                   The MINIMUM DETECTABLE ADVANTAGE: the smallest excess over chance
                      that the pre-registered attack experiment could distinguish from
                      chance at one-sided 95% with 80% power. It is a property of the
                      EXPERIMENT's size, not of the mechanism, and it is what makes a
                      null statement quantitative instead of an absence.
                      ⚠⚠ WITHOUT IT A NULL IS NOT EVIDENCE (standing constraint 7). An
                      attacker measured at chance says "accuracy < 0.5 + delta", and
                      that sentence is only worth something when delta is stated.
                      ✅ PINNED at delta = 0.056 (prereg field 6, signed 2026-10-10).
Computed by:          the design, in closed form. For the BALANCED estimator this
                      project signed, delta = 0.879 / sqrt(n) with n trials PER ARM; for
                      a one-sided rate it is 1.243 / sqrt(n). Both come from
                      (z_0.95 + z_0.80) * sd, with sd = 0.354/sqrt(n) balanced against
                      0.5/sqrt(n) one-sided, because averaging two arms halves the
                      variance. At n = 250 per arm, delta = 0.0556.
                      ⚠ The balanced form is MORE precise per observation, which is why
                      the corrected design needed 25 draws per target rather than 50 to
                      hold delta fixed. That is also how the formula was validated: the
                      one-sided branch reproduces the 1.243 constant already on record.
CHANGES MEANING WITH: the number of trials, and WHICH ESTIMATOR they feed. ⚠ A delta
                      quoted for a one-sided rate does NOT apply to a balanced accuracy
                      and vice versa — the constants differ by 1.41x.
Valid vs:             the same estimator and the same n. ⚠ NEVER compare a delta across
                      arms or tau values; it is a property of the experiment at one
                      configuration.
Status:               PRIMARY — pre-registered, and every null must be reported with it.
                      ⚠ IT DOES NOT COVER THE WHOLE tau GRID, and saying so is required:
                      at tau = 0.05 the bound permits 0.525 while the experiment can
                      only state "< 0.556", so the PROOF carries tau <= 0.10 and the
                      EXPERIMENT carries tau >= 0.25. Reporting a uniform null across
                      all seven points as seven pieces of evidence would overstate it.
Aliases:              "MDE", "minimum detectable effect" — acceptable in prose, but the
                      row label is `delta`. ⚠ NOT a p-value and NOT an effect size that
                      was measured; it is the resolution floor of the instrument.
```

### kl_exact  [privacy] [theory] [evaluation]
```
Is:                   KL(Q_p || Q_q) between the release distributions induced by two
                      input paths, computed EXACTLY with zero Monte-Carlo error.
                      ⚠⚠ IT IS A THEORY QUANTITY, NOT AN ATTACK OUTCOME, and decision L
                      (Dylan, 2026-10-11) makes that explicit. It is NOT a substitute
                      for the pre-registered attacker A2 and must never be reported as
                      one: A2 is exact ENUMERATION of TV over all admissible input
                      pairs, whose positive control already reads "already passes",
                      which is the toy check it was written to be. A2 stays satisfied by
                      the toy; this is its own named quantity at genome scale.
                      ⚠ Reporting it is not selective reporting, and the reason is
                      worth stating: it is DETERMINISTIC given the model and the input
                      pair, so there is exactly one value and nothing to choose among
                      after seeing an outcome. That is what makes it safe to introduce
                      after the pre-registration was signed; an attacker's measured
                      accuracy would not be.
Computed by:          scripts/exact_kl.py. Both E_{Q_p}[u_p] and E_{Q_p}[u_q] are LINEAR
                      functionals of Q_p's posterior state marginals
                      (ffbs.state_marginals), so each is one pass with no sampling and
                      no differencing, and log Z comes from ffbs.log_z_p.
                      ⚠ WHY IT EXISTS AT ALL: sampling cannot resolve this quantity.
                      MEASURED 2026-10-11, chr21, tau = 0.25, D1_extremal: sd(Lambda)
                      per draw is 6.15e-04 against a true value of 1.784e-07, so
                      resolving E[Lambda] by sampling needs n ~ 1.2e7 draws. The
                      pre-registered 250 trials per cell are five orders of magnitude
                      short. This route takes 54 s for all ten targets.
                      ⚠ VALIDATED against brute-force enumeration on the toy at four
                      tau values (tests/test_exact_kl.py), which is the only thing that
                      can check it: at chr21 scale the support cannot be enumerated.
                      That validation caught a real bug -- state_marginals takes
                      (hmm, fwd, potentials) and the first version passed (hmm, psi,
                      fwd), which does not raise because the shapes are compatible.
CHANGES MEANING WITH: the input PAIR, the arm, tau and the cohort. ⚠ KL is NOT
                      symmetric: KL(Q_p||Q_q) != KL(Q_q||Q_p). Always say which
                      direction, and the direction that matters here is from the
                      TARGET's release distribution, because that is the one the
                      attacker observes.
Valid vs:             the same pair, arm, tau and cohort.
Status:               PRIMARY for the Arm 2 privacy statement at genome scale.
                      MEASURED 2026-10-11, chr21, tau = 0.25, D1_extremal: +1.784e-07
                      for HG00096, and POSITIVE for all ten targets (1.38e-07 to
                      1.80e-07). ⚠ KL >= 0 is a theorem, so a negative value falsifies
                      the computation, not the mechanism.
Aliases:              none. ⚠ NOT "the TV" and NOT "A2" — see `acc_ceiling` for what it
                      licenses, and A2 for what it is not.
```

### acc_ceiling  [privacy] [theory] [evaluation] [gate]
```
Is:                   The largest accuracy ANY equal-prior binary attacker can achieve
                      against a given input pair, derived from `kl_exact`:
                          acc_ceiling = (1 + sqrt(kl_exact / 2)) / 2
                      via Pinsker's inequality TV <= sqrt(KL/2) and the Neyman-Pearson
                      identity that optimal accuracy = (1 + TV)/2.
                      ⚠⚠ IT IS STRICTLY STRONGER THAN A MEASURED ATTACK, and that is
                      the point: it covers attackers nobody has written, whereas
                      `attacker_accuracy` is one implemented attacker's result and is
                      only a LOWER bound on distinguishability (constraint 7).
                      ⚠⚠ AND IT IS LOOSE, WHICH MUST BE SAID EVERY TIME. Pinsker is an
                      inequality, not an identity, so the TRUE ceiling is closer to 0.5
                      than this number. Quoting it without that caveat overstates how
                      tight our knowledge is -- and the looseness runs in the direction
                      that flatters us, which is why it is stated rather than implied.
Computed by:          scripts/exact_kl.py, printed as its own row beside
                      `p_succ_bound`.
CHANGES MEANING WITH: everything `kl_exact` changes with. ⚠ Compare it to
                      `p_succ_bound` and nothing else: `p_succ_bound` = (1 + tau)/2 is
                      what the THEOREM promises, `acc_ceiling` is what this cohort and
                      this pair actually permit, and the gap between them is the honest
                      measure of how much slack the bound has at genome scale.
Valid vs:             the same pair, arm, tau and cohort.
Status:               PRIMARY, and a GATE: `acc_ceiling` must not exceed
                      `p_succ_bound`. If it does, either the computation or Theorem 1
                      is wrong, and the computation is the one to suspect.
Aliases:              "the Pinsker ceiling" is acceptable in prose. ⚠ NEVER "the
                      attacker accuracy" — that name belongs to a measured quantity
                      and means close to the opposite thing.
```

### attacker_accuracy  [attack] [evaluation] [gate]
```
Is:                   The empirical success rate of an IMPLEMENTED binary attacker
                      asked to decide which of two admissible INPUT PATHS a release
                      came from, at one release slot, under equal priors. ⚠ RESTATED
                      2026-10-08: the older wording — "came from the target or from a
                      cohort member" — is VACUOUS in ARM 1, where the target IS a
                      cohort member. In Arm 1 the pair is member k's own true diploid
                      path against an admissible decoy at slot k; in Arm 2 it is the
                      external target against any other admissible path. See
                      `p_succ_bound`.
Computed by:          src/attacks.py — BUILT 2026-10-10, and load-bearing: it holds
                      the Neyman-Pearson optimal statistic (lrt_statistic), the exact
                      ceiling (exact_advantage), the BALANCED empirical_accuracy, the
                      arm guard assert_well_posed, AttackSpec with a mandatory
                      positive_control field, and require_signed_prereg. The harness
                      is scripts/arm2_attack.py. ⚠ CORRECTED 2026-10-08: this field
                      once named src/ranpanmixer/attacks.py, a package layout this
                      repo does not use. [INCORRECT] - and it then said "PLANNED"
                      until 2026-10-10, three days after the module was written.
                      ⚠⚠ BOTH ERRORS SURVIVED FOR THE SAME REASON: the contract
                      test's _ABSENT_MARKERS skips any line containing "PLANNED", so
                      a line can claim a file is unbuilt forever without the suite
                      noticing. Recorded here on 2026-10-08 and STILL not fixed when
                      it bit a second time on 2026-10-10 — a known blind spot left
                      open is worse than an unknown one.
CHANGES MEANING WITH: the attacker. It is a LOWER bound on distinguishability and
                      nothing else. It cannot validate the mechanism; it can only
                      falsify an implementation.
                      ⚠ AND WITH THE ATTACKER SET, WHICH IS PRE-REGISTERED
                      (constraint 6) AND IS NOW WRITTEN DOWN —
                      paper/prereg_2026-10-10.md section 4 pins FOUR attackers for
                      both arms: A1 the exact likelihood-ratio distinguisher, A2
                      exact enumeration of TV and ratio_spread, A3 rarity-weighted
                      genotype linkage re-posed as a 1-of-2 discrimination, and A4 a
                      segment-restricted likelihood ratio. Each carries a chance
                      level and a MANDATORY positive control. Database search and
                      membership inference are ill posed in arm 2 and must RAISE
                      rather than return a number (src/attacks.py). Arm 2's target
                      count is settled at all 10, chr21 for the powered attack and
                      genome-wide for exact utility.
                      [INCORRECT] - The superseded text said the file "is PLANNED and
                      does not exist, so there is no pinned attacker set yet for
                      either arm, and Arm 2's attacker set and target count are
                      explicitly OPEN". ⚠ That mattered more than ordinary staleness:
                      this is a GATE whose meaning changes with the attacker set, so
                      telling a future session nothing is pinned invites it to choose
                      an attacker AFTER seeing an outcome — precisely what constraint
                      6 exists to prevent. What IS settled: ARM 1's linkage attack scores "self"
                      from the target's OWN assembly-derived path, defined for all 44
                      cohort members, NOT from its 1000G panel row, which exists for
                      only 39 of the 44 (decision I); PanMixer's panel-row semantics
                      run as a separate LABELLED PARITY COLUMN. The attack database is
                      built ONCE with the target removed AT USE.
                      ⚠ AND WITH THE ARM, for pooling: Arm 1's 44 releases have 44
                      different utility functions, so measured attack OUTCOMES pool
                      across targets and utilities never do (constraint 8).
Valid vs:             p_succ_bound at the SAME tau and in the SAME arm — as a one-sided
                      check that the implementation is not leaking more than the theory
                      allows. Never as evidence the mechanism is "more private than tau".
Status:               PRIMARY — mark * in tables.
Aliases:              "attack accuracy", "identification rate"
```

### tv_empirical  [privacy] [evaluation] [diagnostic]
```
Is:                   An ESTIMATE of TV(Q^tau_p, Q^tau_q) from finite samples.
Computed by:          src/mechanism/exact.py:tv on the toy support, by EXACT enumeration
                      (so on the toy it is not an estimate at all); at scale it is not
                      directly estimable and must be reported as such. Built 2026-09-30.
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

### two_target_indistinguishability  [privacy] [evaluation] [gate]
```
Is:                   The test our guarantee ACTUALLY makes, and the reason the
                      evaluation cannot simply be inherited. Run the mechanism on two
                      admissible INPUT PATHS p and q **at one FIXED release slot** —
                      one D, one R_D, one beta_v, one Y_D — draw R releases of each,
                      and estimate the distance between the two release distributions.
                      The theorem says it is <= tau for EVERY such pair; PanMixer's
                      design cannot bound it at all, because it releases unselected
                      blocks verbatim, so two inputs differing in an unselected block
                      have DISJOINT supports and TV = 1.
                      ⚠ **RESTATED 2026-10-08, and the old wording had no Arm-1
                      instance.** It read "two admissible TARGETS X and Y under
                      identical settings". Under decision E the 44 Arm-1 releases do
                      NOT share settings — per-target leave-one-out gives each its own
                      D, its own K = 86 panel, its own beta_v and its own output
                      support — so "two targets under identical settings" cannot be
                      instantiated in Arm 1 at all. The pair is a pair of INPUTS AT ONE
                      SLOT, which is what Theorem 1 quantifies over, and that version
                      exists in both arms.
Computed by:          tests/test_toy_enumeration.py, EXACTLY, on the toy support
                      (built 2026-09-30): all 3,571,128 pairs of admissible inputs per
                      tau, worst pair reported. ⚠ STILL PLANNED at chr21 scale, and see
                      docs/plan.md R5 before building it — the release barely varies with
                      tau under the current utility, so the measurement would be of a
                      mechanism that does not move. ⚠ **NOTHING in PanMixer's evaluation
                      suite tests this** (established 2026-09-28 by reading every
                      evaluator). Their metrics test their claim.
CHANGES MEANING WITH: the estimator and the support size — see tv_empirical: plug-in TV
                      over a large support is biased UPWARD, so a value above tau at
                      scale is the estimator misbehaving, not a violated theorem. Only
                      exact enumeration on a toy support can falsify Theorem 1.
                      Also with the PAIR: the bound is over all pairs, so a single
                      favourable pair proves nothing. Report the worst pair found.
                      ⚠ AND WITH WHERE THE PAIR IS DRAWN FROM (2026-10-08). Under
                      per-target leave-one-out the SETTINGS are per-target, so a pair
                      drawn ACROSS TWO Arm-1 members is a DIFFERENT and UNBOUNDED
                      quantity: it mixes two cohorts in one comparison, which standing
                      constraint 8 forbids and which tau does not bound. Draw both
                      members of the pair at one slot.
                      ✅ THE DECOY CONSTRUCTION IS PINNED (2026-10-10) in
                      paper/prereg_2026-10-10.md as THREE conditions, reported
                      SEPARATELY and never pooled: D1 extremal, D2 nearest by
                      genotype concordance (deterministic, so it takes no seed), and
                      D3 random. Decoy RNG seed 20261010 via
                      numpy.random.SeedSequence. ⚠ One release is tested against all
                      three, so the three are CORRELATED through the shared draw.
                      [INCORRECT] - The superseded text said the construction "is not
                      yet pinned" and forecast "the EXTREMAL pair PLUS a seeded random
                      sample" — the pinned form has THREE conditions, not two, and
                      the nearest-concordance condition it omits is the hard one that
                      carries the claim.
                      ⚠ The DRAW COUNT is NOT settled: field 6 was reopened
                      2026-10-10 because the design was one-sided and could not
                      measure the quantity p_succ_bound bounds. See that field.
Valid vs:             tau, at matched settings, on a stated support, at one release
                      slot.
Status:               PRIMARY — this is the headline privacy claim.
Aliases:              "the matched-pair experiment", cf. `A-ATK-matched-pair`
```

### joint_release_gap  [privacy] [theory] [evaluation]
```
Is:                   The COMPOSITION GAP in ARM 1: the paper reports 44 releases whose
                      donor panels are each other's true haplotypes, so the PER-RELEASE
                      tau bound holds for every one of them and the JOINT bound over
                      the 44 does NOT follow. Member j's released path is a mosaic of
                      the other 43 members' real haplotypes, so the 44 releases are
                      CORRELATED through the panel, and an attacker who reads all 44
                      together is not the attacker `p_succ_bound` bounds.
                      ⚠ It is a NAMED LIMITATION, not a defect to design around and not
                      a quantity to leave implicit. Standing constraint 5 says so:
                      report it; do not quietly design around it.
Computed by:          NOT COMPUTED, and no bound for it is claimed. What CAN be reported
                      is stated rather than derived: R = 1 release per (target, tau),
                      cached and reused (decision G and constraint 5), so no single
                      genome is re-released; the 3-seed run is a DIAGNOSTIC on one
                      configuration, never 3 releases of a person. ⚠ Any numeric joint
                      bound appearing for this quantity is unverified — none exists.
CHANGES MEANING WITH: the ARM. ARM 1 only: it exists because the released objects ARE
                      cohort members and the donors are each other. ARM 2 has no
                      instance — one external target, one release, the panel is public
                      and is not released by us.
                      ⚠ AND WITH R. It is a composition over RELEASES, so any second
                      release of the same genome at another tau compounds it. A tau grid
                      is an evaluation device, never a release procedure.
Valid vs:             nothing — there is no matched quantity to compare it against. It
                      is reported as a stated limitation beside the per-release bound,
                      with the per-release bound never presented as covering the set.
Status:               PRIMARY — it is a reporting obligation in Arm 1, and an omitted
                      row here would misrepresent what the theorem covers.
Aliases:              "the joint bound gap", "the 44-release composition". ⚠ NEVER call
                      it "the joint bound" — there is no joint bound; that is the point.
```

### private_allele_survival  [privacy] [evaluation]
```
Is:                   ARM 1 ONLY. The share of a target's COHORT-SINGLETON alleles —
                      alleles carried by that member and by no other cohort haplotype —
                      that SURVIVE into the released path. It is the sharpest available
                      measure of whether a release can hand back the one thing that
                      identifies a donor outright.
                      ✅ **OURS IS ZERO BY CONSTRUCTION under `loo_cohort`**, and the
                      construction is the argument: the release is a path through the
                      donor panel, `output_support` Y_D contains only cohort-supported
                      traversals, and at K = 86 the target's own two haplotypes are out
                      of the panel. An allele no remaining donor carries is therefore
                      not emittable at any tau. PanMixer's knapsack carries no such
                      property: it releases UNSELECTED blocks verbatim, so a singleton
                      allele in an unselected block survives intact.
                      ⚠ THE SAME CONSTRUCTION IS WHY IT IS ALSO A UTILITY CEILING: what
                      cannot be emitted cannot be retained. Report it beside
                      `target_fidelity`'s Arm-1 ceiling, never alone — zero survival is
                      both the privacy result and the fidelity cost.
Computed by:          PLANNED — not implemented. It needs the per-target singleton set
                      (derivable from support(v) and the target's own alleles on the LOO
                      panel) and a released path to score against. ⚠ OURS is zero by the
                      construction above and should still be ASSERTED on the real
                      release rather than assumed; PANMIXER'S value on the same 44
                      targets and the same site set is UNMEASURED.
CHANGES MEANING WITH: whether the panel is the LOO panel (K = 86) or the full cohort
                      (K = 88) — at K = 88 the target's own haplotypes are donors, the
                      quantity is not zero by construction, and it stops measuring what
                      it is for. Also with the site set and the build (decision F), and
                      with how a -1 is treated: an UNCALLED target allele is not a
                      singleton and must not be counted as one.
Valid vs:             PanMixer on the identical 44 targets, site set, build and panel
                      definition — this is an ARM 1 head-to-head row. ⚠ Not comparable
                      across arms: ARM 2's target is not in the cohort, so it has no
                      cohort-singleton alleles and the quantity is n/a there.
Status:               PRIMARY in ARM 1 — a mechanism-level claim our construction makes
                      and PanMixer's does not. n/a in ARM 2.
Aliases:              "singleton survival", "private-allele leakage". ⚠ Never "the
                      singleton rate", which is a property of the COHORT, not of a
                      release.
```

---

### panmixer_linkage_score  [attack] [evaluation] [external]

```
Is:                   THE UMBRELLA NAME for the three metrics PanMixer's own privacy
                      evaluators emit, which we now quote because their evaluators have
                      been run UNMODIFIED on both tools' releases. Their names, verbatim
                      from their JSON, because renaming an external tool's output is how
                      a comparison stops being one:
                        * `genotype_g_to_gstar` (`gap_score.py`) and
                          `haplotype_g_to_gstar` (`diploid_gap_score.py`) — the TRUE
                          subject's score MINUS the best candidate's.
                        * `hap_original_rank` out of `hap_original_rank_out_of`
                          (`MIA_privacy.py`) — the true subject's rank among the 44.
                        * `hap_bin_match_frac_vs_original` — how much of the original
                          survives in the release. HIGHER = less obfuscated.
                      ⚠⚠ DIRECTION, or every number reads backwards: a POSITIVE
                      `g_to_gstar`, or a rank of 1, means THE ATTACK IDENTIFIED THE
                      TARGET. Negative, or a rank far from 1, means it did not.
Computed by:          THEIR code, not ours: `external/PanMixer/tools/downstream/privacy/gap_score.py`,
                      `diploid_gap_score.py`, `MIA_privacy.py`, driven by
                      `scripts/pm_eval_run.sh` over the staging built by
                      `scripts/pm_eval_setup.py`; tabled by
                      `scripts/pm_eval_results.py`. ⚠ THREE evaluators, not four:
                      `sampled_gap_score.py` reads a PRE-stacker object and
                      `external/PanMixer/starting_data/chr21/subjects/` is absent from
                      this checkout — see its `docs/data.md` row, state MISSING.
CHANGES MEANING WITH: the candidate set it ranks against (44 here), the site set, and
                      the build — their sweep ran on vg-1.36 data and ours on vg-1.68,
                      so this is NOT a same-build comparison (decision F). Also with
                      WHICH object is scored: their released `obfuscated.vcf.gz`, not
                      the pre-stacker `new_haplotypes.npy` in their sweep output.
Valid vs:             each tool against itself across its own parameter, and across
                      tools only with the cross-build caveat attached.
Status:               DESCRIPTIVE for ARM 1 — a labelled external column, which is what
                      constraint 6 permits without a new pre-registration. ⚠⚠ It is NOT
                      `attacker_accuracy` and must never be reported as one: that term
                      names OUR pre-registered attackers, which remain unrun on Arm 1.
                      ⚠⚠ **THE SAME SCORE IS ALSO THE STATISTIC OF PRE-REGISTERED ATTACK
                      A3**, where it is re-posed as a 1-of-2 discrimination and IS a
                      claimable `attacker_accuracy` — the two uses must not be confused.
                      As PanMixer runs it, the score ranks 44 candidates and the result
                      is their `g_to_gstar`; as A3 runs it, it scores exactly two
                      candidates and the result is a balanced accuracy bounded by
                      `p_succ_bound`. `scripts/arm2_attack_a3a4.py` implements A3 and is
                      held to THEIR code by `tests/test_a3_matches_panmixer.py`, which
                      pins golden values obtained by calling their own function.
                      ⚠⚠⚠ A3's pre-registered positive control CANNOT FAIL, which is a
                      property of this score and not of the harness: it credits only
                      agreement, never penalises disagreement, and its weights are
                      non-negative, so `score(y, q)` over the agreeing sites is a SUBSET
                      SUM of `score(y, y)` and handing the attack the target's own path
                      forces a win. The control tests the plumbing; the power statement
                      is the leak-fraction curve in that script.
                      ⚠⚠⚠ AND A NULL HERE PROVES LESS THAN IT APPEARS TO. Our `tau` = 0
                      release defeats all three, and `tau` = 0 applies NO protection, so
                      the nulls are a property of the RELEASE FORM — a complete mosaic of
                      other people`s haplotypes — not of the tilt. Constraint 7: a
                      failing attacker is a LOWER bound on distinguishability. The upper
                      bound is `acc_ceiling`, which is proof-derived.
Aliases:              "their attacks", "the gap scores", "the MIA rank". ⚠ Never "the
                      attacker accuracy" and never "the privacy result" unqualified.
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
