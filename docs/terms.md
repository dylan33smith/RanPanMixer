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

### phi_v  [utility]
```
Is:                   The local per-position similarity in [0,1] between the target's
                      allele and a candidate donor's allele at chain position v.
                      ✅ **UNDER MODEL C (decided 2026-09-28) THIS COLLAPSES TO AN
                      INDICATOR**: phi_v = 1 if the drawn donor carries the target's
                      allele at v, else 0. There is no block to aggregate over, no
                      per-block normaliser, and no aggregation choice to make.
Computed by:          PLANNED — src/ranpanmixer/utility.py
CHANGES MEANING WITH: the missing_policy (an excluded donor has no phi_v at all), and
                      whether the target's own allele is missing at v — in which case
                      the position contributes no tilt and every donor gets factor 1.
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
Computed by:          PLANNED. ⚠ support(v) needs NO computation — measured 2026-09-28,
                      the pangenome VCF's own INFO/AN field is bit-identical to
                      np.sum(pangenome != -1, axis=(0,2)) across all 340,824 chr21
                      records.
CHANGES MEANING WITH: whether they are uniform (1/T) or support-weighted, and whether
                      length-weighting is ever adopted. The normalization is what keeps
                      u_path inside [0,1]; weights that do not sum to 1 break the
                      calibration, not just the emphasis.
Valid vs:             the same weighting scheme.
Status:               SECONDARY
Aliases:              ⚠ RENAMED from `beta_t` 2026-09-29, same reason as `phi_v`.
⚠ THE BUDGET IS INVARIANT IN T. Since u = sum_v beta_v phi_v with sum beta = 1, u lands
                      in [0,1] however the path is partitioned, and the total tilt
                      between a perfectly-matching and a non-matching path is exp(eta_tau)
                      regardless of T. Partitioning finer does NOT spread the budget
                      thinner: each position moves less and there are proportionally
                      more of them. Measured at tau = 0.5 (eta = 0.5493): per-position
                      tilt 1.0000055 at T = 100,757 (blocks) against 1.0000016 at
                      T = 340,824 (variants), total exp(eta) = 1.7321 either way.
                      Model C is in fact BETTER on utility at matched tau, because a
                      block design must take a whole block from one donor while C can
                      collect partial credit position by position.
⚠ THE REAL CONSEQUENCE OF FINE T IS NUMERICAL. The per-position tilt under C is only
                      **13.5x** the float32 resolution of 1.19e-07, against 45.7x under
                      a block partition. float64 for the stored alphas is load-bearing,
                      not a preference — float32 would consume a meaningful fraction of
                      the target's entire influence.
⚠ THE THIN-BUDGET PROBLEM IS A PROPERTY OF tau, NOT OF T. At tau = 0.5 the best-matching
                      path is only 1.73x more likely than the worst; that is what
                      tau = 0.5 MEANS. exp(arctanh(tau)) = 1.1055 / 1.2910 / 1.7321 /
                      4.3589 at tau = 0.1 / 0.25 / 0.5 / 0.9 (derived). If the achieved
                      utility is too low, the lever is tau — a weaker privacy claim —
                      not a different segmentation.
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
✅ SETTLED 2026-09-28 (was open question 8): **one chain position per VARIANT** —
   model C. Model B (one position per LD block) is the retained fallback.
Is:                   T, the number of ordered positions in the hidden Markov chain —
                      the steps the sampler takes along the chromosome, at each of
                      which it picks a donor haplotype and emits that donor's allele.
Computed by:          PLANNED — one position per VCF record inside `chain_span`.
                      chr21: 305,887 positions (340,824 records less the 34,345 below
                      the coordinate cut).
CHANGES MEANING WITH: the model. C = 340,824 raw / 305,887 after the cut; B = 100,757
                      LD-block entries. Runtime and numerical headroom both scale with
                      it — see beta_t for why finer T costs float32 headroom but NOT
                      privacy budget.
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
                      artifact is indexed by. chr21: **340,824**.
Computed by:          PLANNED — the record count of the deconstructed per-chromosome VCF.
CHANGES MEANING WITH: the VCF. A re-sort, a re-filter or a different `vg` version changes
                      which row is which SILENTLY, because every join in the pipeline is
                      positional. ⚠ Distinct from T: `n_sites` is the axis, T is the
                      number of CHAIN positions, which is `n_sites` restricted to
                      `chain_span` (chr21: 305,887 of 340,824).
Valid vs:             the same VCF, verified by the `sites.tsv` digest.
Status:               SECONDARY
Aliases:              "the site axis", "rows". Never "T".
```

### allele_lengths  [implementation] [dataset]
```
Is:                   The length in base pairs of every declared allele at every record —
                      ragged, `num_alleles[i]` entries for row i. chr21: 767,324 declared
                      alleles over 340,824 records, max 90 at one site.
Computed by:          PLANNED — emitted by our VCF-to-matrix step, which is the only
                      moment REF/ALT are in hand.
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
                      over. ✅ chr21 (decided 2026-09-28, closed at both ends 2026-09-29):
                      **12,968,320 < POS <= 46,680,243** — the interval over which the
                      genetic map carries real markers — keeping **305,887** of 340,824
                      variants (89.75%). Everything outside is excluded from the chain
                      and from the release.
                      ⚠ **PROVISIONAL — Dylan to settle with his PI.** This is the cut
                      that lets v1 run, not a defended answer. Revisit it before any
                      result is published; the full option set is in the REVISIT block.
Computed by:          PLANNED — a fixed coordinate rule in the preprocessing config.
CHANGES MEANING WITH: the chromosome and the genetic map. It is a property of the
                      ASSEMBLY and the MAP, never of the target — which is what makes it
                      safe: it is fixed before the target is opened, so R_D stays a
                      function of D alone and Theorem 1 is untouched. ⚠ NEVER let an
                      exclusion become target-dependent (e.g. "drop positions where the
                      target is missing"); that would leak.
Valid vs:             the same span. Any utility or fidelity number must name it, since
                      the denominator changes with it.
Status:               PRIMARY — a release that does not state its span is not quotable.
Aliases:              "the coordinate cut"
```

**WHAT THE p ARM IS, and why chr21's is excluded.** Every chromosome has a pinch point,
the **centromere**, dividing it into a short arm — **p**, from French *petit* — and a long
arm, **q** (simply the next letter). On most chromosomes both arms carry genes. But five
human chromosomes — **13, 14, 15, 21 and 22** — are **acrocentric**: the centromere sits
very near one end, so the p arm is small and consists almost entirely of repetitive DNA plus
long tandem arrays of **rDNA** genes (the ones encoding ribosomal RNA). Those arrays are
near-identical *between* the five acrocentric chromosomes, which is why assemblers cannot
place them uniquely and why GRCh38 fills them with N. chr21's p arm is the 0–~12 Mb region.

**The live defect the cut fixes.** `np.interp` clamps by default, so all 26,018 variants
below the map's start (10,326,676) receive the map's first value; and the map's own first
interval (10,326,676 → 12,968,320, spanning the centromere) carries **dcM = 0.000000**.
Together that is **34,345 consecutive positions with identical cM**, so `Delta_x = 0`,
`P(switch) = 0`, and the sampler emits **one donor's real haplotype verbatim across
7.25 Mb** — departure (b) at roughly 200x the scale of the worst single block. This is a
privacy failure, not a utility one, and it is silent.

**Five independent reasons, each sufficient** (all measured or cited 2026-09-28):

| leg | evidence |
|---|---|
| the sequence is not real | GRCh38 chr21 5,010,000–10,814,560 is placeholder model sequence: 25 contigs separated by fabricated N-gaps, 22 of them exactly 50 kb. Base-pair distance there is fictitious |
| unmappable | 1000G's own strict accessibility mask rules **26,014 of 26,018** (99.98%) of those variants inaccessible to short reads |
| our own data is thin | cohort matrix is **36.76%** missing below the cut against **6.57%** above |
| no measured recombination | deCODE's pedigree map (4.5M observed crossovers) assigns **0.0 cM/Mb to every 1 Mb bin from 0–13 Mb** |
| none observed directly | a four-generation CEPH pedigree: *"Not a single allelic recombination was observed on the p-arm"* — 107 transmissions across ~38 Mbp of acrocentric short arms |

⚠ **REVISIT — this is a pragmatic cut to get v1 running, not a settled answer.** Recorded
here so it can be picked up once the base mechanism works:
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
306,479 = **0.19%**. The invariant *every chain position carries a genuine interpolated cM*
is now literally true, which is what makes the `assert no NaN` gate meaningful.

⚠ **ONE ZERO-DISTANCE CASE REMAINS, and it is correct to leave.** **618 interior map
intervals have dcM = 0**, covering 2.71 Mb, with ~1,650 kept variants inside them. These are
REAL map plateaus, and `P(switch) = 0` between two variants with no genetic separation is
the CORRECT model, not a defect.

✅ **WHAT CLOSING BOTH ENDS SIMPLIFIES.** Two planned artifacts fall away. `genetic_pos_valid`
is unnecessary — inside the span every position has a real cM by construction. And explicit
chain SEGMENTS are unnecessary: a long variant-free stretch inside the span produces a large
`Delta_x`, hence a switch probability approaching uniform, which is already the correct
behaviour. The chain does not need to be broken by hand anywhere inside `chain_span`.

⚠ **INTERPOLATE WITH `left=nan, right=nan`**, then assert no NaN survives the cut. That one
change converts a silent 34,345-position failure into a crash.

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
                      ⚠ GUARD REQUIRED: measured 2026-09-28 on chr21, 1 site of 340,824
                      has all 88 donors missing and 89 sites leave <= 1 donor. The
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

### reason_array  [implementation] [dataset]
```
Is:                   An auxiliary array row- and slot-aligned to the allele matrix,
                      holding a small integer code per entry saying WHY that entry is
                      `-1`. Shape `(88, n_sites)` int8, code 0 = called.
Computed by:          PLANNED — our VCF-to-matrix step, which is the only moment the
                      INFO column is in hand.
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
Valid vs:             itself across capacities on one cohort. It is a COST, unbounded
                      above — it is NOT a similarity and must not be used as phi_t
                      without the normalization recorded under phi_t.
Status:               DIAGNOSTIC
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
                      stratum is 0/0. Also with the row count: comparing a 44-row
                      cohort AF to a 45-row one shifts every frequency by ~1/45.
Valid vs:             the other arm at the identical site set and row count.
Status:               PANMIXER-ONLY. ⚠ **This is not a metric this project computes.**
                      Its signal is the change in the COHORT's allele frequencies caused
                      by overwriting one member's genotype column. We never modify the
                      cohort, so the quantity is not defined for our release object —
                      running it on our output would measure the distortion of pretending
                      our external target were cohort member k. Our utility axes are
                      `target_fidelity` (always against its ceiling) and `utility_retained`.
Aliases:              "U_AF"
```

### ld_loss  [utility] [evaluation]
```
Is:                   PanMixer's linkage-disequilibrium utility loss U_LD: the change
                      in pairwise LD among nearby sites induced by a release.
Computed by:          external/PanMixer/tools/downstream/utility_in/ld_loss.py:16-123
CHANGES MEANING WITH: the window — the shipped KB_WINDOW_SIZE = 5 is compared against
                      raw bp, so the window is 5 BASE PAIRS, not the paper's 5 kb; and
                      the site set (pangenome_mask.npy again), which under the build
                      mismatch makes it return (0.0, 0) rather than raising.
Valid vs:             the other arm at the identical window and site set. Run it both
                      as-shipped (for parity with the paper) and corrected, reporting both.
Status:               PANMIXER-ONLY — same as `af_loss`: a cohort-distortion metric that
                      is not defined for a release which never edits the cohort. Not a
                      metric this project computes. The honest analogue, if one is wanted,
                      is the LD DECAY CURVE of our release ensemble against the cohort's
                      over R draws at real distances — a different quantity, not `ld_loss`.
Aliases:              "U_LD"
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
Computed by:          PanMixer writes it from tools/panmixer/stacker.py. OPTIONAL for us.
CHANGES MEANING WITH: the site axis. It is positional — a row-order mismatch silently
                      scores the wrong variants rather than erroring. ⚠ The shape
                      assertion in PanMixer's gap-score evaluators is VACUOUS (both sides
                      equal site_mask.sum() by construction), so nothing upstream catches
                      this. Ship a `sites.tsv` beside any emitted copy and verify it.
Valid vs:             another new_haplotypes over the identical site axis and cohort.
Status:               COMPARISON-ARM ONLY. Emitting it is the cheapest way to run
                      PanMixer's ATTACK evaluators on our release unmodified, and it is
                      worth doing for the head-to-head. It is NOT a requirement of the
                      mechanism, and emitting it does not make PanMixer's cohort-level
                      utility metrics applicable to us — `af_loss` and `ld_loss` would
                      run and return meaningless numbers.
Aliases:              "the release", "the obfuscated haplotypes"
```

### loo_cohort  [dataset] [method]
```
Is:                   The leave-one-out cohort D = HPRC \ {target}, with every derived
                      artifact recomputed without the target: donor panel, allele
                      frequencies, support counts, the site list, and the block assignment.
Computed by:          PLANNED
CHANGES MEANING WITH: which artifacts were actually recomputed. Recomputing only the
                      donor panel (what PanMixer does) is NOT a loo_cohort.
Valid vs:             the same held-out individual in the other arm.
Status:               COMPARISON-ARM ONLY. ⚠ **Not a concept our mechanism needs.** Our
                      target is external to `G`, `D` and every panel, so there is nothing
                      to hold out: `K = 88`, unconditionally. `loo_cohort` exists solely
                      for the head-to-head arm, where a cohort member stands in as target
                      and externality must therefore be SIMULATED. See `target_fidelity`,
                      which reasons correctly about a genuinely external target throughout
                      and is the model for how these entries should read.
Aliases:              "LOO". ⚠ Even a full loo_cohort is an APPROXIMATION of the
                      external-target setting: the PGGB graph TOPOLOGY was built from all
                      HPRC assemblies, so G itself saw the stand-in. That is a limit of
                      the comparison arm, not of our design — with a real external target
                      the graph never saw them.
```

### target_fidelity  [utility] [evaluation]
```
Is:                   How much of the TRUE target the release retains — as opposed to
                      cohort-level utility, which a pure prior draw satisfies perfectly
                      while carrying zero target information.
Computed by:          external/PanMixer/tools/downstream/privacy/MIA_privacy.py:80-83,140
                      already computes match fractions of a release against the target's
                      ORIGINAL haplotypes, mechanism-agnostically.
CHANGES MEANING WITH: the site set it is measured over, and whether a half-missing
                      genotype is treated as hom-ref (MIA_privacy.py:49-51 does).
⚠ CEILING:            In the EXTERNAL-target setting this quantity has a maximum BELOW 1
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
Valid vs:             the other arm at matched empirical privacy, and only against a
                      stated ceiling computed on the same graph.
Status:               PRIMARY — without this axis our tau = 0 point looks like a free lunch.
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

### two_target_indistinguishability  [privacy] [evaluation] [gate]
```
Is:                   The test our guarantee ACTUALLY makes, and the reason the
                      evaluation cannot simply be inherited. Run the mechanism on two
                      admissible targets X and Y under identical settings, draw R
                      releases of each, and estimate the distance between the two
                      release distributions. The theorem says it is <= tau for EVERY
                      such pair; PanMixer's design cannot bound it at all, because it
                      releases unselected blocks verbatim, so two inputs differing in
                      an unselected block have DISJOINT supports and TV = 1.
Computed by:          PLANNED — ours to write. ⚠ **NOTHING in PanMixer's evaluation
                      suite tests this** (established 2026-09-28 by reading every
                      evaluator). Their metrics test their claim.
CHANGES MEANING WITH: the estimator and the support size — see tv_empirical: plug-in TV
                      over a large support is biased UPWARD, so a value above tau at
                      scale is the estimator misbehaving, not a violated theorem. Only
                      exact enumeration on a toy support can falsify Theorem 1.
                      Also with the PAIR: the bound is over all pairs, so a single
                      favourable pair proves nothing. Report the worst pair found.
Valid vs:             tau, at matched settings, on a stated support.
Status:               PRIMARY — this is the headline privacy claim.
Aliases:              "the matched-pair experiment", cf. `A-ATK-matched-pair`
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
