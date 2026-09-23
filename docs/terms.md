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
⚠ UNSETTLED: what counts as a position is open question 8 in docs/plan.md. Do not
   hard-code a segmentation until that is decided.
Is:                   T, the number of ordered positions in the hidden Markov chain —
                      the steps the sampler takes along the chromosome, at each of
                      which it picks a donor haplotype and emits that donor's allele.
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
CHANGES MEANING WITH: whether the scoring panel includes the target (as shipped it
                      DOES, obfuscate.py:130,132, which pins eps_j near log(2N)); and
                      whether the block took the HMM or the allele-frequency branch —
                      the thresholds for scoring and sampling differ (>=1 vs <=1 anchor).
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
Status:               SECONDARY — a COHORT-level metric, near-blind to target fidelity.
                      Never quote it as evidence a release preserved the TARGET.
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
Status:               SECONDARY — cohort-level, same blindness as af_loss.
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
Is:                   THE integration artifact: the released allele vector for one
                      target, shape (n_sites, 2) int16, -1 = missing, row-aligned to
                      that chromosome's pangenome VCF record order. Every PanMixer
                      attack and utility evaluator consumes exactly this file.
Computed by:          PLANNED for us; PanMixer writes it from tools/panmixer/stacker.py
CHANGES MEANING WITH: the site axis. It is positional — a row-order mismatch silently
                      scores the wrong variants rather than erroring.
Valid vs:             another new_haplotypes over the identical site axis and cohort.
Status:               PRIMARY — our sampler MUST emit this, or the benchmark is not shared.
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
Status:               PRIMARY
Aliases:              "LOO". ⚠ Even a full loo_cohort is an APPROXIMATION of our
                      threat model: the PGGB graph TOPOLOGY was built from all HPRC
                      assemblies, so G itself saw the target. Never claim full
                      externality on PanMixer's data.
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
Valid vs:             the other arm at matched empirical privacy.
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
