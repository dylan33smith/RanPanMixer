## Section 4 — Our Mechanism

This section describes the mechanism this project is building. Its working title is **"Private Genome Path Release against a Public Pangenome Graph"**, and the durable source document for everything in Part A is the proposal PDF in this repository at `paper/Private_Genome_Path_Release.pdf`. Part A builds the theory from scratch. Part B compares it, step by step, with PanMixer (the published obfuscation tool dissected in Section 3; for readers arriving here first, PanMixer is an existing, published method that protects people who *contributed* to a pangenome graph by rewriting some of their haplotype blocks inside the released cohort file).

Two housekeeping notes before we start.

**Where the numbers come from.** Every number in this section is either (i) quoted from a published paper, (ii) measured in this project, or (iii) arithmetic done here from (i) or (ii). Measured-in-this-project numbers are attributed to durable locations in this repository: `docs/memory.md` (the permanent dated ledger of what was run and found), `docs/bugs.md` (the defect register), `docs/data.md` (the artifact registry), and the raw run logs under `logs/`. The PanMixer source we read and executed is the pinned checkout `external/PanMixer` at commit `c182c38d5bc8bb6f00f4b0b101207c4a009ca045`. The published PanMixer paper is `paper/s41467-026-77591-0_reference.pdf`; the superseded preprint is `archive_docs/Blindenbach2026_PanMixer.pdf`, and where the two disagree we use the published version.

**Symbol reset.** Section 1 introduced the document-wide symbol table, including the warning that PanMixer's paper and our proposal reuse the same Greek letters for unrelated things. That warning belongs early and is repeated here in full, because Section 4 is where both sets of symbols appear on the same page. Every symbol below is also re-glossed in words at the point it is used.

| Symbol | Read it as | Meaning in Section 4 |
|---|---|---|
| G | "the graph" | The fixed, public pangenome graph |
| D | "the cohort" | The public set of genomes the graph was built from, D = {g_1, …, g_n} |
| g | "the target" | One external genome, **not** in D, whose privacy we protect |
| p | "the input path" | The target's path through G — the private input |
| q | "another input path" | A second, hypothetical input path, used only inside proofs and attacker arguments |
| y | "an output path" | A candidate sanitized path we might release |
| P(G) | "the admissible inputs" | The set of paths that could legitimately be given to the mechanism |
| Y_D | "the output support" | The fixed set of paths the mechanism is allowed to emit |
| R_D | "the baseline" | A probability distribution over Y_D built only from public data |
| u(p, y) | "utility" | A number in [0, 1]: how well output y serves input p |
| tau | "the privacy parameter" | A number in [0, 1) chosen in advance; smaller = more private |
| eta_tau | "the tilt strength" | arctanh(tau); how hard we are allowed to bend the baseline |
| epsilon_tau | "the LDP budget" | 2·arctanh(tau); the log-likelihood-ratio budget |
| Q_p | "the release distribution" | The distribution of outputs when the input is p |
| Z_p | "the normalizer" | The sum that makes Q_p add up to 1 |
| T | "number of blocks" | How many ordered segments the path is cut into |
| K | "number of states" | How many cohort haplotypes the model can copy from |
| Z_t | "the state at block t" | Which cohort haplotype the model is copying at block t |
| h_{t,j} | "the local traversal" | The actual stretch of sequence that state j supplies at block t |
| beta_t, phi_t | "block weight, block similarity" | The two pieces of the per-block utility |
| psi_t(j) | "the local potential" | The per-block tilt factor, exp(eta_tau·beta_t·phi_t) |
| ‖ | "concatenate" | Glue sequences end to end, in order |

**Collisions to hold in your head.** These are genuine clashes of notation between documents, not different views of one quantity.

| Letter | Meaning A | Meaning B | Meaning C |
|---|---|---|---|
| eta | **Our** eta_tau = tilt strength (this section) | PanMixer's per-move **utility loss** eta_j (Section 3) | — |
| epsilon | **Our** epsilon_tau = a local-DP budget (this section) | PanMixer's per-block **privacy risk** eps_j, and its operating point eps_private = 0.002 (Section 3) | The HMM's per-anchor **mismatch probability**, 1e-4, in PanMixer's code (Section 3, defect (c)) |
| p, q | **Our** two input paths (this section) | **Allele frequencies** in the Hardy–Weinberg genotype model (Section 3) | **Stay and switch transition probabilities** in the Li–Stephens model (§4.13, §4.14 below) |
| r | The recombination-rate constant **r = 1.26** in PanMixer's published Fig. 6 | The squared LD correlation **r²** between two sites (Section 1) | — |
| n | Cohort size, n = 44 individuals (this section, §4.2) | Number of **donor states** in the stay/switch formula: 86 in PanMixer's code, K = 88 for us (§4.13) | Number of **anchors in a block**, in the defect (c) cap (§4.13) |

Wherever one of these appears below, the surrounding sentence says which meaning is intended.

---

## Part A — The Theory

### 4.1 What kind of guarantee do we actually want?

Section 3 described how PanMixer is evaluated: build a specific attacker (the "gap-score" linkage attack, and in the published version an additional membership-inference analysis), run it against the obfuscated output, and show that it fails. That is an **empirical** guarantee. It is genuinely informative — an attack that succeeds is a proof of insecurity — but it has a structural weakness: it tells you nothing about the attacker you did not think of. A gap-score attack failing does not imply that a k-mer-based attack fails, or that an attacker with a better reference panel fails, or that an attacker five years from now with a larger public cohort fails. This is not a hypothetical worry in genomics: the result usually credited to **Homer and colleagues (2008)** showed that even *aggregate* allele frequencies — a summary that everyone had assumed was safe — could reveal whether a known individual was in a pool. That is standard background in the field rather than anything measured in this project, and it is the canonical example of an assumed-safe release later turning out to be attackable.

We want the other kind of guarantee: a **worst-case, attacker-independent** one. The trick for getting such a thing — and it is the central trick of the whole privacy literature — is to stop reasoning about attackers and start reasoning about **distributions**.

Here is why that works. Suppose the release mechanism is random. Then for each possible private input there is a *distribution* over what might come out. Call the distribution induced by input p the release distribution, written Q_p. An attacker who sees one released output and tries to decide "did this come from input p, or from input q?" is, mathematically, running a **binary hypothesis test** between two distributions, Q_p and Q_q. (A hypothesis test is just a rule that maps each possible observation to one of the two answers; its accuracy is the probability that the answer is right.) And there is a classical fact that caps how well *any* such rule can possibly do.

That fact is stated in terms of **total variation distance** (TV). For two probability distributions P and Q over the same finite set of outcomes Y_D,

TV(P, Q) = (1/2) · sum over y in Y_D of |P(y) − Q(y)|.

*In words:* line the two distributions up outcome by outcome, take the absolute gap at each outcome, add all the gaps up, and halve. TV is 0 when the two distributions are identical and 1 when they are completely disjoint (no outcome has positive probability under both). The proposal defines TV in its more general form as the supremum — the least upper bound — over all subsets A of Y_D of |P(A) − Q(A)|; on a finite outcome set like ours the supremum is attained, so "supremum over events" and "maximum over events" mean the same thing here, and both equal the half-sum-of-gaps formula above.

**Micro-example.** Three possible outputs y1, y2, y3. P = (0.5, 0.3, 0.2) and Q = (0.4, 0.4, 0.2). The gaps are 0.1, 0.1, 0.0; their sum is 0.2; TV = 0.1. Now be the attacker. You see one output and must guess which distribution produced it, knowing each was equally likely a priori (a **prior** is the attacker's belief before seeing any data; "equal priors" means 50/50). The best rule is to name whichever distribution gives the observed outcome higher probability: on y1 guess P (0.5 > 0.4), on y2 guess Q (0.4 > 0.3), and on y3 it is a tie (0.2 = 0.2), so any rule will do. Ties genuinely do not matter: whatever you do on y3, you are right with probability 0.5 · 0.2 / (0.5 · 0.2 + 0.5 · 0.2) = 1/2 of the times y3 occurs, so y3 contributes the same amount however you break it. Assigning y3 to P for bookkeeping, your success probability is 0.5 · P({y1, y3}) + 0.5 · Q({y2}) = 0.5 · 0.7 + 0.5 · 0.4 = **0.55**. And 0.55 = (1 + 0.1)/2 = (1 + TV)/2. This is not a coincidence; it is the general identity.

**The identity we will lean on:** under equal priors, the best achievable accuracy of *any* binary test distinguishing Q_p from Q_q is exactly (1 + TV(Q_p, Q_q))/2.

So if we can *prove* TV(Q_p, Q_q) ≤ tau for every pair of possible private inputs p and q, we have proved that no attacker of any kind, with unlimited compute and full knowledge of our algorithm, can beat accuracy (1 + tau)/2 at telling the target's release apart from anyone else's. That is the guarantee we are after. It is the content of Definition 1 in the proposal, which states it in two strengths: the weaker form asks only that the target be hard to distinguish from each of the n cohort members (max over i of TV(Q_p, Q_{p_i}) ≤ tau), and the stronger form asks it for *every* pair of admissible inputs. Our mechanism achieves the stronger form.

### 4.2 The setting: what is public and what is private

Be precise about the threat model, because the guarantee is only as good as its assumptions.

- **G**, the pangenome graph, is **public**. It was built from the cohort **D = {g_1, …, g_n}**, which is also public — these are the HPRC assemblies, already released. In the data this project actually works with, the released cohort VCF carries 45 samples; one of them is **CHM13**, which is a *human* cell line derived from a complete hydatidiform mole (an abnormal human conception in which all chromosomes come from the father, so the resulting line is homozygous nearly everywhere — an excellent assembly substrate, but a reference assembly rather than a consenting cohort donor whose privacy is at stake). The pipeline drops CHM13 and the X chromosome, leaving **44 individuals = 88 haplotypes**, which is also the cohort size stated in the published PanMixer paper (measured in this project and cross-checked against the paper; see `docs/memory.md`, entries dated 2026-09-18).
- Each cohort genome g_i has a path p_i = Map(g_i, G) through the graph. These paths are **public** too: they are literally in the graph.
- **g** is an external target genome that is **not** in D. This is the crucial structural difference from PanMixer and we will return to it repeatedly. The target did not contribute to the graph; the graph knows nothing about them.
- **p = Map(g, G)** is the target's mapped path. **This, and only this, is private.**
- The attacker knows G, D, all the cohort paths, our entire algorithm, our utility function, and the value of tau. The only thing they do not know is the outcome of our random number generator. (This is Kerckhoffs's principle from cryptography, imported wholesale: assume the enemy knows the system, keep only the key — here, the randomness — secret.)

Two more objects:

- **P(G)**, the *admissible input-path space*. The proposal introduces this set without pinning it down, so state a concrete reading and flag it: for our purposes an admissible path is **any valid traversal of G** — a sequence of nodes in which every consecutive pair is joined by an edge — regardless of whether any cohort member happens to carry it. This matters because Theorem 1 is quantified over all p, q in P(G); the larger P(G) is, the stronger the theorem, and the reader cannot judge its strength without knowing what is in there. **Open question:** whether the implementation should restrict P(G) further (for example to paths a real mapper could actually produce), and whether doing so would weaken the guarantee's interpretation.
- **Y_D**, the *output support*: the fixed set of paths the mechanism is permitted to emit, a subset of P(G). The proposal allows several choices — the complete cohort haplotype paths, synthetic paths from a cohort haplotype model, mosaics of cohort-supported local traversals, or graph-valid paths generated by an HMM. ("Mosaic" here means a path stitched together from segments of several different cohort haplotypes, the way a real chromosome is a mosaic of its ancestors' chromosomes because of recombination.) The requirement is that **Y_D is chosen without looking at the target.**

### 4.3 The central idea: do not edit, resample

PanMixer's release rule is *editing*: take the target's own haplotype and overwrite some blocks. Our release rule is fundamentally different:

> We never edit the target's path. We **draw a brand-new path** from a probability distribution over cohort-supported paths, and we release that.

The target's path influences *which* path gets drawn — otherwise the output would be useless — but only by **re-weighting** the probabilities, and only by a **bounded** amount. Nothing of the target is ever copied through unchanged by construction; every symbol in the output is a symbol that came from a cohort-supported candidate.

This is what makes a worst-case bound reachable. If the output is always drawn from the same fixed menu, then two different inputs produce two distributions over *the same menu*, and we can hope to bound the distance between them. If instead the output contains verbatim copies of the input, then two inputs that differ produce outputs that differ with certainty, and no bound is possible. (Part B, §4.12, makes this argument precise against PanMixer.)

### 4.4 Ingredient 1: the target-independent baseline, R_D

**R_D** is a probability distribution over the output set Y_D, built **only** from the public graph G and the public cohort D. Think of it as a generative model of "what a plausible human chromosome looks like, according to the cohort" — the same kind of object as the Li–Stephens copying HMM that PanMixer uses (Section 3), but put to a different purpose. Two requirements from the proposal:

1. **R_D(y) > 0 for every y in Y_D.** No output the mechanism could emit is assigned probability zero. (If some output had probability zero under the baseline but positive probability under a tilt, likelihood ratios would blow up.)
2. **R_D does not depend on the target in any way.**

**Why requirement 2 matters — a worked counterexample.** Imagine we cheated and let R_D depend on the target: say we built the baseline HMM from D *plus the target*, or we pruned Y_D to "paths within 5% of the target." Now consider two possible targets, Alice and Bob, whose paths differ at exactly one rare variant. Alice's pruned menu contains paths carrying that rare allele; Bob's does not. An attacker who sees an output carrying the rare allele knows instantly it was Alice. The supports are disjoint on that slice, TV = 1, and the guarantee is gone — not weakened, *gone*. This is exactly the failure mode the proposal's "common support" condition rules out.

Note a pleasant consequence of our threat model: because the target is **external** to D, every quantity derived from the graph (allele frequencies, variant support counts, cohort haplotypes, LD blocks) is automatically target-independent. We get condition 2 nearly for free. PanMixer cannot get it for free, because its target is *inside* the cohort — which is precisely the origin of released-code defect (c), where the privacy scorer scores the target against a panel that contains the target. We return to the measured size of that effect in §4.13.

### 4.5 Ingredient 2: a bounded utility function, u(p, y)

**u(p, y)** is a number between 0 and 1 saying how well output y serves input p. Larger is better. u = 1 means "perfect for this target"; u = 0 means "useless for this target."

The proposal's blockwise construction: cut the path into **T** ordered blocks (LD blocks, graph blocks, or anchor sites — the block segmentation described in Section 3 works fine here), then

u(p, y) = sum over t = 1..T of beta_t · phi_t(p_t, y_t)

where:

- **p_t** is the input path restricted to block t; **y_t** is the output path restricted to block t.
- **phi_t(p_t, y_t)** is a per-block similarity in [0, 1]. 1 = the output block is perfect for the target at this block, 0 = worthless.
- **beta_t ≥ 0** is a per-block importance weight, and the weights sum to 1 across all T blocks: sum of beta_t = 1.

*In words:* score each block separately on a 0-to-1 scale, then take a weighted average of those scores, with the weights saying how much each block matters.

Because each phi_t is in [0, 1] and the beta_t are non-negative and sum to 1, u is a **weighted average of numbers in [0, 1]**, so u itself is automatically in [0, 1]. That is the whole point of the construction — boundedness is structural, not something we have to check case by case.

**A concrete phi_t.** The proposal suggests a population-weighted k-mer agreement. A **k-mer** is a length-k substring of DNA; K(a) denotes the set of k-mers appearing in sequence a. The notation 1[x in K(a)] is an **indicator**: it is 1 if k-mer x occurs in a, and 0 otherwise. Then

phi_t(a, b) = [sum over x of w(x)·min(1[x in K(a)], 1[x in K(b)])] / [sum over x of w(x)·max(1[x in K(a)], 1[x in K(b)])]

*In words:* the min is 1 only when the k-mer is in *both* sequences, and the max is 1 when it is in *either*, so the numerator is the weighted size of the intersection and the denominator the weighted size of the union. This is a **weighted Jaccard index**: weighted-intersection over weighted-union, which is 1 for identical k-mer sets and 0 for disjoint ones. The weight w(x) downweights k-mers that appear in everybody (uninformative) and upweights population-informative ones.

**Micro-example.** K(a) = {x1, x2, x3}, K(b) = {x2, x3, x4}, with weights w(x1) = 3, w(x2) = 1, w(x3) = 1, w(x4) = 5. The intersection is {x2, x3} with total weight 2. The union is {x1, x2, x3, x4} with total weight 10. So phi = 2/10 = 0.2.

**A concrete beta_t.** Section 3 described PanMixer's per-variant utility weight: each variant v carries weight 1/support(v), where support(v) is the number of cohort haplotypes traversing v — so rare variants, traversed by few haplotypes, are worth more. §4.13 below turns that raw weight into our beta_t and phi_t.

### 4.6 The exponential tilt

Now combine them. Fix tau in [0, 1) and define the **tilt strength**

**eta_tau = arctanh(tau) = (1/2)·log[(1 + tau)/(1 − tau)]**.

(arctanh is the inverse hyperbolic tangent: the function that undoes tanh. For now just treat it as a number that grows smoothly from 0 to infinity as tau grows from 0 to 1. §4.7 explains why *this* function and no other.)

The release distribution for input path p is

**Q_p(y) = R_D(y) · exp(eta_tau · u(p, y)) / Z_p**, for y in Y_D,

where **Z_p = sum over z in Y_D of R_D(z)·exp(eta_tau · u(p, z))** is the normalizer that makes the probabilities sum to 1.

**In plain English:** start with the public baseline probability of each candidate output. Multiply it by a bonus factor that is larger for candidates that are more useful to the target. Renormalize so the numbers are probabilities again. That is all.

The bonus factor exp(eta_tau · u) ranges over exactly [1, exp(eta_tau)], because u ranges over [0, 1] and exp is increasing. So the *worst* candidate keeps its baseline weight untouched and the *best* candidate is multiplied by at most exp(eta_tau). The tilt cannot change any candidate's weight relative to any other by more than a factor of exp(eta_tau). **That single sentence is the entire privacy argument**; everything else is bookkeeping.

The two extremes are worth internalizing:

- **tau = 0** gives eta_tau = 0, so the bonus factor is exp(0) = 1 for every candidate, so Q_p = R_D for *every* input p. The output is a pure public cohort sample that does not depend on the target at all. Perfect privacy, zero target-specific utility.
- **tau → 1** gives eta_tau → infinity, so the mechanism concentrates arbitrarily hard on the highest-utility candidate. Maximum utility, no privacy guarantee.

tau is the dial between those two, and it is chosen **before** looking at any data.

### 4.7 Why arctanh? Walking through the proof

The proof in the proposal is four lines. Here it is slowly.

**Step 1 — the bonus factor is bounded.** Since 0 ≤ u(p, y) ≤ 1, we have 1 ≤ exp(eta_tau · u(p, y)) ≤ exp(eta_tau). *(Exponentiating a number between 0 and eta_tau gives a number between 1 and e^eta_tau.)*

**Step 2 — the normalizer is bounded.** Z_p is a weighted average of those bonus factors, weighted by R_D (whose weights sum to 1 because it is a probability distribution). An average of numbers all lying in [1, exp(eta_tau)] must itself lie in [1, exp(eta_tau)]. So **1 ≤ Z_p ≤ exp(eta_tau)** for *every* input p. Consequently, for any two inputs p and q, the ratio Z_q/Z_p is at most (biggest possible Z_q)/(smallest possible Z_p) = exp(eta_tau)/1, and at least 1/exp(eta_tau): exp(−eta_tau) ≤ Z_q/Z_p ≤ exp(eta_tau).

**Step 3 — write down the likelihood ratio.** For any output y,

Q_p(y) / Q_q(y) = exp(eta_tau · [u(p, y) − u(q, y)]) · (Z_q / Z_p).

The R_D(y) factor is identical in numerator and denominator and cancels — *this is why the baseline must be target-independent.* If R_D depended on the target, the two R_D factors would be different functions and would not cancel, and there would be no bound at all.

**Step 4 — bound it.** Since both utilities lie in [0, 1], their difference lies in [−1, 1], so the first factor lies in [exp(−eta_tau), exp(eta_tau)]. Multiplying by the bound on Z_q/Z_p from Step 2:

**exp(−2·eta_tau) ≤ Q_p(y)/Q_q(y) ≤ exp(2·eta_tau)**, for every output y.

**Step 5 — convert a likelihood-ratio bound into a TV bound.** This is the standard step and it is where arctanh earns its keep. If two distributions have pointwise likelihood ratio everywhere within [exp(−s), exp(s)], then TV ≤ (exp(s) − 1)/(exp(s) + 1), and that expression is exactly tanh(s/2). Here s = 2·eta_tau, so TV ≤ tanh(eta_tau). And we *chose* eta_tau = arctanh(tau), so tanh(eta_tau) = tanh(arctanh(tau)) = **tau**.

That is the answer to "why arctanh": the likelihood-ratio-to-TV conversion is a tanh, so to land on exactly tau we invert the tanh. The calibration is not a heuristic; it is the exact inverse of the conversion.

**Theorem 1 (from the proposal).** For any two admissible input paths p and q in P(G): TV(Q_p, Q_q) ≤ tau, and moreover exp(−epsilon_tau) ≤ Q_p(y)/Q_q(y) ≤ exp(epsilon_tau) for every y, with **epsilon_tau = 2·arctanh(tau) = log[(1 + tau)/(1 − tau)]**.

**Corollary (attacker success).** Combining with the identity from §4.1: under equal priors, no binary attacker can exceed accuracy **(1 + tau)/2** at distinguishing the target's release from any other input's release.

**Corollary 1 of the proposal (cohort mixture).** Let pi = (pi_1, …, pi_n) be any public prior over the n cohort members (non-negative numbers summing to 1) and let Q_D be the mixture distribution "pick a cohort member at random according to pi, then run the mechanism on their path": Q_D = sum_i pi_i · Q_{p_i}. Then TV(Q_p, Q_D) ≤ tau as well. The proof uses the **convexity of total variation distance**, which is the inequality

TV(P, sum_i pi_i·Q_i) ≤ sum_i pi_i·TV(P, Q_i).

*In words:* blending distributions together cannot push the blend further from P than the worst ingredient is. Since each TV(Q_p, Q_{p_i}) ≤ tau, the weighted average of those is also ≤ tau. **Interpretation:** the target's release is hard to distinguish not only from any *specific* cohort member's release but also from a release generated by a randomly chosen cohort member — i.e., the target blends into the cohort as a whole.

**The dial, numerically** (arithmetic computed here from the formulas above):

| tau | eta_tau = arctanh(tau) | epsilon_tau = 2·arctanh(tau) | exp(epsilon_tau) = (1+tau)/(1−tau) | Max attacker accuracy (1+tau)/2 |
|---|---|---|---|---|
| 0.01 | 0.0100 | 0.0200 | 1.0202 | 0.505 |
| 0.05 | 0.0500 | 0.1001 | 1.1053 | 0.525 |
| 0.10 | 0.1003 | 0.2007 | 1.2222 | 0.550 |
| 0.25 | 0.2554 | 0.5108 | 1.6667 | 0.625 |
| 0.50 | 0.5493 | 1.0986 | 3.0000 | 0.750 |
| 0.90 | 1.4722 | 2.9444 | 19.000 | 0.950 |

**Relation to local differential privacy.** *Differential privacy* (DP) is the standard framework in which a mechanism's output distribution is required to change by at most a factor exp(epsilon) when the input changes. *Local* DP (LDP) is the variant where the randomization is applied by each individual to their own data before it ever leaves their hands — which is exactly our situation: one target, one path, one randomized release, no trusted curator. The pointwise ratio bound in Theorem 1 *is* an epsilon_tau-LDP guarantee over the admissible input-path domain. It is strictly stronger than the TV statement (pointwise ratio bounds imply TV bounds, but not conversely). The proposal nevertheless leads with tau rather than epsilon_tau because tau has the direct operational reading "attacker accuracy is at most (1 + tau)/2," whereas an epsilon is notoriously hard for non-specialists to interpret. Remember the collision: this epsilon_tau has nothing to do with PanMixer's per-block privacy risk eps_j or its operating point eps_private = 0.002.

### 4.8 The conditions, and what breaks if each fails

The proposal's §2.5 lists **three** conditions that the theorem's proof actually depends on. There is a fourth requirement that the proposal states elsewhere as practical guidance (in its §2.6 on sampling versus Viterbi, and its §2.8 on repeated releases) rather than in the numbered list. We keep it as row 4 and label it as ours, because an implementer who ignores it loses the guarantee just as surely.

| # | Condition | What it means concretely | What breaks if violated |
|---|---|---|---|
| 1 | **Common support** (proposal condition 1) | Y_D fixed before looking at the target; no target-specific candidate lists | A candidate present for input p and absent for input q gives disjoint supports on that outcome; the likelihood ratio is infinite and TV jumps to 1 |
| 2 | **Target-independent baseline** (proposal condition 2) | R_D learned only from G and D | R_D(y) no longer cancels in Step 3 of the proof; there is no bound on the ratio at all |
| 3 | **Bounded utility** (proposal condition 3) | 0 ≤ u(p, y) ≤ 1 | Step 1 fails. If the raw utility has range Δu = u_max − u_min, the proposal's fix is to recalibrate: eta_tau = arctanh(tau)/Δu. If the range is *unbounded*, no calibration exists and there is no guarantee |
| 4 | **One release per genome, and sample rather than maximize** (our framing of the proposal's §2.6 and §2.8 guidance) | Do not independently resample the same genome for every query; never take the argmax | Two independent releases of the same genome give an attacker two draws. TV does not compose for free: the standard subadditivity bound for k independent releases is k·tau, which degrades fast and is only an upper bound, not a promise. Taking the argmax makes the mechanism deterministic, so TV becomes 0 or 1 and the guarantee is vacuous |

Condition 4 has a practical implication the proposal states explicitly: a deployed system should **sample one sanitized path and cache it**, serving that same path for every subsequent query about that genome.

Two further caveats from the proposal, which are limitations rather than conditions. First, the guarantee applies **only to the released sanitized path**; it says nothing about the raw genome data, which must be protected by other means. Second, the utility function must be chosen *before* inspecting target-specific attack outcomes — tuning phi_t against measured attack success on a real target would make the mechanism target-dependent through the back door. The proposal also notes that any target-dependent pruning of the HMM state space would require its own separate privacy analysis.

### 4.9 Implementation: the baseline as a hidden Markov model

We cannot enumerate Y_D — the number of cohort-supported paths through a chromosome is astronomically large. So we need R_D and the tilt to have a structure that permits exact sampling without enumeration. Hidden Markov models provide exactly that.

A **hidden Markov model (HMM)** has a sequence of unobserved ("hidden") states, each depending only on the previous one, and each emitting an observation. Here:

- The path is cut into **T** ordered blocks, t = 1, …, T.
- At block t the hidden state **Z_t** takes one of **K** values, indexing which cohort haplotype (or which cohort-supported local traversal) the model is copying at that block.
- **h_{t,j}** is the actual local sequence/traversal that state j supplies at block t.
- **rho(j)** is the probability of starting in state j; **A_t(i, j)** is the probability of moving from state i at block t−1 to state j at block t. Transitions with probability zero encode graph-incompatible concatenations — pairs of local traversals that cannot be glued together because the graph has no edge joining them.

The baseline is then R_D(z_1..z_T) = rho(z_1) · product over t = 2..T of A_t(z_{t−1}, z_t). A state sequence z_1..z_T induces a released path by concatenation:

**Phi(z_1..z_T) = h_{1,z_1} ‖ h_{2,z_2} ‖ … ‖ h_{T,z_T}**,

where **‖ means concatenation** — glue the block sequences end to end, in order, to get one long sequence.

This is precisely the Li–Stephens "copying" model that Section 3 described: the released chromosome is a **mosaic** that copies one cohort haplotype for a stretch, occasionally switches to another (a modelled recombination), copies that one for a stretch, and so on. The whole model is built from the cohort and does not mention the target.

Now the key algebraic fact. Because u is **additive across blocks**, the exponential tilt factorizes — a sum inside an exponential becomes a product of exponentials:

exp(eta_tau · u(p, z_1..z_T)) = exp(eta_tau · sum_t beta_t · phi_t(p_t, h_{t,z_t})) = product over t of **psi_t(z_t)**, where **psi_t(j) = exp(eta_tau · beta_t · phi_t(p_t, h_{t,j}))**.

So the tilted distribution is

Q_p(z_1..z_T) proportional to rho(z_1)·psi_1(z_1) · product over t = 2..T of A_t(z_{t−1}, z_t)·psi_t(z_t).

**This is still an HMM.** The transitions are unchanged public cohort transitions; the private target enters *only* through the per-block local potentials psi_t, which are bounded by construction: 1 ≤ psi_t(j) ≤ exp(eta_tau · beta_t). In HMM language, psi_t plays exactly the role an emission probability plays — the target's block is the "observation" and psi_t measures how compatible each state is with it. So we can sample from Q_p with standard HMM machinery.

### 4.10 Forward-filtering / backward-sampling (FFBS), from scratch

We need to draw a *whole state sequence* Z_1..Z_T from the tilted distribution, exactly — not approximately, and not the most likely one. FFBS does this in two passes.

**Pass 1: forward filtering (left to right).** Maintain, for each block t, a length-K vector alpha_t where alpha_t(j) is the probability that the chain is in state j at block t *given everything seen up to and including block t*.

- Initialize: r_1(j) = rho(j)·psi_1(j); c_1 = sum_j r_1(j); alpha_1(j) = r_1(j)/c_1.
- For t = 2..T: r_t(j) = psi_t(j) · sum_i alpha_{t−1}(i)·A_t(i, j); c_t = sum_j r_t(j); alpha_t(j) = r_t(j)/c_t.

*In words:* to get the unnormalized weight of state j at block t, first ask "how likely is it that the chain arrives in state j at all?" (that is the sum over the previous state i of where we were times the probability of moving i → j), then multiply by how well state j suits the target at this block. The division by c_t rescales the vector so its entries sum to 1.

The per-block renormalization by c_t prevents numerical underflow (multiplying T ≈ 100,000 small numbers together would otherwise round to zero in floating point). As a free by-product, log Z_p = sum over t of log c_t.

**Pass 2: backward sampling (right to left).** Now actually draw a sequence.

- Sample Z_T from the **categorical** distribution alpha_T. (A categorical draw just means "pick one of K options with the listed probabilities.")
- For t = T−1 down to 1, given the already-sampled Z_{t+1} = j, sample Z_t from weights proportional to alpha_t(i)·A_{t+1}(i, j), normalized over i.

**Why this is exact.** The joint distribution factorizes as P(Z_T) · P(Z_{T−1} | Z_T) · P(Z_{T−2} | Z_{T−1}) · … — each backward step is the exact conditional of Z_t given Z_{t+1} and all the evidence, because in a Markov chain, conditioning on Z_{t+1} makes Z_t independent of everything to the right. So the sequence produced is an exact draw from Q_p. The forward pass computes the conditionals; the backward pass consumes them.

**Worked example.** T = 3 blocks, K = 2 states. Start rho = (0.5, 0.5). Stay/switch transitions with stay probability a = 0.9 and switch probability b = 0.1 (with K = 2 these sum to 1). Take **tau = 0.5**, so eta_tau = arctanh(0.5) = 0.5493, and equal block weights beta_t = 1/3, so eta_tau·beta_t = 0.1831. Suppose the target matches state 1 in blocks 1 and 2 and state 2 in block 3:

| Block t | phi_t(state 1) | phi_t(state 2) | psi_t(1) | psi_t(2) |
|---|---|---|---|---|
| 1 | 1.0 | 0.0 | 1.2009 | 1.0000 |
| 2 | 0.8 | 0.2 | 1.1578 | 1.0373 |
| 3 | 0.0 | 1.0 | 1.0000 | 1.2009 |

Forward pass:

- t = 1: r_1 = (0.5·1.2009, 0.5·1.0) = (0.6005, 0.5000); c_1 = 1.1005; **alpha_1 = (0.5457, 0.4543)**.
- t = 2: prediction = 0.1·1 + 0.8·alpha_1 = (0.5365, 0.4635); r_2 = (1.1578·0.5365, 1.0373·0.4635) = (0.6212, 0.4807); c_2 = 1.1019; **alpha_2 = (0.5637, 0.4363)**.
- t = 3: prediction = 0.1 + 0.8·alpha_2 = (0.5510, 0.4490); r_3 = (1.0·0.5510, 1.2009·0.4490) = (0.5510, 0.5393); c_3 = 1.0902; **alpha_3 = (0.5054, 0.4946)**.

(The "prediction" line uses the stay/switch shortcut derived in §4.10.2: with two states, b·(alpha(1)+alpha(2)) + (a − b)·alpha(j) = 0.1·1 + 0.8·alpha(j).)

log Z_p = log(1.1005) + log(1.1019) + log(1.0902) = 0.2792, so Z_p = 1.322 — which sits inside the proof's interval [1, exp(eta_tau)] = [1, 1.7321], exactly as Step 2 of §4.7 promised.

Backward pass. Draw Z_3 from (0.5054, 0.4946); say we get **Z_3 = 2**. Then for t = 2 the weights are alpha_2(1)·A(1→2) = 0.5637·0.1 = 0.0564 and alpha_2(2)·A(2→2) = 0.4363·0.9 = 0.3926; normalized, **(0.126, 0.874)** — state 2 is strongly favoured, because switching is expensive. Say we draw **Z_2 = 2**. Then for t = 1 the weights are 0.5457·0.1 = 0.0546 and 0.4543·0.9 = 0.4089; normalized, **(0.118, 0.882)**.

Notice what the example teaches: at block 1 the tilt genuinely favours state 1 (alpha_1 puts 0.546 on it, above the 0.5 it would have with no tilt), but once the backward pass commits to state 2 at block 3, the stay probability drags the whole sequence toward state 2. The public transition model and the private tilt are genuinely competing, and at these settings the transition model usually wins. That is what "the target may only bend the baseline a little" looks like in practice.

### 4.10.1 Why sampling, not argmax

It is tempting to run **Viterbi** instead of sampling. Viterbi is the standard HMM dynamic-programming algorithm that returns the single most probable state sequence — an **argmax**, meaning the one highest-scoring answer, returned deterministically, with no randomness at all. **Using it here would destroy the guarantee.** A deterministic mechanism has a **point-mass** output distribution — probability 1 on one outcome and 0 on everything else. For two inputs whose argmaxes differ, the two point masses sit on different outcomes, so their supports are disjoint and TV(Q_p, Q_q) = 1: the worst possible value, and an attacker accuracy bound of (1 + 1)/2 = 1, i.e. no protection. Randomness is not a nuisance to be minimized here; it *is* the privacy. The proposal states this explicitly when it contrasts forward sampling with Viterbi decoding.

One more piece of good news: the released path is p̃ = Phi(Z_1..Z_T), a deterministic function of the sampled states. **Deterministic post-processing can never increase TV** — if you could raise the distance by relabelling outcomes, you could build a better test on the relabelled outcomes and beat the optimal test on the originals, a contradiction. So TV(Law(p̃ | p), Law(p̃ | q)) ≤ TV(Q_p, Q_q) ≤ tau. The guarantee survives the conversion from states to an actual path.

### 4.10.2 Cost

With a **dense** K-by-K transition matrix, the forward recursion costs O(K²) per block (each of K new states sums over K old states), hence O(T·K²) overall; backward sampling costs O(K) per block, hence O(T·K); computing the potentials costs O(T·K). The dense forward pass dominates.

But PanMixer's transition model is **stay/switch**: A_t(i, j) = a_t when i = j and b_t when i ≠ j — one probability for staying on the current donor haplotype and one for switching to any other. That structure collapses the sum. Writing S_{t−1} = sum_i alpha_{t−1}(i) for the total of the previous forward vector,

r_t(j) = psi_t(j)·[b_t·S_{t−1} + (a_t − b_t)·alpha_{t−1}(j)].

*Where that comes from:* every previous state contributes b_t to arriving at j, except state j itself which contributes a_t. So add b_t times everything, then correct the single diagonal term by (a_t − b_t). One shared scalar plus one per-state correction gives all K forward values in **O(K)** per block, so **O(T·K)** for the whole forward pass. The backward weights collapse identically: given a sampled next state j, the weight is a_{t+1}·alpha_t(j) for i = j and b_{t+1}·alpha_t(i) otherwise. No dense matrix is ever formed. **Total: O(T·K) time.** The proposal also notes a middle case: if only a sparse edge set E_t of state pairs can legally be concatenated at block t, the forward pass costs O(|E_t|) there, for a total of O(T·K + sum_t |E_t|).

Memory: storing every alpha_t is O(T·K).

**Our numbers on chr21.** On this project's rebuilt chr21 — the HPRC graph deconstructed to 340,824 VCF records, with LD blocks computed on the 1000 Genomes **30x GRCh38** panel rather than the panel the shipped pipeline wires in — the block dictionary has **T = 100,757 entries** and the state set is **K = 88 haplotypes** (44 individuals × 2 haplotypes). These are measured in this project and recorded in `docs/memory.md` under 2026-09-18, with the run logged in `logs/pm_chr21.log`. So T·K ≈ 8.87 million operations per haplotype — trivially cheap, a fraction of a second of arithmetic. Storing all forward vectors at 64-bit precision is 100,757 × 88 × 8 bytes ≈ 71 MB, which fits comfortably in memory. For comparison, our measured PanMixer run takes **68 s and 1.7 GB peak resident memory per subject-chromosome on chr21** (measured in this project; `docs/memory.md` 2026-09-18, raw figures in `logs/pm_obf.log`), so the sampler itself will not be the bottleneck. The proposal separately notes that for a diploid cohort of n = 500 individuals, K would be at most about 2n = 1,000 states before de-duplication — still cheap.

### 4.10.3 Diploid releases

Humans are **diploid**: every cell carries two copies of each autosome, one from each parent. To be precise about the counting, because it is easy to garble: a human karyotype has **23 chromosome pairs**, but there are **24 distinct chromosome sequences** (22 autosomes plus X plus Y). The often-quoted figure of about **3.1 billion base pairs is the haploid genome size** — one copy of each chromosome; a diploid cell contains roughly twice that much DNA. For our purposes the consequence is simply that the target supplies two paths, p = (p⁽¹⁾, p⁽²⁾), not one.

The proposal's safest formulation treats the pair as a single private input with joint utility

u_diploid(p, y) = [u(p⁽¹⁾, y⁽¹⁾) + u(p⁽²⁾, y⁽²⁾)]/2.

Averaging two numbers in [0, 1] stays in [0, 1], so condition 3 still holds and **the same tau guarantee covers the complete diploid release** — not tau per haplotype, tau for the pair. This matters: releasing two haplotypes each with its own tau budget would *not* give tau for the pair, for the same compositional reason that two independent releases of one genome do not (condition 4).

### 4.10.4 An honest tension: the budget is thin

Here is arithmetic worth confronting early. The per-block tilt factor is at most exp(eta_tau · beta_t), and the beta_t must sum to 1 across all T blocks. If we weight blocks uniformly, beta_t = 1/T, and on our rebuilt chr21 with T = 100,757 (measured in this project; `docs/memory.md`, 2026-09-18):

| tau | eta_tau | Uniform per-block max tilt exp(eta_tau/T) |
|---|---|---|
| 0.10 | 0.1003 | 1.0000010 |
| 0.50 | 0.5493 | 1.0000055 |
| 0.90 | 1.4722 | 1.0000146 |

At that level the tilt is numerically indistinguishable from nothing: the output is, to six decimal places, a pure draw from R_D. Concentrating the weight helps — with beta spread uniformly over only m blocks and zero elsewhere, the per-block factor is exp(eta_tau/m), which at tau = 0.5 is 1.7321 for m = 1, 1.0565 for m = 10, and 1.0055 for m = 100 — but concentration means abandoning target-specific utility everywhere else.

This is not a flaw in the proof; it is the honest price of a **whole-chromosome worst-case** guarantee, which must hold even for the adversarially chosen pair of inputs that differ at every single block. Several **open questions** follow directly, and they are research questions for this project rather than settled facts:

- Should the guarantee be stated per chromosome, per region, or per some smaller unit, and what is the right composition story across units?
- Should beta_t be concentrated on a small set of clinically or analytically important blocks, and if so chosen how — noting that any *target-dependent* choice would violate condition 2 and void the theorem?
- What tau actually yields acceptable downstream utility (read mapping, genotyping, allele-frequency fidelity) on real data? **We have not measured this.** No tau has yet been run end to end in this project.

### Step table — Part A, the mechanism end to end

| Step | Input | What happens | Output | Why it matters |
|---|---|---|---|---|
| A0. Choose tau | A policy decision | Pick tau in [0, 1); compute eta_tau = arctanh(tau) and epsilon_tau = 2·arctanh(tau) | Two scalars | Fixes the privacy/utility dial *before* any data is seen; the whole theorem hangs on tau being chosen in advance |
| A1. Build the baseline R_D | Public graph G, public cohort D | Fit a cohort haplotype HMM: initial distribution rho, stay/switch transitions A_t, local traversals h_{t,j} | A target-independent generative model over T blocks with K states | Supplies the "menu" of plausible human paths; target-independence is what makes R_D cancel in the likelihood ratio (§4.7, Step 3) |
| A2. Fix the output support Y_D | The HMM from A1 | Declare the set of emittable paths — cohort haplotypes, HMM mosaics, graph-valid traversals — without consulting the target | A fixed support | Condition 1; a target-dependent support sends TV to 1 |
| A3. Map the target | External genome g, graph G | p = Map(g, G) | The private input path p, split into blocks p_1..p_T | The only private object in the system |
| A4. Compute potentials | p, the traversals h_{t,j}, weights beta_t, similarity phi_t | psi_t(j) = exp(eta_tau·beta_t·phi_t(p_t, h_{t,j})) for every block t and state j | A T-by-K table of numbers in [1, exp(eta_tau·beta_t)] | This is the *only* place the private target touches the model, and its influence is bounded by construction |
| A5. Forward filtering | rho, A_t, psi_t | Left-to-right recursion with per-block renormalization; the stay/switch collapse makes each block O(K) | All T forward vectors alpha_1..alpha_T, plus log Z_p = sum of log c_t | Computes the exact conditionals needed for exact sampling; storing *all* alphas is what makes a backward pass possible |
| A6. Backward sampling | alpha_1..alpha_T, transitions | Draw Z_T from alpha_T; then for t = T−1 down to 1 draw Z_t from weights alpha_t(i)·A_{t+1}(i, Z_{t+1}) | One exact draw Z_1..Z_T from Q_p | Sampling (not argmax) is what keeps the output distribution non-degenerate and therefore bounded in TV |
| A7. Emit the path | Z_1..Z_T, traversals h_{t,j} | Concatenate: p̃ = h_{1,Z_1} ‖ … ‖ h_{T,Z_T} | One sanitized path | Deterministic post-processing cannot increase TV, so the guarantee carries through |
| A8. Cache and serve | p̃ | Store the single released path; serve it for all future queries about this genome | A frozen release | Condition 4; independent re-releases would multiply the attacker's evidence, with only a loose k·tau bound to fall back on |

### In plain words (Part A)

We want a privacy guarantee that does not depend on guessing which attacker will show up, so instead of testing attacks we reason about probability distributions. The mechanism never edits the target's genome path; it draws an entirely new path from a distribution over paths that the public cohort already supports. That distribution has two parts: a public baseline model, R_D, built only from the public graph and cohort, and a bounded utility score u(p, y) in [0, 1] saying how useful a candidate output y would be for the target's real path p. We multiply the baseline probability of each candidate by exp(eta_tau · u), where eta_tau = arctanh(tau) and tau is a privacy dial we choose in advance, then renormalize. Because u can only range over [0, 1], the bonus factor can only range over a factor of exp(eta_tau), which caps how much any private input can bend the public baseline. That cap turns into a total variation distance of at most tau between the output distributions of *any* two possible inputs, which in turn caps any attacker's distinguishing accuracy at (1 + tau)/2 — so tau = 0.1 means at most 55% accuracy on a task where guessing gets 50%. Implementing it is standard machinery: the baseline is a Li–Stephens-style copying HMM over T blocks and K cohort haplotype states, the tilt becomes a per-block bonus factor, and we sample exactly with forward-filtering/backward-sampling in O(T·K) time — about 8.9 million operations on our rebuilt chr21, where T = 100,757 block-dictionary entries and K = 88 haplotypes. The proposal lists three conditions the proof depends on — the output menu must be fixed without looking at the target, the baseline must be learned from public data only, and the utility must be bounded — and we add a fourth operational one, that only one release per genome may be made and that we must *sample* rather than take the most likely path, because a deterministic mechanism has total variation distance 1 and therefore no guarantee at all. The honest catch is that a whole-chromosome worst-case budget spread over 100,757 blocks leaves almost no tilt per block, so how to allocate that budget is an open research question rather than a solved one, and no value of tau has yet been run end to end here.

---

## Part B — Comparison with PanMixer

### 4.11 The conceptual comparison

| Dimension | PanMixer (Section 3) | Our mechanism |
|---|---|---|
| **Who is protected** | An individual who **is a contributor** to the released pangenome graph | An **external** individual, not in the cohort and not in the graph |
| **What is private** | That contributor's haplotypes, which are already inside the public graph | The target's mapped path p through an unchanged public graph |
| **What is released** | The **cohort VCF with the contributor's genotype column rewritten** — i.e. a modified reference resource | A **standalone sanitized path**; nothing is ever written back into the reference |
| **Release rule** | Per block, generate one candidate; run a knapsack-style linear-programming optimization to decide replace-or-keep; apply the selected replacements | Resample the **entire path** from a tilted distribution; no per-block keep/replace decision exists |
| **Blocks not selected** | Released **verbatim** | There is no such category — every block is resampled |
| **Privacy quantity** | Per-block Pointwise Mutual Information (PMI) as named in the code, which for an obfuscated block reduces to the self-information −log p(block); "high score = high risk" | tau, a worst-case bound on the total variation distance between output distributions |
| **Type of guarantee** | Average-case, empirical: a measured per-block privacy score plus demonstrated failure of specific attacks (gap-score linkage; membership inference in the published version) | Worst-case, analytic: holds for all admissible input pairs and all attackers, given the conditions of §4.8 |
| **The control knob** | A utility-loss **capacity** (the fraction of a utility-loss denominator the user is willing to spend) | tau, chosen in advance, with the reading "attacker accuracy ≤ (1 + tau)/2" |
| **Role of randomness** | Randomness is present (candidate sampling) but **no formal claim is proved from it**. Note this is *not* the same as saying the seed does not matter: seed secrecy remains operationally necessary for PanMixer, because an attacker who knew the seed could replay the algorithm on candidate inputs until the output matched. Section 2 makes that operational argument; this row is only about what the *formal* claims rest on | Load-bearing in both senses: the randomness *is* the guarantee, and the draw must stay secret |
| **Where the private input enters** | Through the per-block privacy score (the scorer reads the target's own haplotype) and hence the knapsack objective and selection, and through the stacker, which copies the target's unselected blocks. Note that in the **released code** the candidate generator itself is target-independent given the panel: `sample_block_prior` is an unconditional prior draw over the without-target panel, with no comparison to the original and no rejection loop. The published paper describes a constraint the released code does not implement (measured in this project by reading the pinned checkout; see `docs/memory.md`, 2026-09-18) | Only the bounded local potentials psi_t |
| **Operating point** | Chosen empirically: the published paper reports eps_private = 0.002, defined there as the **largest** privacy-risk value at which all targets could no longer be linked (the preprint's wording was the reverse; use the published definition) | Chosen a priori as a policy parameter |

**The two privacy quantities are formally incomparable.** PanMixer's per-block score is an *average-case information quantity* attached to a pair of haplotypes; our tau is a *worst-case bound on a whole output distribution*. Neither implies the other, and there is no conversion formula. A head-to-head comparison between the two systems can only be made on **measured empirical attack success** — running the same corrected attack against both outputs and comparing — never by comparing a privacy score to tau directly.

### 4.12 Why PanMixer's design cannot satisfy our guarantee

This is not a criticism of PanMixer; it is solving a different problem and never claimed a pairwise TV bound. But the reason is instructive, and it is exactly the reason our release rule looks the way it does.

Recall Section 3's **stacker** step: it copies the target's original haplotypes, then overwrites only the blocks the optimizer selected. Every unselected block is republished **byte-for-byte as it was in the target**. We verified this directly in the pinned checkout: `stack()` starts from `original_haplotypes.copy()` and writes into it only at positions where the solution vector equals 1.

Now run the argument. The output alphabet is paths over all T blocks. Consider a block b that the optimizer does not select — and there are many. The optimizer's decision vector has one variable per (block, haplotype) pair, so on our chr21 it has 100,757 × 2 = **201,514** decisions. At capacity 0.1 the shipped run selects **47,070** of them, i.e. about 23%, leaving roughly **77%** of the target's block-haplotypes republished verbatim; the stacker's own accounting on the same run reports **68,335 of 650,541 allele slots changed (10.50%)**, so about 89.5% of individual allele slots come through untouched. (Measured in this project, one subject HG00438 on chr21, rebuilt on the 30x GRCh38 panel; raw figures in `logs/pm_obf.log`, summarized in `docs/memory.md` under 2026-09-18.)

Take two admissible inputs p and q that are identical everywhere except at block b, where they carry different alleles, and suppose b is unselected under both — easy to arrange: pick a block of common alleles with negligible privacy score, which the optimizer will never prioritize.

Then:

- Every output in the support of Q_p carries p's alleles at block b, with probability 1.
- Every output in the support of Q_q carries q's alleles at block b, with probability 1.
- The two supports are **disjoint**.
- Therefore TV(Q_p, Q_q) = **1**, its maximum possible value.
- Therefore the attacker bound (1 + TV)/2 = **1**: an attacker who simply reads block b identifies the input with certainty.

No tuning of the capacity fixes this. Raising the capacity selects more blocks but never all of them, and in fact the knapsack **saturates below its budget**: at capacity 0.5 the same run makes 177,865 moves for 33.1% utility loss rather than the requested 50%, because every move with positive value already fits (measured in this project; `logs/pm_obf.log`, `docs/memory.md` 2026-09-18). And even at full capacity the *conditional* structure remains: which blocks get replaced is itself a function of the target, so the selection pattern leaks.

The general lesson, and the design rule we derive from it: **any mechanism that passes part of its private input through unchanged has total variation distance 1 between some pair of inputs.** A bounded-TV mechanism must resample everything. That is why our release rule cannot be "edit some blocks" and must be "draw a whole path."

**The symmetric honesty.** Our mechanism cannot produce PanMixer's deliverable either. PanMixer's output is a *reference resource* — a pangenome VCF or graph that third parties will genotype against — and it has to keep 44 individuals' columns mutually consistent inside one file. Our output is a single sanitized path for one external individual. Neither tool subsumes the other; they address related but distinct scenarios, exactly as the proposal's Discussion states.

### 4.13 The reuse map

Keyed to the pipeline stages walked through in Section 3. Verdicts: **REUSE AS IS** / **ADAPT** / **REPLACE** / **EVALUATION ONLY**.

Two things to keep straight throughout this table, because several numbers depend on them:

- **Two different 1000 Genomes releases are in play.** The shipped pipeline downloads 1000 Genomes **Phase 3** (2,504 samples, build **GRCh37**, release 20130502). This project substituted the **30x GRCh38** panel (3,202 samples; 1,002,753 chr21 records in the copy registered in `docs/data.md`). Every panel-derived number below says which panel produced it.
- **Block counts are entries, not variants.** On our rebuilt chr21 the block dictionary has **100,757 entries** covering **340,824 variants**. Of those entries, **89,087 (88.4% of entries)** are singletons holding exactly one variant each. Those singletons account for 89,087 of 340,824 variants = **26.1% of variants**. Separately, **10.4% of entries** carry two or more anchor variants and therefore take the HMM code path, and those entries hold **71.4% of variants**. The 100,757 entries are an expansion of **14,137 LD blocks produced by PLINK's `--blocks` routine**. Never say "88.4% of variants": always name the denominator. (All measured in this project; `docs/memory.md` 2026-09-18, and re-verified by counting `blocks_dict.json` directly.)

#### Preprocessing

| Stage (Section 3) | Verdict | Reason and required change |
|---|---|---|
| Graph → VCF (`vg deconstruct -a`, GRCh38 backbone, all nested variants retained) | **REUSE AS IS** | Identical data model. The 340,824 chr21 records are our starting point too |
| VCF → haplotype matrix (`pangenome.npy`, positions, subjects; drop CHM13 and chrX; 44 individuals = 88 haplotypes) | **REUSE AS IS** | This is exactly the K = 88 state set our HMM needs. CHM13 is dropped because it is a reference assembly from a human cell line rather than a cohort donor, not because it is non-human |
| Allele-frequency table (`allele_frequencies.npy`) | **ADAPT** | Must be built against a **build-matched (GRCh38)** panel. The released pipeline matches against the **GRCh37 Phase 3** panel: on chr21 that gives **1,025** exact (POS, REF, ALT) matches out of 340,824, versus **258,610** against the **30x GRCh38** panel — a 252-fold gap. This is defect (a). (Measured in this project; `docs/bugs.md`, `docs/data.md`, `docs/memory.md` 2026-09-18) |
| LD block boundaries (`plink --blocks`) | **ADAPT** | Same build-mismatch fix: PLINK runs on the 1000 Genomes panel, so the boundaries inherit whichever panel is wired in. Boundaries determine T, and T determines the beta_t budget of §4.10.4 |
| Block dictionary / segmentation (`blocks_dict.json`) | **REUSE AS IS structurally, recompute** | The structure — 14,137 PLINK LD blocks expanded to 100,757 dictionary entries, of which 88.4% of **entries** are singletons — is exactly the T-block segmentation we want. Recompute it downstream of the corrected panel |
| Anchor mapping (`pangenome_to_thousand_g_alignments.pickle`) | **ADAPT** | In the released code the variable named `thousand_g_alignments` is bound to `PG.vcf.gz`, i.e. the **PanGenie** callset, not the 1000 Genomes set the paper defines as anchors. This is defect (d): the name is a trap. Because PanGenie is GRCh38, the anchor set itself is build-consistent and the HMM-vs-allele-frequency routing is *not* degraded by the build mismatch; we must simply choose the anchor definition deliberately and name it correctly |

#### Release (this is where the mechanism diverges)

| Stage (Section 3) | Verdict | Reason and required change |
|---|---|---|
| Per-variant utility weight, `utility_loss[v] = 1/support(v)` | **ADAPT** | The raw quantity is an *unbounded cost*, but condition 3 needs a *bounded similarity*. The natural conversion: set W_t = sum over variants v in block t of 1/support(v); then **beta_t = W_t / (sum over all blocks of W_t)** and **phi_t(p_t, y_t) = [sum of 1/support(v) over variants in block t where y agrees with p] / W_t**. *In words:* beta_t is block t's share of the total utility weight, and phi_t is the fraction of block t's own weight that the output preserves. Both live in [0, 1] and the beta_t sum to 1 by construction, so u = sum beta_t·phi_t is exactly "fraction of PanMixer's own utility weight preserved", and it is in [0, 1] by construction rather than by inspection |
| Per-block privacy score (`forward_algorithm` self-information, with an allele-frequency fallback for singleton and zero-anchor blocks) | **REPLACE** | Our guarantee needs no per-block privacy score at all — the privacy is in the calibration, not in a measurement. This also sidesteps defect (c) entirely. **Defect (c), stated carefully:** the scorer scores the target against a panel that *includes* the target (`forward_algorithm` reads `self.pangenome_haplotypes`, all 88 rows, while normalizing transitions by the 86-row without-target count). Because the target always matches itself, the model can always explain the observation by copying the target: probability at least 1/(2N) = 1/88 for picking that state, times (1 − mismatch) at each anchor, where the code's per-anchor mismatch probability is 1e-4. So the score −log p is bounded by about **log(2N) = log(88) = 4.4773 nats plus a small additive term proportional to the number of anchors in the block** — it is not a hard cap at 4.4773. That is why our measured maximum of **4.508** can slightly exceed log(88) without contradicting the argument. **Measured in this project** by executing the shipped `HaplotypeHMM` class on real chr21 data for subject HG00438 over 150 real HMM-path blocks: shipped scores span 0.071–4.508 with median 1.242; re-scored with the target removed from the emission panel, the median is 1.299 but the maximum rises to **78.1**, and for the 9 blocks sitting at the ceiling the leave-one-out values have median 12.06 and maximum 78.12, i.e. understated by up to 17.4-fold. So the effect is a ceiling on the rarest blocks, not a flattening of all of them. Recorded in `docs/memory.md` under 2026-09-18 (Correction 3) and in `docs/bugs.md`. Note that in *our* setting the target is external, so this failure mode cannot arise by construction |
| Candidate resampling (`sample_block_prior`, one candidate per block) | **REPLACE** | Replaced by the tilted FFBS over the whole path. Also note the released sampler is effectively inert. Four compounding defects, collectively defect (b): `EFFECTIVE_N` is set to 1/10,000, the **reciprocal** of the paper's Ne = 10,000; the anchor "distance" fed to the transition model is a difference of VCF **row indices** rather than centiMorgans; the forward pass and the sampler use distance values a factor 4 apart (the published Fig. 6 has no factor 4, so the sampler is the deviation); and the non-anchor fill keys a dictionary by row index while iterating base-pair positions, so its "is this an anchor" guard never fires. **Measured in this project:** on **300 of 300 randomly chosen HMM-path blocks** on real chr21 for subject HG00438, a sampled block equalled **one donor haplotype copied verbatim** — the mosaic never recombines. Recorded in `docs/memory.md` under 2026-09-18 and in `docs/bugs.md` |
| Stay/switch forward recursion (algebra) | **REUSE THE ALGEBRA, ADAPT THE CONSTANTS** | The O(K)-per-block collapse of the stay/switch recursion is precisely what makes our sampler O(T·K), and we inherit it. But the constants must be **re-derived from the published methods, not copied from the code**. The published specification is: effective population size Ne = 10,000; the standard recombination-rate constant **r = 1.26** (this r is unrelated to the LD correlation r² of Section 1); Δx the inter-site distance **in centiMorgans**, read off a genetic map; **d = Δx · Ne · r**, where *d is the expected number of ancestral recombination events separating the two sites*; then stay probability **p = e^(−d/n) + (1 − e^(−d/n))/n** and switch probability **q = (1 − e^(−d/n))/n**, where **n is the number of donor states in the copying model** — 86 in PanMixer's code, which excludes the target's own two haplotypes, and K = 88 for us, since our target is external. *In words:* with probability e^(−d/n) no recombination happened and we stay on the current donor; otherwise we re-draw a donor uniformly from all n, which lands back on the current one 1/n of the time. (Check: p + (n−1)q = 1.) **Beware the collision:** this p and q are transition probabilities, not the input paths p and q of §4.2 and not the allele frequencies p and q of the Hardy–Weinberg model in Section 3. Re-deriving these requires a genetic map to convert base pairs to centiMorgans, and requires dropping the code's spurious factor of 4 |
| Forward pass storage | **ADAPT** | The released `forward_algorithm` keeps only a running `log_alpha`, overwrites it at every step, and returns a single scalar log-likelihood. A backward sampling pass needs **every** alpha_t. The change is to store all T forward vectors: 100,757 × 88 × 8 bytes ≈ 71 MB on chr21, entirely affordable |
| Backward sampling pass | **NEW — nothing to reuse** | There is **no backward pass and no FFBS anywhere in the PanMixer codebase** (checked across the pinned checkout). This must be written from scratch |
| Knapsack / LP relaxation (OR-Tools GLOP, continuous variables in [0, 1] then rounded) | **REPLACE** | We have no selection step. tau replaces the capacity knob entirely. As a by-product we shed both the LP-relaxation-then-round approximation and the observed saturation behaviour of §4.12 |
| Stacker (apply selected moves, keep the rest verbatim) | **REPLACE** | This is the precise step that makes the disjoint-support argument of §4.12 bite. It cannot survive in any form |
| VCF writer (rewrite the subject's genotype column in the cohort VCF) | **REPLACE / REDIRECT** | Our deliverable is a standalone sanitized path for an external individual, not a column inside the cohort file. The output plumbing changes accordingly |

A note on the single-file standalone, since it is the most convenient entry point into their code. `tools/panmixer/obfuscate.py` is **not** a partial reimplementation of a numbered subset of pipeline steps: its own docstring says it "Runs the full pipeline: PMI/utility computation -> optimize -> stack -> VCF output" for one subject on one chromosome, and it **requires pre-existing inputs** including `blocks_dict.json` and `allele_frequencies.npy`. So it consumes the preprocessing artifacts rather than producing them. It is the file we executed to obtain the measurements cited above.

#### Evaluation

| Stage (Section 3) | Verdict | Reason and required change |
|---|---|---|
| Gap-score linkage attack | **EVALUATION ONLY** — must be fixed | This is our best instrument for the only legitimate head-to-head comparison (measured attack success). But the released implementation takes a maximum over the **entire** attack database with no exclusion of the target, while the published equation specifies a maximum over **non-target** individuals. That is defect (e), and it must be corrected before any comparison |
| Diploid gap score | **EVALUATION ONLY** | Directly relevant given the diploid formulation of §4.10.3 |
| Membership-inference analysis (added in the published version; absent from the preprint) | **EVALUATION ONLY** — semantics change | Membership inference asks "was this individual in the training set?" For us, the answer is **no by construction** — the target is external. **Open question:** what the correct analogue is in our setting, and whether a held-out-target protocol would make the two systems comparable on this axis |
| Allele-frequency loss and LD loss | **EVALUATION ONLY** — adapt site sets | Both are exactly the right utility instruments, and they are the axis on which the published paper reports its headline fidelity numbers: at eps_private = 0.002, Wasserstein AF divergence **0.006 for all alleles, 0.006 for SNPs only, and 0.005 for rare SNPs (population MAF < 0.05)**, against a "Removed" baseline — target individuals excluded from the graph entirely — of **0.027**; pairwise LD off by 0.003 on average versus 0.017 for the removal baseline. Always attach those labels; the three AF numbers are not interchangeable. In the released pipeline the site sets for these metrics are drawn from the build-mismatched **GRCh37 Phase 3** panel and must be recomputed on a build-matched one. **We have reproduced none of these published numbers** |
| Beagle refinement and imputation-accuracy statistics | **EVALUATION ONLY** | Reusable essentially as-is as a downstream utility probe. The repository pins the jar `beagle.27Feb25.75f.jar` — refer to it as the **27Feb25 build**; do not attach a marketing version number to it |
| Read mapping (`vg_prep`, `quick_align`, filtered and personalized graph mapping with Giraffe) | **EVALUATION ONLY** — and more central for us | Since our deliverable *is* a mapped path, mapping-quality metrics are arguably the most important utility signal we have. The published baseline to beat, on the full graph: perfectly aligned reads **77.82%** for PanMixer versus **77.83%** for the unmodified graph, gapless **95.61%** versus **95.62%**, and MAPQ 60 **77.01%** versus **77.01%**. Again, unreproduced here |
| Plotting and results aggregation (`plot_util/`, `gather_results`) | **EVALUATION ONLY** | Reusable essentially wholesale |

**Overall reuse assessment — stated as an estimate, because it is one.** The pinned checkout contains about **8,661 lines of Python** (counted across `.py` files in the checkout). Earlier drafts put the reusable fraction at "55–65%" and the reusable share of the *mechanism* at "under 10%". Those percentages have **no line-by-line counting method behind them** and should be read as rough impressions, not measurements. The qualitative claim is the one we can stand behind: **the data model, the preprocessing chain, the block segmentation, the stay/switch forward algebra and essentially the whole evaluation suite are reusable; the sampler, the privacy score and the release rule are not**, and the backward sampling pass at the heart of our method does not exist anywhere in their codebase. **Open item:** if a precise figure is wanted, it needs an actual per-file accounting.

### 4.14 The four concrete adaptations, spelled out

Collecting the engineering consequences in one place:

1. **Recompute every artifact that currently includes the target.** PanMixer's scorer is bounded near log(88) because the target sits inside its own reference panel. In our setting the target is external, so every cohort-derived artifact — allele frequencies, support counts, the HMM state set, the attack database — is target-free by construction; but the code paths that quietly re-include the target must be audited and rebuilt rather than trusted.
2. **Re-derive the transition constants from the published parameters, not from the code.** Use Ne = 10,000; the recombination-rate constant r = 1.26 (again, unrelated to LD r²); Δx in centiMorgans from a real genetic map; d = Δx · Ne · r, the expected number of ancestral recombinations between the two sites; stay probability p = e^(−d/n) + (1 − e^(−d/n))/n and switch probability q = (1 − e^(−d/n))/n, with n the number of donor states (K = 88 for us). Do not carry over `EFFECTIVE_N = 1/10,000`, the row-index "distances", or the factor of 4. Without this the chain never recombines — which we measured on 300 of 300 blocks — and R_D degenerates to "copy one donor haplotype for the whole chromosome", which would make our tilt meaningless just as it makes theirs inert.
3. **Normalize the per-variant utility weight per block into a bounded similarity.** Convert 1/support(v) into the (beta_t, phi_t) pair defined in §4.13 so that u lands in [0, 1] by construction rather than by inspection. This is condition 3, and it is the one condition an implementer can silently break: a utility function that is merely "usually small" does not satisfy it, and the theorem gives nothing at all if the range is unbounded.
4. **Extend the forward pass to store all its intermediate vectors.** The existing forward algorithm discards alpha_t as it goes because it only needs a scalar likelihood. Backward sampling requires all of them. This is a small change with a large consequence: it is the single code change that converts PanMixer's scorer into the first half of our sampler.

### 4.15 Side-by-side flow

```
PANMIXER (Section 3)                       OURS (Section 4)
====================                       ================

public graph G  (target IS inside)         public graph G  (target is OUTSIDE)
        |                                          |
   [ deconstruct to VCF ]  -- shared ----->  [ deconstruct to VCF ]
        |                                          |
   [ LD blocks, AF table ] -- shared ----->  [ LD blocks, AF table ]
        |          (both rebuilt on the build-matched 30x GRCh38 panel)
        |                                          |
   [ per-block privacy score ]               [ per-block potential psi_t ]
        |  target IS scored, against              |  target only re-weights,
        |  a panel containing itself              |  and by a bounded factor
        v                                          v
   [ sample 1 candidate per block ]          [ FORWARD FILTER, keep all alpha_t ]
        |  (unconditional prior draw,              |     -- shared algebra --
        |   no tilt; inert as shipped)             |
        v                                          v
   [ KNAPSACK: replace or keep? ]            [ BACKWARD SAMPLE whole path ]
        |  capacity budget                         |  tau budget
        v                                          v
   [ STACK: unselected blocks VERBATIM ]     [ concatenate h_{t,Z_t} ]
        |  ~77% of block-haplotype                 |
        |  decisions on our chr21 run              |
        v                                          v
   cohort VCF with target's column           ONE standalone sanitized path
   rewritten                                 (graph untouched)
        |                                          |
   TV between two inputs = 1                 TV between ANY two inputs <= tau
   (disjoint supports on any                 attacker accuracy <= (1+tau)/2
    unselected differing block)
        |                                          |
        +---- shared evaluation suite -------------+
              (gap score, AF/LD loss, Beagle
               27Feb25 build, read mapping,
               plotting)
```

### Step table — where each of PanMixer's steps goes

| Step | Input | What happens | Output | Why it matters |
|---|---|---|---|---|
| Preprocessing (graph → VCF → matrix → blocks) | HPRC PGGB graph | Reused essentially unchanged; allele-frequency and block artifacts rebuilt against the build-matched 30x GRCh38 panel instead of the shipped GRCh37 Phase 3 one | 88 cohort haplotypes over 340,824 chr21 records, segmented into 14,137 PLINK LD blocks expanded to 100,757 dictionary entries | Gives us K = 88 states and T = 100,757 blocks for free; the build fix is a prerequisite for both projects |
| Per-variant utility weight | Haplotype matrix | 1/support(v) is aggregated per block into W_t, normalized into beta_t (block t's share of total weight) and phi_t (fraction of block t's weight preserved) | u(p, y) in [0, 1] | Condition 3; without boundedness there is no calibration and no theorem |
| Per-block privacy score | Target haplotypes, panel including the target | **Deleted.** Privacy comes from calibration, not measurement | — | Removes the panel-includes-target failure mode, whose measured effect is a ceiling near log(88) that truncates the rarest blocks by up to 17.4-fold |
| Candidate sampling | Cohort haplotypes, anchors | **Replaced** by a tilted forward pass; transition constants re-derived from the published Li–Stephens parameters with distances in centiMorgans | A T-by-K potential table plus all T forward vectors (~71 MB) | The tilt is the only channel through which the private target influences the output, and it is bounded by construction |
| Knapsack selection | Block scores and costs, capacity | **Deleted.** tau replaces capacity as the control knob | — | Removes both the LP-relaxation approximation and the structural source of disjoint supports |
| Stacking | Target haplotypes, selected moves | **Replaced** by backward sampling of the entire state sequence and concatenation of local traversals | One exact draw from Q_p, post-processed into one path | The replaced step is precisely the one that makes TV = 1; whole-path resampling is what makes TV ≤ tau achievable |
| Output writing | Sanitized path | **Redirected**: emit a standalone path rather than rewriting a column in the cohort VCF | The release | Different deliverable for a different threat model; the public reference is never modified |
| Evaluation | Both systems' outputs | Reused nearly wholesale, with the gap-score attack corrected to exclude the target from its maximum | Comparable empirical attack-success and utility numbers | Empirical attack success is the *only* axis on which PanMixer's privacy score and our tau can be compared |

### In plain words (Part B)

PanMixer and our mechanism look superficially similar — both cut a chromosome into linkage-disequilibrium blocks, both use a Li–Stephens haplotype hidden Markov model over the same 88 cohort haplotypes from 44 individuals, both trade privacy against utility — but they protect different people against different threats and they end at completely different places. PanMixer protects someone who is *inside* the public graph, and its output is the cohort file with that person's genotype column rewritten; ours protects an *external* person's mapped path and never touches the public reference at all. PanMixer picks one candidate per block and then runs a knapsack-style optimization to decide which blocks to replace, leaving the rest exactly as they were; ours has no selection step and resamples the whole path from a tilted distribution. That one difference is decisive: because PanMixer republishes unselected blocks verbatim — about 77% of the 201,514 block-by-haplotype decisions on our chr21 run at capacity 0.1, corresponding to about 89.5% of individual allele slots untouched — two inputs that differ at an unselected block produce outputs whose supports do not overlap at all, so their total variation distance is 1 and no bounded-distinguishability guarantee is possible. That is not a bug in PanMixer, which never claimed such a guarantee; it is a structural fact about any mechanism that passes part of its private input through unchanged, and it is the reason our release rule has to resample everything. The two privacy numbers are also formally incomparable — PanMixer's per-block score is an average-case information quantity and our tau is a worst-case bound on a whole output distribution — so any head-to-head must be run on measured attack success, using a corrected version of their gap-score attack that excludes the target from its maximum. In terms of code, the pinned checkout is about 8,661 lines of Python, and while the data model, preprocessing, block segmentation, stay/switch forward algebra and almost the whole evaluation suite are reusable, the sampler, the privacy score and the release rule are not, and the backward sampling pass at the heart of our method does not exist anywhere in their codebase; any percentage attached to those statements would be a guess until someone does the per-file accounting. The four concrete adaptations are: rebuild every artifact so it never contains the target, re-derive the transition constants from the published Li–Stephens parameters with distances in centiMorgans rather than copying the code's constants, normalize the per-variant utility weight into a bounded per-block similarity so that u stays in [0, 1] by construction, and extend the forward pass to keep all 100,757 of its intermediate vectors so a backward sampling pass becomes possible. Do those four things and PanMixer's biologically grounded cohort model becomes our target-independent prior, with a formal pairwise identifiability guarantee bolted on top — a guarantee that, it must be said, we have specified and proved on paper but not yet run on data.