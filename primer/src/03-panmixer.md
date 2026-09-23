## Section 3 — PanMixer: the published tool, dissected

Everything in this section is about one specific tool: **PanMixer**, published in *Nature Communications* (DOI 10.1038/s41467-026-77591-0, "Article in Press"), with source code released under the MIT licence.

Two sources, kept strictly apart:

- Where I quote a **design decision**, I take it from the **published** paper. Our copy of record is `paper/s41467-026-77591-0_reference.pdf` in this project's repository. An earlier **preprint** exists (`archive_docs/Blindenbach2026_PanMixer.pdf`); it is superseded and several of its numbers were revised, so I flag it whenever the two disagree and never quote it as current.
- Where I describe **a file, an array shape, or a line of logic**, I take it from the released code, pinned at commit `c182c38d5bc8bb6f00f4b0b101207c4a009ca045` in this project's `external/PanMixer` checkout.

Those two things are not always the same, and the subsection "Where the implementation departs from the paper" is entirely about where they diverge.

**A note on provenance, which matters more here than anywhere else in this primer.** Numbers in this section come from three different places, and I label every one:

| label | meaning | where it is recorded |
|---|---|---|
| **published** | printed in `paper/s41467-026-77591-0_reference.pdf` | that PDF |
| **measured in this project** | produced by us, by executing the pinned code or inspecting our rebuilt data | `docs/memory.md` (the permanent dated ledger), `docs/data.md` (the artifact registry), `docs/bugs.md` (defects), and the run logs under `logs/` |
| **estimate** | a rough figure with no line-by-line counting method | flagged in the sentence itself |

Anything not carrying one of those labels is either a definition or a statement about code you can read yourself at the cited path. I do not cite scratch or temporary files anywhere, because a primer is a thing you come back to and a temporary path does not survive.

**One structural warning before any symbol appears.** This document uses several Greek letters in more than one sense, because PanMixer's paper and our own proposal (`paper/Private_Genome_Path_Release.pdf`) independently chose the same letters for different things. The table below is the authoritative disambiguation for Section 3 and should have been in front of you from Section 1 onward. Read it now; I will re-gloss each symbol at its point of use as well.

| symbol | meaning **in this section** | meaning **elsewhere** | reminder |
|---|---|---|---|
| **η** (eta) | PanMixer's **utility loss** — how much the release degrades the data | our proposal's **η_τ** is a *tilt strength*, a completely different object introduced in Section 4 | η here costs; η_τ there steers |
| **ε** (epsilon) | PanMixer's **privacy risk**. Two flavours: **ε_j**, the per-block risk, and **ε_private**, an empirical threshold | Section 1's **ε** is a differential-privacy budget; our proposal's **ε_τ = 2·arctanh(τ)** is a likelihood-ratio bound | none of these are interchangeable |
| **ε as a mismatch rate** | in the HMM's emission model only, a per-anchor error probability of 10⁻⁴. I write it **μ_mis** from here on to stop the collision | — | renamed on purpose |
| **p, q** | at a biallelic site, **allele frequencies** (p = reference-allele frequency, q = 1 − p) | in Section 4, **p** and **q** are *input paths*; in the Li–Stephens transition formula they are *stay* and *switch probabilities* | three unrelated uses; always check the sentence |
| **r** | the Li–Stephens **recombination-rate constant, r = 1.26** (published, Fig. 6) | **r²** is the LD correlation statistic of Section 1 | *unrelated quantities that happen to share a letter* |
| **N_e** | **effective population size**, fixed at 10,000 (published) | — | defined in §A.3 |

The section has two parts. **Part A** is the idea, stated without any code. **Part B** is the tool as a data pipeline: every step, its inputs, what happens to them, its outputs, and what you should expect to see if it worked.

---

## Part A — The theory

### A.1 Who is being protected, and from what

PanMixer protects a **contributor** to the pangenome graph.

To unpack that, three definitions, in order.

A **pangenome graph** is a data structure that represents many human genomes at once. Instead of one linear string of DNA letters (a **reference genome** — a single agreed-upon sequence that everyone's coordinates are stated relative to, e.g. GRCh38), you have a graph whose **nodes** are stretches of sequence and whose **edges** connect stretches that occur consecutively in at least one real genome. Each contributing person's chromosome is then a **path** through that graph — the specific sequence of nodes their DNA actually spells out.

A person has two copies of each **autosome** (a non-sex chromosome), one inherited from each parent, so each person contributes **two paths per chromosome**. Each such single-parental-copy sequence is a **haplotype**. When you know which variants sit on which parental copy, the data is **phased**; when you only know the unordered pair, it is **unphased**, and the unordered pair is called a **genotype**.

PanMixer works on the **variant-centric** view of this graph rather than on nodes and edges directly. A **variant** is a site where individuals differ from the reference: a **SNP** (single-nucleotide polymorphism, one letter swapped), an **indel** (a short insertion or deletion), or a **structural variant** (a large insertion, deletion, inversion, and so on). Each possible sequence at a variant site is an **allele**. A site with exactly one alternate allele — two alleles in total — is **biallelic**; a site with two or more alternates is **multi-allelic**. The standard text format for "here are the variant sites and here is which allele each sample carries" is a **VCF** (Variant Call Format) file: one row per variant site, one column per sample, and the genotype field written as `0|1` — allele index 0 (the reference allele) on the first haplotype, allele index 1 (the first alternate allele) on the second, with `|` meaning phased. In this view a haplotype path is just a vector of allele indices, one per variant site.

Now the privacy problem. If you are one of the 44 people whose genomes were used to build the graph, then your two allele vectors are *literally published inside the released VCF*, under your sample name. Anyone who separately obtains a piece of your DNA — a discarded sample, a consumer-genetics upload, a clinical record — can check whether it matches. This is a **re-identification** or **linkage attack**: an attacker holds a database of named genomes and asks which row of the published data corresponds to which person. It works so well because rare variants are nearly unique fingerprints. Sharing a common variant with someone is uninformative, because a large fraction of the planet shares it; sharing three or four rare ones is essentially a signature.

It is worth knowing, as **standard background rather than anything measured here**, that the vulnerability is not limited to per-individual records. The result usually credited to **Homer and colleagues (2008)** showed that even *aggregate* allele frequencies — a table of "this fraction of the cohort carried the alternate allele at each site", with no individual rows at all — can be enough to tell whether a known individual was part of the pool. That result is a standard piece of field background; I state it as attribution, not as something this project verified, and a reader who wants to rely on it should read the original.

PanMixer's answer is: **rewrite the contributor's genotype column** so that it is no longer a faithful copy of the real person, while keeping the released graph scientifically useful. The release is still "the cohort VCF", with the same sites and the same 44 sample columns; only the target's column has been altered.

### A.2 Why LD blocks are the unit of action

PanMixer does not edit one variant at a time. It edits **LD blocks**.

**Linkage disequilibrium (LD)** is the statistical fact that nearby variants are not inherited independently. **Recombination** — the shuffling of the two parental chromosomes during the production of eggs and sperm — happens at particular rates along the chromosome, and the genetic distance over which it happens is measured in **centiMorgans (cM)**: 1 cM is the distance over which a recombination occurs in about 1% of meioses. Note that a centiMorgan is *not* a physical distance: the conversion between cM and base pairs varies along the chromosome, because recombination is concentrated in hotspots. Stretches of DNA shorter than the spacing between hotspots tend to be passed down intact, so the alleles inside them travel together as a unit. An **LD block** is a contiguous run of variants that behaves like such a unit.

This matters for privacy for a blunt reason: **if you edit one variant but leave its LD neighbours alone, you have not hidden anything.** The neighbours still spell out which allele you "really" had, because in this population that combination essentially always co-occurs. Editing must happen at the granularity at which the data is actually correlated, which is the block. It also matters for utility: a half-edited block is a haplotype no human has, which corrupts exactly the LD structure that downstream users want.

PanMixer computes block boundaries with **plink** (a long-established open-source genetics toolkit; its `--blocks` flag runs a confidence-interval LD-block caller) on a large external **reference panel** — a published collection of phased haplotypes from many people, used as a statistical stand-in for "the population". It then assigns every pangenome variant to whichever block interval contains its coordinate. Variants that fall outside all block intervals become **blocks of their own**, holding exactly one variant each.

A number to carry, **measured in this project** on our rebuilt chr21 (recorded in `docs/memory.md`, entry dated 2026-09-18): plink produced **14,137 LD blocks**, but the working block dictionary has **100,757 entries**, of which **89,087 — that is 88.4% of *entries*** — hold exactly one variant. Because singleton entries hold one variant each by definition, those 89,087 entries account for 89,087 of 340,824 variants, or **26.1% of *variants***. **These two denominators must never be confused**: "88.4%" is a statement about entries, and the corresponding statement about variants is "about 26%". Separately, **10.4% of entries** contain two or more **anchor** variants (defined in §A.3) and therefore take the HMM path, and those entries hold **71.4% of all variants**.

The shape to internalise: *most entries are trivial, but most of the genome lives in the non-trivial minority.*

### A.3 What an "obfuscation move" is

For each block, and for each of the target's two haplotypes, PanMixer constructs exactly **one candidate replacement**: a plausible alternative haplotype block, sampled from a model of the cohort. That candidate, together with the decision to install it, is an **obfuscation move**. With *k* blocks and 2 haplotypes there are *2k* possible moves, and each is a binary yes/no.

An **anchor** is a variant the tool treats as trusted enough to model explicitly — concretely, a pangenome variant that also appears in a chosen external callset. (A **callset** is the set of variants a particular method reported from a particular dataset.) Which callset supplies the anchors is exactly where the code and the paper part company; see defect (d).

How the candidate is produced depends on how much structure the block has:

- **Blocks with fewer than two anchors**: sample each variant's new allele independently from that site's **allele frequency** (AF) — the fraction of haplotypes in a reference population carrying that allele. This ignores LD, which is acceptable precisely because a one-variant block has no internal LD structure to preserve.
- **Blocks with two or more anchors**: sample from a **Li–Stephens** model. This is a **hidden Markov model (HMM)** — a probabilistic model with an unobserved ("hidden") state that evolves along a sequence and emits observations. In Li–Stephens the hidden state at each anchor is "which cohort haplotype am I currently copying from", the transition is "stay on the same donor, or switch to another (a recombination)", and the emission is "copy that donor's allele, with a small probability of mismatch". A haplotype produced this way is a **mosaic**: a patchwork of segments copied from different donors. The whole point of the mosaic is that it is a *new* combination — not any single real person's data — while still looking like a real human haplotype, because every piece of it came from one.

The switch rate is where the population genetics enters. The published model (Fig. 6) sets

> **d = Δx · N_e · r**

Reading that in words: **d is the expected number of ancestral recombination events separating two neighbouring anchors**. It is the genetic distance **Δx** between them (in centiMorgans), scaled by how much ancestral shuffling the population has accumulated (**N_e**) and by a standard rate constant (**r = 1.26**, published). **N_e**, the **effective population size**, is the size of an idealised, randomly mating population that would show the same amount of genetic drift — and the same rate of shared ancestry — as the real one; it is much smaller than the census population because of past bottlenecks, and humans are conventionally modelled with **N_e = 10,000** (published). Again: this **r** is the recombination-rate constant and has nothing to do with the LD correlation **r²**.

From d, with **n** donor states available, the two transition probabilities are

> **stay: p = e^(−d/n) + (1 − e^(−d/n))/n**  
> **switch: q = (1 − e^(−d/n))/n**

In words: with probability e^(−d/n) no recombination occurred and we stay on the current donor; otherwise we re-draw a donor uniformly from all n, which lands back on the current donor 1/n of the time. (Here **p** and **q** are transition probabilities, *not* the allele frequencies p and q of §A.7 and Step 21, and *not* the input paths p and q of Section 4.) In PanMixer's sampler **n = 86**, because the target's own two haplotypes are excluded from the donor set. Note that the two lines sum correctly: p + (n−1)q = 1.

Each move buys some privacy — the block no longer reveals the target's real alleles — and costs some utility, because the cohort's allele frequencies and LD structure shift a little.

### A.4 The privacy quantity: PMI, and why it collapses to self-information

PanMixer's privacy number is called **PMI** — **pointwise mutual information** — both in the paper and in the code. Since the primer has so far used the abbreviation without defining it, here is the definition:

> **PMI(x, y) = log [ P(x, y) / (P(x)·P(y)) ]**

In words: **how much more often do x and y occur together than they would if they were independent?** If x and y are statistically independent, P(x,y) = P(x)P(y), the ratio is 1, and PMI = 0. If observing y makes x much more likely than chance, PMI is large and positive. It is the pointwise (per-outcome) version of mutual information, which is the average of PMI over the joint distribution.

PanMixer applies this to the pair (original path *h*, released path *h\**) and calls it privacy risk **ε**. The reasoning is direct: a released path that is statistically *dependent* on the original leaks information about the original; a released path that is *independent* of it leaks nothing. PMI(h, h\*) = 0 is the ideal.

The paper defines the joint probability blockwise over blocks (b₁, …, b_k), which is where the useful collapse comes from. Two cases, and this is the whole argument:

- For a block you **do not** obfuscate, the released block *is* the original block. So the joint probability equals the marginal, p(h_bj, h\*_bj) = p(h_bj), and the block contributes its full PMI.
- For a block you **do** obfuscate, PanMixer samples the replacement *independently of the original*. So the joint factorises, p(h_bj, h\*_bj) = p(h_bj)·p(h\*_bj), and PMI for that block is log 1 = 0.

Therefore obfuscating block *b_j* **removes** exactly the amount of risk that block was contributing, and the paper's Supplementary Lemma 2 gives that removed amount as

> **ε_j = −log p(h_bj)**

where p(h_bj) is the population probability of the target's *own* haplotype block. In words: **how surprising is the target's real block under a model of the general population?**

That quantity is the **self-information** (also called surprisal) of an event. The intuition is clean. If your block is the most common haplotype in the population, p is large, −log p is small, and disclosing it tells an attacker almost nothing — millions of people have it. If your block is one that essentially nobody else carries, p is tiny, −log p is large, and disclosing it identifies you. Because log is monotone, ranking blocks by ε_j is exactly ranking them by rarity. **Rarer blocks score higher and get obfuscated first.**

**Units.** Self-information is measured in **bits** when the logarithm is base 2 and in **nats** when it is the natural logarithm. PanMixer's code uses `np.log`, the natural log, so everything in this section is in **nats**. One nat ≈ 1.44 bits.

**Micro-example.** Suppose a block reduces to a single variant with two alleles, and the reference allele has frequency 0.9, the alternate 0.1. A target carrying the reference allele scores −log(0.9) = 0.105 nats. A target carrying the alternate scores −log(0.1) = 2.303 nats — about 22 times more. Now suppose the target carries a **novel** variant, one seen in no external database at all. Step 13 gives such a variant a frequency with a deliberately large denominator: 2 × (N_panel + N_pangenome). On our rebuilt chr21, using the **30x GRCh38 panel (3,202 samples)**, that denominator is 2 × (3,202 + 44) = **6,492**, so a variant carried on one haplotype scores −log(1/6492) = **8.78 nats**. (Had the shipped pipeline's **Phase 3 panel (2,504 samples)** been used, the denominator would be 2 × (2,504 + 44) = 5,096 and the score 8.54 nats. Always say which panel produced the number.) The optimizer will reach for that block long before it touches the common one.

For multi-anchor blocks, p(h_bj) is not a product of independent allele frequencies — that would double-count correlated variants. It is computed with the **forward algorithm** on the same Li–Stephens HMM. The forward algorithm is the standard dynamic program answering "what is the total probability that this HMM emits the sequence I observed", summing over all hidden state paths in O(T·K) time for T positions and K states, by carrying forward one number per state instead of enumerating K^T paths. Using it here means a block whose alleles are individually uncommon but *jointly* a well-known haplotype gets a correctly modest score.

Two properties are worth naming now, because our own mechanism (Section 4 — the cohort-tilted τ-private path release set out in `paper/Private_Genome_Path_Release.pdf`) departs from both:

1. ε_j is a property of **the target's own data** evaluated under a population model. It is an *average-case information quantity about one observed value* — not a bound on any distribution the mechanism outputs.
2. ε_j is defined *per block*, and the total risk is the sum across blocks. Blocks you do not obfuscate contribute their full ε_j to the residual risk.

### A.5 The utility cost: one over the support

Each move's cost is

> **η_j = Σ_{v ∈ Δ} 1 / support(v)**

where **Δ** is the set of variants the move actually changes and **support(v)** is the number of haplotypes that have a *called* (non-missing) allele at variant v. In words: **for every variant this move rewrites, add one divided by how many haplotypes were actually observed at that variant; the sum is the move's cost.** (Reminder from the symbol table: this η is PanMixer's utility loss, not our proposal's tilt strength η_τ.)

Why 1/support? Because to a good approximation it is exactly the change in allele frequency you cause. Take a single variant with support *n* — that is, *n* haplotypes have a called allele there. Its allele frequency is (count of haplotypes carrying that allele)/n. If you flip one haplotype from reference to alternate, the numerator changes by 1, so the frequency changes by 1/n. The **L1 distance** (sum of absolute differences) between the before and after frequency vectors therefore changes by 1/n for the allele lost and 1/n for the allele gained. So η is, up to a constant factor of two, the induced L1 change in allele frequency, summed over the variants you touched.

Two consequences follow, and both are important:

- **Cost is flat in allele frequency.** Changing a rare variant costs the same 1/n as changing a common one. Privacy *value*, by contrast, scales steeply with rarity. So the value-per-unit-cost ratio is enormously better for rare blocks, and the optimizer naturally strips out exactly the identifying material first. That is the behaviour you want, and it is not an accident of the implementation — it falls out of the asymmetry between a log and a constant.
- **Cost is roughly linear in the number of edits.** A block containing 200 variants of which 40 change costs about 40/n. So "utility loss" and "amount of the genome rewritten" track each other.

The paper calls the aggregate of these costs the **Weighted Path Edit Distance (WPED)**: the minimum total cost of the edits needed to turn the original path into the released path, by analogy with edit distance on strings.

### A.6 The knapsack and the LP relaxation

Now the optimization. You have *2k* items. Item *j* has a **value** ε_j (privacy removed) and a **weight** η_j (utility destroyed). You have a **budget** ΔU on total weight. You must choose a subset.

This is literally the **0–1 knapsack problem**: choose x ∈ {0,1}^{2k} to

> **maximize Σ_j ε_j x_j subject to Σ_j η_j x_j ≤ ΔU.**

In words: **take the most privacy you can buy without spending more utility than your budget allows.**

0–1 knapsack is NP-hard in general, and 2k here is on the order of 200,000 items per chromosome (twice the 100,757 block-dictionary entries measured on our rebuilt chr21). PanMixer therefore solves the **LP relaxation** — "LP" being **linear program**, an optimization with a linear objective and linear constraints, solvable in polynomial time. The relaxation replaces x_j ∈ {0,1} with x_j ∈ [0,1], i.e. it permits "take 37% of this move", and the solution is then rounded back to 0 or 1.

This is not a heuristic hack, and it is worth understanding why. The LP relaxation of a knapsack with a **single** capacity constraint is the *fractional* knapsack, whose optimum is obtained by the greedy rule: sort items by value density ε_j/η_j descending, take them in order until the budget runs out, and take a fraction of the one item that straddles the boundary. **At most one variable is fractional at the LP optimum.** So rounding perturbs the solution by at most one item out of roughly 200,000, and the LP optimum upper-bounds the true integer optimum. In practice the relaxation is essentially exact here, and the published paper reports that its results came from solutions the solver reported as optimal.

The dual framing — fix the privacy level and minimize utility loss — is the same LP with value and weight swapped and the inequality flipped, which is why the tool exposes `--fixed_param utility` and `--fixed_param privacy`.

### A.7 What kind of guarantee this is — and what it is not

This is the single most important conceptual point in the section.

PanMixer's headline privacy number is **ε_private = 0.002** (**published**). The published definition: *the largest privacy risk value at which all obfuscated target individuals could no longer be linked.* (The superseded preprint gave 0.001 and defined it as the *smallest* such value; the published wording is the one to use.) It is found empirically, by running an actual linkage attack against an actual external genotype database at decreasing privacy-risk settings and observing where every target's **gap score** goes negative. Published equation (7) is

> **GapScore(g_target) = L(g_target, g\*_target) − max_{i ≠ target} L(g_i, g\*_target)**

where L is a rarity-weighted matching score between a database entry and the released obfuscated genome. In words: **how much better does the true person match the release than the best impostor in the database?** Positive gap score = the true person is still the best match = attack succeeded. Negative = someone else scores higher = attack failed.

**A caution carried forward from §2.4.1, repeated here because anyone re-implementing from this section alone would otherwise not know it exists.** The published equation (5) as typeset reads L = Σ_{j ∈ S} log f_j, where S is the set of loci at which the two genotypes are identical and f_j is that genotype's population frequency. Read literally, that is a *negative* number that gets *more* negative for rarer matches, which is the opposite of the "emphasizes rare allele matches" reading the surrounding text gives it. The natural reading — and the one the code's behaviour is consistent with — is the negative sum, −Σ log f_j, i.e. a positive rarity-weighted score. **The sign discrepancy is unresolved and is flagged as an open question, not settled here.**

So ε_private is:

- **empirical** — it is measured by running one attack, not derived;
- **attacker-specific** — a threshold against *this* scoring function and *this* attack;
- **database-specific and cohort-specific** — the paper says so explicitly: "ε_private is not universal as it depends on the specific graph, external genotype database, and target individuals";
- **an average-case statement about observed scores**, not a bound on a distribution.

At that threshold the published paper reports an **average utility loss of 0.28**, on a scale whose maximum is 1.

What ε_private is **not** is a worst-case distributional guarantee. In the privacy literature the gold standard is a statement of the form "for any two neighbouring inputs, the output distributions are close", which bounds *every possible attacker at once* — **differential privacy** is the canonical example, and **total variation distance** (the largest difference in probability that two distributions assign to any single event) is the natural way to measure "close". PanMixer makes no such claim, and structurally cannot: because unselected blocks are copied through **verbatim**, two inputs that differ in an unselected block produce outputs that differ there with certainty, so the two output distributions have **disjoint support** and total variation distance 1. That is not a flaw in PanMixer's execution of its own goal; it is a statement that PanMixer's goal is a different one. It is also the precise gap our own mechanism is designed to close, which is Section 4.

One clarification that has caused confusion elsewhere in this primer. PanMixer's mechanism **does** use randomness — the candidate for each block is a random draw — and that randomness is operationally load-bearing: if the seed leaked, an attacker could re-run the sampler on candidate inputs until the output matched. But no **formal** claim of the paper is *proved from* that randomness; ε_private is measured against realised outputs. Both statements are true and they are not in tension: *seed secrecy is necessary in practice; randomness is not what the guarantee rests on, because there is no distributional guarantee to rest.*

### In plain words — Part A

PanMixer is for people whose genomes are *inside* the published pangenome. Their data sits in the released file under their own name, and rare variants make that data a fingerprint. PanMixer chops the genome into LD blocks — runs of variants that are inherited together and therefore have to be edited together or not at all — and for each block it prepares one plausible substitute block, sampled either from population allele frequencies (for blocks with fewer than two anchors) or from a Li–Stephens HMM that stitches together pieces of other cohort members' haplotypes (for blocks with two or more). Each substitution is a "move". A move's *value* is how surprising the target's real block was under a population model, so rare blocks are worth a lot and common blocks almost nothing; a move's *cost* is one over the number of haplotypes observed at each variant it changes, summed up, which is essentially the allele-frequency distortion it inflicts. Because value grows with rarity while cost does not, the best moves are exactly the ones that strip out identifying rare material. Choosing which moves to take under a utility budget is a 0–1 knapsack problem, solved by its linear-programming relaxation, which for a single-constraint knapsack is near-exact — at most one item comes out fractional. The resulting privacy claim is empirical: ε_private = 0.002 is the largest privacy-risk setting at which a specific linkage attack against a specific database failed for every target. It is a measured threshold against one attacker, not a worst-case bound against all attackers.

---

## Part B — The tool, as a data pipeline

The pipeline splits cleanly in two. **Preprocessing** (`starting_data/run_get_starting_data_pipeline.sh`) runs once per pangenome and per target and is the expensive part. The repository's README says so in as many words, at line 78: *"These commands take about 1–2 days to run end-to-end."* **Experiment runtime** (`main.py`) is cheap and is run many times at different budgets.

Numbers below labelled "chr21" are **measured in this project** on our rebuilt chr21 unless stated otherwise, and are recorded in `docs/memory.md` (entries dated 2026-09-18), `docs/data.md` and `docs/bugs.md`, with run logs under `logs/` (`pm_chr21.log`, `pm_obf.log`, `pm_prep1.log`). Array shapes are given because they are the fastest way to check whether a step did what it should.

**Which panel.** This matters throughout. The shipped pipeline downloads **1000 Genomes Phase 3** — 2,504 samples, **GRCh37** coordinates. This project **substituted** the **1000 Genomes 30x GRCh38 panel** — 3,202 samples, 1,002,753 chr21 records — by re-pointing the hardcoded filename slot, so that the pinned checkout runs unmodified on build-correct data. Every panel-derived number below says which panel produced it.

**Vocabulary used from here on, glossed once.** A **read** is a short DNA fragment output by a sequencing machine (typically 100–300 bp). **Alignment** is finding the place in a reference sequence that a read most closely matches, allowing a few mismatches. An **assembly** is a genome sequence reconstructed from sequencing fragments directly, rather than by aligning them to a reference. A **contig** is a named contiguous reference sequence — for a human reference, usually one chromosome; the contig name is the key that every coordinate is relative to, which is why a naming-convention mismatch breaks joins. **`vg`** (short for "variation graph") is the standard open-source toolkit for building, indexing and aligning against pangenome graphs. **bcftools** is the standard command-line toolkit for filtering, subsetting and indexing VCFs. **SLURM** is the job scheduler used on academic compute clusters; a SLURM "array job" launches many independent copies of the same script with different indices. A variant is **polymorphic** in a panel when both alleles actually occur in that panel; a site where everyone is identical carries no information and is dropped. **kb** = 1,000 base pairs, **Mb** = 1,000,000. **30x** means each base was covered by about 30 sequencing reads on average — deeper, and therefore more accurate, than the roughly 4x of the older Phase 3 release.

There is also a single-file standalone, `tools/panmixer/obfuscate.py`. Its docstring reads: *"Runs the full pipeline: PMI/utility computation -> optimize -> stack -> VCF output."* It runs the **full per-subject pipeline in one process** for one subject and one chromosome, and it **requires `blocks_dict.json` and `allele_frequencies.npy` as pre-existing inputs** — it consumes the preprocessing products rather than producing them, and it replaces the experiment-directory machinery with command-line arguments. **Measured in this project**: 68 s wall-clock and 1.7 GB peak resident memory for one subject on chr21 (`docs/memory.md`, 2026-09-18; log `logs/pm_obf.log`). For calibration, the **published** Table 1 reports **5 min / 14 GB** for the longest-running chromosome on HPRC v1 with 44 subjects, all 22 chromosomes run in parallel. Those are not the same measurement — ours is chr21, theirs is the longest chromosome — but scaling chr21 to chr2 by record count gives roughly 6.2 min / 9.3 GB, which is the same order.

### Step 1 — Download the pangenome VCF
**File:** `starting_data/scripts/get_pangenomes.sh`
**Input:** a URL on the Human Pangenomics AWS Open Data S3 bucket.
**What happens:** `wget` fetches `hprc-v1.0-pggb.grch38.1-22+X.vcf.gz` and its index. This is the HPRC v1.0 PGGB graph deconstructed to a VCF with `vg deconstruct -a` against a GRCh38 backbone. `vg deconstruct` walks the graph's bubbles and writes each out as a VCF record; the `-a` flag retains **nested** variants — a variant occurring inside the span of a larger one, e.g. a SNP inside a long insertion — rather than collapsing them.
**Output:** `starting_data/pangenome.vcf.gz` (+ `.tbi` index).
**Expect:** contig names of the form `grch38#chr1` … `grch38#chrX` (the graph's PanSN path-naming convention, not bare `chr1`), and **45 samples** in the `#CHROM` header line — 44 HPRC individuals plus `chm13`.
**Why it matters:** this file *is* the thing being protected. Every downstream array is indexed by its row order.

### Step 2 — Download the PanGenie callset
**File:** `starting_data/scripts/get_pangenie_alignments.sh`
**Input:** Zenodo record 7669083, file `grch38_all-samples_bi_all.vcf.gz`.
**Output:** `starting_data/PG.vcf.gz` (+ `.tbi`), saved under the short name `PG`.
**Expect:** GRCh38 coordinates, bare `chr21`-style contig names, and **368 samples** (**measured in this project**, recorded in `docs/data.md`'s row for `external/PanMixer/starting_data/PG.vcf.gz`); chr21 has **384,720 records**.
**Why it matters:** the paper uses PanGenie frequencies for structural and nested variants that a SNP panel does not contain. In the released code this file also, unexpectedly, defines the HMM's anchor set — see defect (d).

### Step 3 — Download the phased reference panel
**File:** `starting_data/scripts/get_1000g_phased.sh`
**Input, as shipped:** `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{N}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`.
**Output:** `starting_data/chr{N}/1000g_phased.vcf.gz`.
**Expect, as shipped:** 1000 Genomes **Phase 3**, 2,504 samples, phased, **GRCh37 coordinates**. This is the build-mismatch problem — see defect (a). **In our rebuilt chr21 this filename slot instead holds the 30x GRCh38 panel: 1,002,753 records, 3,202 samples**, array shape `(3202, 1002753, 2)`. The original GRCh37 file is retained beside it, renamed, so the deviation is visible and reversible (`docs/data.md`).
**Why it matters:** this panel supplies the reference haplotypes used for LD block boundaries, for common-SNP allele frequencies, and as the attack database in the evaluation. Three load-bearing roles, one file.

### Step 4 — Remove chrX
**File:** `scripts/remove_X.sbatch`
**Input:** `pangenome.vcf.gz`.
**What happens:** `bcftools view --regions grch38#chr1,…,grch38#chr22`.
**Output:** `pangenome_no_X.vcf.gz`.
**Expect:** identical sample set, autosomes only.
**Why it matters:** humans have 23 chromosome **pairs**, i.e. **24 distinct chromosome sequences** — 22 autosomes plus X and Y. (For scale: a **haploid** human genome is about 3.1 billion base pairs, and a diploid cell contains about twice that, since it carries two copies of each autosome.) X has sex-dependent ploidy — one copy in males, two in females — which breaks the uniform "every sample has exactly two haplotypes" assumption that every later array depends on. Dropping it is the cheapest way to keep that invariant.

### Step 5 — Remove the chm13 sample
**File:** `scripts/remove_chm13.sbatch`
**Input:** `pangenome_no_X.vcf.gz`.
**What happens:** `bcftools view -s ^chm13`.
**Output:** `pangenome_no_X_no_chm13.vcf.gz`.
**Expect:** **44 samples = 88 haplotypes.**
**Why it matters:** **CHM13 is human.** It is a human cell line derived from a **complete hydatidiform mole** — an abnormal human conception in which all chromosomes come from the father, so the resulting cell line is effectively homozygous throughout. That homozygosity is precisely what makes it an excellent substrate for building a reference assembly, and CHM13 is used here as a **reference assembly**, not as a cohort donor whose privacy is at stake. It should not be described as non-human. Leaving it in the cohort would pollute the population model with a sample that is not a typical individual and is not a person to protect. The published paper reports the cohort as **44 individuals**, matching this count exactly. (The superseded preprint said 47; that discrepancy is preprint-only.)

### Step 6 — Split per chromosome
**File:** `scripts/split_data.sbatch`
**Input:** the 44-sample pangenome and `PG.vcf.gz`.
**Output:** `chr{N}/pangenome.vcf.gz` and `chr{N}/PG.vcf.gz`, each indexed.
**Expect:** chr21 pangenome has **340,824 records** (**measured in this project**, `docs/data.md`).
**Why it matters:** every later stage is embarrassingly parallel across chromosomes, and the paper notes the optimization is done per chromosome for efficiency. Note the different region syntax for the two files (`grch38#chr21` versus `chr21`) — a standing reminder that these datasets do not share a contig-naming convention.

### Step 7 — Count possible alleles per site
**File:** `scripts/src/get_possible_alleles.py`
**Input:** `chr{N}/pangenome.vcf.gz`.
**What happens:** for each record, `len(ALT.split(',')) + 1` — the number of comma-separated alternates, plus one for the reference.
**Output:** `chr{N}/num_alleles.npy`, shape `(340824,)`, dtype int64.
**Expect:** every entry ≥ 2; many entries > 2, because the graph is deeply multi-allelic.
**Why it matters:** it bounds the allele-index space and is used as an assertion in the allele-frequency step — if frequency mass ever lands on an allele index beyond what the site declares, two datasets have been mis-joined and the assertion fires immediately rather than silently corrupting the model.

### Step 8 — Compute LD blocks with plink
**File:** `scripts/get_blocks.sbatch`
**Input:** `chr{N}/1000g_phased.vcf.gz`.
**What happens:** `plink --vcf … --make-bed` converts to plink's binary genotype format (`.bed` genotypes, `.bim` variants, `.fam` samples), then `plink --blocks no-pheno-req` runs plink's confidence-interval block-detection procedure and writes the block intervals.
**Output:** `chr{N}/blocks.blocks` and `chr{N}/blocks.blocks.det` (a whitespace table with `BP1`, `BP2` start/end base-pair coordinates per block).
**Expect:** **14,137 blocks** on our rebuilt chr21 from the **30x GRCh38 panel** (**measured in this project**, `docs/memory.md`, 2026-09-18).
**Why it matters:** this defines the unit of editing. Note that the boundaries are derived from the *reference panel*, not from the pangenome cohort — 44 people is far too few to estimate LD structure reliably, so using an external panel of thousands is the right call. It also means the panel's genome build determines where the boundaries land, which is why defect (a) reaches into this step.

### Step 9 — Convert VCFs to numpy
**File:** `tools/common/VCFtoNP.py`
**Input:** the three per-chromosome VCFs. `pangenome` and `1000g_phased` are read as phased (`|`), `PG` as unphased (`/`).
**What happens:** two passes — count rows and samples, then fill the array.
**Output, and this is the layout everything else assumes:**

> `X[num_samples, num_sites, 2]`, dtype `int16`.
> `X[i, s, h]` = the **allele index** carried by sample *i*, at site *s*, on haplotype *h* ∈ {0,1}.
> `0` = the reference allele. `1, 2, 3, …` = the first, second, third… ALT allele **as listed in that VCF row**.
> **`-1` = missing** — either the VCF wrote `.` for that allele, or the genotype was haploid so there was no second entry.

Plus `*_positions.npy`, shape `(num_sites,)` int64, holding the base-pair POS column; and `*_subjects.npy`, the sample names in column order.
**Expect, chr21 (measured in this project):** `pangenome.npy` `(44, 340824, 2)` int16; `PG.npy` `(368, 384720, 2)`; `1000g_phased.npy` `(3202, 1002753, 2)` in our 30x GRCh38 rebuild.
**Why it matters:** VCF parsing is slow and everything downstream is array arithmetic. **The `-1` is load-bearing**: it is *not* the reference allele, and treating it as 0 would silently invent data — turning "we do not know" into "this person matches the reference", which is exactly the kind of error that inflates both frequency estimates and apparent matches. Every later step masks on `!= -1` before counting.

### Step 10 — Compute variant mappings between datasets
**File:** `scripts/src/get_mappings.py`
**Input:** the three VCFs and their `_positions.npy`.
**What happens:** builds, for each VCF, a dictionary keyed on the triple **(POS, REF, full ALT string)** → row index. Then intersects: for every pangenome key present in a target VCF, record pangenome-row → target-row. It additionally builds a **relaxed** mapping keyed only on (POS, REF), taking the first match.
**Output:** `pangenome_to_thousand_g_phased.pickle` (to the phased panel), `pangenome_to_thousand_g_alignments.pickle` (to PanGenie — note the name; see defect (d)), the paired index arrays `pangenome_mask.npy` / `thousand_g_phased_mask.npy`, and the `_posref` relaxed variants.
**Expect (measured in this project, `docs/memory.md` 2026-09-18 and `docs/data.md`):** on our rebuilt chr21 against the **same-build 30x GRCh38 panel**, **258,610** strict matches and **268,851** relaxed — 75.88% of the 340,824 pangenome records. Against the panel the released script actually downloads (**Phase 3, GRCh37**), the identical computation yields **1,025** strict and **2,176** relaxed — 0.30%. A **252× gap**. That single number is how you diagnose defect (a), and an independent `bcftools` comparison predicted the same figure, so two routes agree.
**Why it matters:** these mappings are the join keys for everything downstream — which pangenome sites have trustworthy external frequencies, which sites are scorable by the attack, and which sites count as anchors.
**Pedagogical note:** the strict key requires exact match on the *entire* ALT string. In a deeply multi-allelic graph VCF, two callsets that discovered the same SNP but a different set of neighbouring alternates produce different ALT strings and will **not** match. That is why the relaxed (POS, REF) key exists as a fallback, and it is also why the residual ~24% of non-matching records under the correct build should not be assumed to be a problem: the pangenome carries graph-only and structural variants that a SNV/indel/SV panel simply does not have. That residual has not been characterised here and should not be assumed benign either — it is an open question.

### Step 11 — Build the biallelic SNP mask
**File:** `scripts/src/build_biallelic_snp_mask.py`
**Input:** `chr{N}/pangenome.vcf.gz`, `pangenome_positions.npy`.
**What happens:** marks `True` where the record has exactly one ALT and both REF and ALT are a single base in {A, C, G, T}. It asserts that record *i*'s POS equals `positions[i]`, so any row-order drift is caught loudly rather than silently.
**Output:** `chr{N}/biallelic_snp_mask.npy`, shape `(340824,)`, dtype bool.
**Expect, chr21 (measured in this project):** **274,458 biallelic SNPs, 80.5%** of records.
**Why it matters:** the linkage attack's scoring model is a biallelic **Hardy–Weinberg** genotype model, which only makes sense on biallelic sites.

> **Box — Hardy–Weinberg, since the entire empirical privacy claim rests on it.** Let **p** be the reference-allele frequency at a biallelic site and **q = 1 − p** the alternate-allele frequency. If a person's two chromosome copies were drawn independently from the population, then
> **P(0/0) = p²,  P(0/1) = 2pq,  P(1/1) = q².**
> (The factor 2 on the heterozygote is because there are two ways to get one of each.) This is **Hardy–Weinberg equilibrium**: the assumption that alleles pair at random. It is what turns an *allele* frequency, which is what the panel gives you, into a *genotype* frequency, which is what the attacker's score needs. It fails in real populations under inbreeding, population structure and selection, so it is an approximation — but a standard and generally serviceable one. Reminder: **this p and q are allele frequencies**, not the stay/switch probabilities of §A.3 and not the input paths of Section 4.

### Step 12 — Refine plink blocks into the block dictionary
**File:** `scripts/src/get_blocks_simple.py`
**Input:** `blocks.blocks.det` and `pangenome_positions.npy`.
**What happens:** reads the (BP1, BP2) intervals, drops degenerate ones where start == end, sorts by start, and for every pangenome variant position uses a binary search (`np.searchsorted`) to find the interval with the greatest start ≤ position; the variant is *in* that block if it is also ≤ that interval's end. Variants inside no interval each receive **their own unique block id**.
**Output:** `chr{N}/simple_blocks_idx.npy` (shape `(340824,)`, one block id per variant) and `chr{N}/blocks_dict.json`, a JSON map `block_id → (list_of_positions, list_of_row_indices)`.
**Expect, chr21 (measured in this project, `docs/memory.md` 2026-09-18):** **100,757 dictionary entries** from only 14,137 plink blocks. Of those entries, **89,087 — 88.4% of entries — hold exactly one variant**, and since each holds one variant they account for 89,087 of 340,824 variants, i.e. **26.1% of variants**. Say which denominator, every time. (Most of these singleton entries arise because the variant fell outside every plink interval and was given a block of its own; an interval that happens to contain exactly one pangenome variant would also count as one. We have not separated the two sources, and the distinction does not affect anything downstream.) Separately, **10.4% of entries contain ≥2 anchors** and therefore take the HMM path, and those entries hold **71.4% of all variants**.
**Why it matters:** this is the index structure the whole mechanism iterates over. Internalise its shape: *most entries are trivial, but most of the genome lives in the non-trivial minority.* Note also the ordering convention — `block[0]` is **positions** (base pairs) and `block[1]` is **row indices**. Confusing these two is the root of defect (b).

### Step 13 — Compute allele frequencies
**File:** `scripts/src/get_af.py`
**Input:** `pangenome.npy`, `1000g_phased.npy`, `PG.npy`, both mapping pickles, `num_alleles.npy`.
**What happens:** for each pangenome site *i*, a three-way priority, documented in the file's own docstring as "the 1000g phased dataset is the most reliable, the alignments dataset second, the pangenome least":
1. if *i* maps into the phased panel, count alleles **from the panel only**;
2. else if *i* maps into PanGenie, count from PanGenie **and** from the pangenome;
3. else (a **novel** variant, in neither), count from the pangenome only, then pad the reference-allele count so the row sums to `2 × (N_panel + N_pangenome)`.

Then divide each row by its own sum.
**Output:** `chr{N}/allele_frequencies.npy`, shape `(340824, max_allele + 1)` float64 — `(340824, 90)` in our build, i.e. up to 90 distinct alleles at a single site.
**Expect:** every row sums to 1; the assertion that no probability mass sits beyond `num_alleles[i]` holds for every site.
**Why it matters:** this array *is* the population model. Case 3 deserves attention: a variant that nobody outside the graph has ever reported is given a frequency with a *large* denominator — on our 30x GRCh38 rebuild, 2 × (3,202 + 44) = **6,492**; on the shipped Phase 3 panel it would be 2 × (2,504 + 44) = 5,096 — so it comes out very rare, so its −log is very large, so the optimizer prioritises it. That is deliberate and conservative: it is the mechanism by which novel variants and under-represented populations get protected first.
**Caution for our own reuse:** in cases 2 and 3 the counts **include the target**. A model that includes the person it is scoring is not a population model in the sense §A.4 needs (see defect (c) for the same error in a more damaging place).

### Step 14 — Precompute per-block privacy, utility, and the sampled replacement
**Files:** `scripts/src/get_support_and_pmi.py` with `scripts/src/hmm.py`. This is the heart of the tool. Run as a SLURM array of 968 tasks = 44 subjects × 22 chromosomes.
**Input:** `pangenome.npy`, `allele_frequencies.npy`, `blocks_dict.json`, `pangenome_to_thousand_g_alignments.pickle`, and a subject id.

It does four things.

**(14a) The utility vector.** In code, `not_missing = np.sum(pangenome != -1, axis=(0,2))`. In words: **for each site, count how many of the 88 haplotype slots have a called allele** — the sum runs over the sample axis and the haplotype axis, leaving one number per site. Then `utility_loss[s] = 1/not_missing[s]` where that is positive, else 0. Output `utility_loss.npy`, shape `(340824,)`. This is the `1/support(v)` of §A.5. Note precisely what is counted: **haplotypes with any called allele**, not haplotypes carrying a particular allele.

**(14b) Per-block privacy (`pmi_blocks`).** First a per-site self-information. The code one-hot encodes the target's allele and multiplies by `−log(AF)`; in words, **for each site, pick out the −log frequency of the specific allele the target carries at that site**. Missing sites get 0. Then per block:
- block has exactly 1 variant → sum the per-site values;
- block has ≥2 variants but **0 anchors** → sum the per-site values;
- otherwise → `−forward_algorithm(haplotype, block)` for each of the two haplotypes.

Output `subjects/{id}/pmi_blocks.npy`, shape `(num_blocks, 2)` float64. Column *h* is the privacy value of obfuscating haplotype *h* of that block.

**Worth noticing while reading:** the forward algorithm is used whenever a multi-variant block has **≥1** anchor, but the *sampler* below switches to the HMM only at **≥2** anchors. So a block with exactly one anchor is *scored* by the HMM and *resampled* from allele frequencies. Nothing depends on this being wrong, but the asymmetry is real and it is easy to misread the code otherwise.

**A second thing worth noticing**, because it bites the rarest variants: an allele whose frequency comes out exactly 0 gives `−log 0 = +inf`. The code detects values above 10²⁰⁰ and sets them to **0** — the minimum — rather than the maximum. The pipeline copy carries a literal `#FIX ME` comment at that line. So the most identifying alleles of all can be scored as worth nothing. This is recorded in `docs/bugs.md` as a related finding under the ε_pmi entry.

**(14c) Sample a replacement for every block, both haplotypes.** Start from an array filled with −1.
- Singleton blocks, or blocks with ≤1 anchor: draw each variant's new allele independently from that site's normalised allele-frequency row.
- Otherwise: `sample_block_prior` walks the Li–Stephens chain over the block's anchors, picking a start donor uniformly among the **86** non-target haplotypes and then stepping donor to donor; anchor sites take the chosen donor's allele; non-anchor sites are meant to inherit from the donor chosen at the nearest anchor.

Output `subjects/{id}/new_haplotypes.npy`, shape `(340824, 2)`. **Only one candidate is generated per block per haplotype** — this is not a menu the optimizer chooses from, it is a single take-it-or-leave-it offer. The script raises `RuntimeError("OBSFUCATION FAILED")` if the candidate is byte-identical to the original everywhere.

Two consequences of "one unconditional draw", both of which matter later. First, the draw is made **without comparing to the original** — there is no rejection loop forcing the candidate to differ, despite the paper describing such a constraint. Second, therefore, a move can be selected, charged a privacy credit ε_j > 0, and change nothing at all, at zero utility cost. That is a free lunch the LP will happily take.

**(14d) Per-block utility cost (`support_blocks`).** For a multi-variant block: sum `utility_loss[s]` over the sites where the candidate **actually differs** from the original, separately per haplotype. For a single-variant block it assigns that site's `utility_loss` to both columns unconditionally — i.e. a singleton's cost is charged whether or not the resampled allele happened to equal the original.
Output `subjects/{id}/support_blocks.npy`, shape `(num_blocks, 2)` float64.

**Expect if all is well:** `pmi_blocks` and `support_blocks` have identical shapes and are elementwise non-negative (the optimizer asserts both). `new_haplotypes.npy` differs from the original at a healthy fraction of called sites.
**Why it matters:** after this step the biology is gone. What remains is two arrays of numbers — a value vector and a weight vector — and a pre-computed answer for every move. Everything downstream is arithmetic.

### Step 15 — Total utility normaliser
**File:** `scripts/src/get_total_utility_loss.py`
**Input:** every subject's `utility_loss.npy`, across all 22 chromosomes.
**What happens:** per subject, sums `2 × Σ utility_loss` — the factor 2 because both haplotypes can be edited — and writes one number per subject.
**Output:** `starting_data/utility_loss.json`.
**Why it matters, and a point of genuine confusion the code forces us to state plainly.** Step 15 produces a *genome-wide* per-subject total, but **Step 17's optimizer recomputes its own denominator from the single chromosome it is solving**: `max_utility_loss = 2 × Σ utility_loss` over that chromosome's sites only. So the capacity budget as executed is **per chromosome**, and the genome-wide JSON is a reporting aggregate rather than the constraint the LP sees. Concretely, our chr21 denominator is 9,288.85 — a chr21 number, not a genome number.

**Name discipline, because "maximum utility loss" means two different things.** Call Step 15/17's `2 × Σ utility_loss` the **budget denominator**. Call the paper's *"Removed"* baseline — where the target individual is excluded from the graph entirely, which the published text describes as "the maximum possible utility loss for a given individual" — the **removal baseline**. They are not the same quantity, and the budget denominator is not even achievable: `support_blocks` only charges sites that actually *change*, so total achievable loss is strictly less than the denominator. **This is exactly why the knapsack in Step 17 saturates below its nominal budget**, and that causal link is the thing to carry, not the numbers.

### Step 16 — `experiment_starter`
**File:** `tools/common/experiment_starter.py`
**Input:** a capacities file (one float per line) and a subjects file (one sample name per line).
**What happens:** finds the next free `experiments/exp_N/`, takes the cross product subjects × capacities into a row table, attaches demographic columns, and pre-creates `data/chr{1..22}/{row_id}/` plus `data/all/{row_id}/`.
**Output:** `experiments/exp_N/data.csv` and the directory skeleton.
**Expect:** one row per (subject, capacity) pair; baseline modes instead write a single sentinel capacity (0.0 unedited, 1.0 empty, −1.0 remove-unique).
**Why it matters:** `row_id` is the primary key for the entire run — every later stage is addressed by `(chromosome, row_id)`.

### Step 17 — The optimizer
**File:** `tools/panmixer/optimizer.py`
**Input:** for one (chromosome, row_id): `support_blocks.npy` (weights), `pmi_blocks.npy` (values), `utility_loss.npy`, and the row's capacity.
**What happens:** computes `max_utility_loss = 2 × Σ utility_loss` for that chromosome and `target = capacity × max_utility_loss`. Flattens the two `(num_blocks, 2)` arrays to length 2k, builds a **GLOP** model — GLOP is Google's linear-programming solver, reached through the OR-Tools library — with one **continuous** variable per move in [0,1], objective Σ value·x maximised, one constraint Σ weight·x ≤ target. Solves, exits with an error unless the solver reports `OPTIMAL`, rounds the solution to integers, reshapes.
**Output:** `xsol.npy`, shape `(num_blocks, 2)` int, and `optimizer_stats.json`.
**Expect (measured in this project on HG00438 / chr21, 30x GRCh38 rebuild; `docs/memory.md` 2026-09-18, logs `logs/pm_obf.log`):**
- capacity 0.1 → **47,070 moves**, utility loss **928.92 / 9,288.85 = 10.000%** — the constraint is exactly tight, as it must be when the budget binds;
- capacity 0.5 → **177,865 moves**, utility loss **33.1%** — the knapsack **saturates below its budget**, because once every move with positive value has been taken there is nothing left to buy. This is the Step 15 point made concrete: the budget denominator is not reachable.
- Same seed → byte-identical `new_haplotypes.npy` and `xsol.npy`; different seed → different output; **zero** allele-bound violations.
**Why it matters:** this is the only place the privacy/utility trade-off is actually decided. Everything before it prices the options; everything after it executes the decision.

### Step 18 — The stacker
**File:** `tools/panmixer/stacker.py`, strategy `to_best`.
**Input:** `xsol.npy`, `new_haplotypes.npy` (the precomputed candidates), `pangenome.npy`, `blocks_dict.json`.
**What happens:** copies the target's original haplotypes, then for every `(block index i, haplotype j)` with `xsol[i,j] == 1`, writes the candidate's alleles into every site of block *i* on haplotype *j*. **Blocks not selected are left exactly as they were.**
**Output:** `exp_N/data/chr{N}/{row_id}/new_haplotypes.npy`, shape `(340824, 2)`, plus `stacked.json`.
**Expect (measured in this project, HG00438 / chr21 / capacity 0.1):** **68,335 of 650,541** called alleles changed.
**Why it matters:** this is where PanMixer's guarantee structure is fixed for good. The "leave unselected blocks verbatim" behaviour is what makes the release cheap and high-utility, and it is also exactly what makes a worst-case pairwise guarantee impossible (§A.7). Also note the alternative strategies: `to_empty` blanks the subject entirely (the strongest baseline and the paper's "Removed" arm), `to_unedited` passes the original through (the null baseline), `to_random` an unguided draw.

### Step 19 — Write the VCF
**File:** `tools/common/convert_2_vcf.py`
**Input:** the stacked `new_haplotypes.npy` and the original `chr{N}/pangenome.vcf.gz` as a template.
**What happens:** streams the template line by line; on the header it locates the target's column index; on each data line it replaces that one field with `g0|g1`, writing `.` for −1. Any other sample field lacking a `|` is rewritten to `.|.`. For the `to_empty` baseline it drops the column entirely. It asserts that the number of data lines equals the number of genotype rows. Then `bgzip` (block-gzip, which keeps the file randomly seekable) and `bcftools index`.
**Output:** `new_haplotypes.vcf.gz` (+ index).
**Expect:** same site count, same sample list, one column's contents different.
**Why it matters:** **this is the released artefact.** Note exactly what it means: the release is the *whole cohort VCF*, with one person's column rewritten in place. The graph topology, the sites, and the other 43 people are untouched. This is the precise sense in which PanMixer "modifies the pangenome" — it does not modify the graph.

### Step 20 — Combine chromosomes
**File:** `tools/common/combine_vcfs.py` — `bcftools concat` of chr1…chr22 into `data/all/{row_id}/merged.vcf.gz`.
**Why it matters:** the re-identification attack and the read-mapping evaluation operate genome-wide or need a single index.

### Step 21 — Downstream: the linkage (gap-score) attack
**Files:** `tools/downstream/privacy/gap_score.py` (single) and `diploid_gap_score.py` (all variants of the attack).
**Input:** the released `new_haplotypes.npy`, the masked attack database (`1000g_*_masked_posref.npy`), the biallelic SNP mask, allele frequencies.
**What happens:** restrict to sites that are biallelic SNPs *and* **polymorphic** in the attack panel (both alleles actually present — a site where every panel member is identical carries no information); binarise alleles to 0/1; compute each genotype as the sum of the two haplotypes (0, 1 or 2); model each genotype's probability as p², 2pq or q² under Hardy–Weinberg (see the box at Step 11); score a database entry by summing **−log(probability)** over the loci where its genotype *equals* the released genotype. Then compute self-score minus the maximum over the database. Four scoring variants are run: genotype, haplotype, both-haplotypes-equal, and unweighted match count.
**Output:** `genotypes_scores.npy` (one score per database entry), `genotypes_scores_self.npy`, and `gap_score_all.json` with `genotype_g_to_gstar`.
**Expect:** a negative `g_to_gstar` means the attack failed for that subject at that capacity. Sweeping capacity and finding where every subject goes negative is how ε_private is produced.
**Why it matters:** this is the actual empirical privacy measurement. The whole published privacy claim rests on this one function's behaviour, which is why defect (e) — the missing self-exclusion — is as consequential as it is. Note also the sign caution repeated from §A.7: the published equation (5) as typeset and the "emphasizes rare matches" reading differ by a sign, and that is unresolved.

### Step 22 — Downstream: membership inference
**File:** `tools/downstream/privacy/MIA_privacy.py`
**Input:** the released haplotypes, the original, the other 43 pangenome samples, and optionally the external panel.
**What happens:** counts matching haplotype alleles and matching genotypes between the release and each comparison set, at all sites and at panel-overlapping sites.
**Output:** `MIA_privacy.npz` per (chromosome, row_id).
**Why it matters:** **membership inference** is a weaker but often easier attack than re-identification. The adversary already has the target's full genome and only wants to know *whether the target is in the cohort at all*. That alone can be sensitive — if cohort membership implies a diagnosis, then confirming membership discloses the diagnosis. This analysis appears in the **published** version of the paper and not in the superseded preprint.

### Step 23 — Downstream: allele-frequency loss
**File:** `tools/downstream/utility_in/af_loss.py`
**What happens:** recomputes allele frequencies over the cohort before and after, takes the L1 difference per site, and reports it overall and stratified by **minor allele frequency (MAF)** bands (<0.01, 0.01–0.05, 0.05–0.1, 0.1–0.5) and by SNP versus complex. MAF is the frequency of the *less* common allele at a site, so it always lies in [0, 0.5] and is the conventional axis for "how rare is this variant".
**The metric.** The paper reports a **Wasserstein divergence** between the before and after allele distributions. For alleles modelled as Bernoulli random variables, the published Methods state that this equals **two times the mean absolute difference** between allele frequencies, and they give the reason: *"The factor of two arises because each diploid individual carries two alleles per variant."* If you took the earlier, general "earth-mover's distance" description literally you would compute |p − q| and get half this; the factor of two is the Bernoulli-specific specialisation, and the two descriptions are consistent once you know that.
**Expect (published, at ε_private = 0.002), with labels attached:** Wasserstein divergence **0.006 for all alleles**, **0.006 for SNPs only**, **0.005 for rare SNPs (population MAF < 0.05)**. The paper's **"Removed" baseline** — the target excluded from the graph entirely — is **0.027**, close to 5× larger, and that contrast is the calibration that makes the 0.006 meaningful. (The superseded preprint's values were 0.004 / 0.004 / 0.002; do not quote them.)
**Why it matters:** the allele-frequency spectrum — the histogram of how many variants sit in each frequency band — is the single most-used summary in population genetics, and it is the main input to an **association study**, which tests variant by variant whether carrying one allele correlates with a trait or disease across a large sample. If PanMixer distorted it badly, the released graph would be useless for exactly the science it exists to serve.

### Step 24 — Downstream: LD loss
**File:** `tools/downstream/utility_in/ld_loss.py`
**What happens:** for each SNP, gathers the SNPs within a **5 kb** window, computes **D = P(A and B) − P(A)·P(B)** — the covariance of allele co-occurrence, zero exactly when the two sites are independent — and then

> **r² = D² / [p(1−p)·q(1−q)]**

which is D rescaled to sit between 0 and 1: **the squared correlation between the two loci treated as 0/1 variables.** r² = 1 means knowing one allele tells you the other exactly; r² = 0 means it tells you nothing. Computed before and after obfuscation, with the L1 difference accumulated. (Here p and q are the two loci's allele frequencies; and again, this r² is unrelated to the recombination constant r = 1.26.)
**Expect (published, at ε_private):** the LD matrix was on average **0.003** off, against **0.017** for the complete-removal baseline.
**Why it matters:** this is the direct check that the block-level design worked. If PanMixer had edited variants individually rather than in whole LD blocks, r² would be destroyed, and this metric is precisely where you would see it.

### Step 25 — Downstream: Beagle reconstruction attack
**Files:** `tools/beagle/beagle_refinement.py`, `accuracy_stats.py`
**What happens:** runs **Beagle** — a standard **imputation** tool, which fills in or corrects genotypes by matching a target against a reference panel of haplotypes — on the released genotypes, using the external panel as reference with the targets excluded. Then compares the refined output to the *original* genotypes and to a naive "most frequent allele" baseline. The repository pins the jar **`beagle.27Feb25.75f.jar`** (downloaded by `external_tools/get_tools.sh`), i.e. the **27Feb25 build**; that is the version this project can verify. Note that the repo's `environment.yaml` ships no Java, so the jar has no declared runtime.
**Why it matters:** it tests whether obfuscation is *reversible*. If an off-the-shelf imputation tool can restore the original from the release, the privacy was cosmetic. The mosaic construction is what is supposed to prevent that: an imputer will happily snap the released haplotype onto *some* plausible reference haplotype, but the right one should no longer be distinguished. It is worth seeing that imputation is the **same Li–Stephens machinery as §A.3 run in the other direction** — there, sampling a plausible haplotype from a panel; here, inferring the most plausible haplotype given partial observations. That is why the attack is natural and why it is the right test of this particular mechanism.

### Step 26 — Downstream: read mapping with vg giraffe
**Files:** `tools/downstream/utility_out/vg_prep.py` and `quick_align.py`
**What happens:** `vg_prep` normalises the released chr21 VCF (`bcftools view -c 1` to drop sites now carried by nobody, `norm -N -f REF`, `--trim-alt-alleles`, rename contigs, restrict to SNPs, sort), then `vg autoindex --workflow giraffe` builds the graph indexes. `quick_align` runs `vg giraffe` on chr21 **FASTQ** reads — FASTQ is the standard text format holding raw sequencing reads plus a per-base quality score — from five 1000 Genomes individuals **not in the graph**, then runs `vg stats -a` on the resulting **GAM** (Graph Alignment/Map) file.
**Output:** per-sample `*.giraffe.stats.txt`, parsed for percent of reads aligned perfectly, gapless (no insertions or deletions in the alignment), and at **MAPQ 60** — mapping quality 60 being the aligner's maximum-confidence value, meaning it considers the placement essentially unambiguous.
**Expect (published), with labels attached:** obfuscated graph **77.82% perfect / 95.61% gapless / 77.01% MAPQ 60**, against the baseline graph's **77.83% / 95.62% / 77.01%**. The differences are in the second decimal place.
**Why it matters:** this is the *external* utility test, and it is the one that matters most for the pangenome's reason to exist. A pangenome exists to reduce **reference bias** — the tendency of an aligner to place reads better when they happen to resemble the single reference, which systematically disadvantages people genetically distant from whoever the reference was built from. If obfuscation degraded mapping for outside individuals, the graph would have lost the thing it was built for. Note that **this project has not run this evaluation**: the five read FASTQs and the reference FASTAs it needs are not downloadable from the repository and no script in the checkout produces them (`docs/data.md`).

### Step 27 — `gather_results`
**File:** `tools/common/gather_results.py` — walks every `(chromosome, row_id)` directory, collects each JSON, and aggregates according to the `AGGREGATION_DICTIONARY` in `constants.py`: privacy and utility summed across chromosomes, rate-like metrics averaged.
**Why it matters:** privacy and utility are budgeted per chromosome (Step 15) but the claim is genome-wide, so the sums have to be reassembled. Summing a budget that was applied chromosome by chromosome is only meaningful because the blocks are treated as independent; that independence is an assumption, not a measurement.

---

## Where the implementation departs from the paper

Five findings, each established by reading the cited code and by executing it. **These are defects in the released code relative to the published methods.** We have **not** reproduced the paper's figures, and nothing here asserts that the paper's published results are wrong — the repository ships an *unwired* GRCh38 script and contains code branches that prefer 3,202-sample artifacts, which is consistent with the authors having run a path the release does not wire in. That is an inference, not a measurement. Each finding is included because it teaches something real about how the method is supposed to work. All five are recorded with evidence in `docs/bugs.md`, with the supporting measurements dated 2026-09-18 in `docs/memory.md`.

### (a) The LD blocks and mappings are computed against the wrong genome build

`get_1000g_phased.sh` downloads the 1000 Genomes **Phase 3** release `20130502`, which is in **GRCh37** coordinates. The pangenome is **GRCh38**. Steps 8 and 10 then intersect them by base-pair coordinate.

**Measured in this project** on chr21: the strict (POS, REF, ALT) intersection gives **1,025** matches (0.30% of 340,824 records) against the GRCh37 Phase 3 panel, versus **258,610** (75.88%) against the same-build 30x GRCh38 panel — a **252× gap**. Two independent routes agree: PanMixer's own `get_mappings` and a direct `bcftools` comparison. The published Data availability section pins the same GRCh37 URL while the Methods state a GRCh38 backbone three times and the full published text contains zero mentions of GRCh37, liftover, 30x, high-coverage or NYGC — so the inconsistency is in the paper as well as the code, and the published version sharpens it rather than resolving it (the preprint had no Data availability section at all). There is no easy fix available upstream: the official GRCh38 liftover of Phase 3 was withdrawn in 2021, with 1000 Genomes redirecting users to the 30x GRCh38 collection — which is the panel we substituted, and the one the repo's own unwired `get_blocks_grch38.sh` downloads.

**What it teaches.** Genome coordinates are only meaningful relative to a named build. GRCh37 and GRCh38 differ by insertions, deletions and rearrangements, so the same base-pair number denotes different DNA. Joining two datasets on coordinate without checking the build is the genomics equivalent of joining two tables on an integer key that means different things in each — and, crucially, it fails *silently*, producing a small valid-looking join rather than an error. Concretely, it means block **boundaries** are placed using positions that do not correspond to the pangenome's sites, and the **attack database** and the **AF/LD utility site sets** are built on a near-empty intersection. A useful canary: print `len(pangenome_to_thousand_g_phased)` after `get_mappings`; near-zero confirms it.

### (b) The Li–Stephens sampler never recombines, so the mosaic is inert

Four things compound in `starting_data/scripts/src/hmm.py` (and identically in the standalone `tools/panmixer/obfuscate.py`):

1. `effective_N = 1.0 / 10_000.0` — the **reciprocal** of the paper's stated N_e = 10,000. The published Fig. 6 gives d = Δx · N_e · r, with N_e in the numerator.
2. `get_anchor_snps` returns pairs `(snp, block[1][j])` where `snp` iterates over `block[1]` — the list of **VCF row indices**, not `block[0]`, the base-pair positions. So *both* columns of the returned array are row indices, and `get_transitions(pos_1, pos_2)` differences **row indices**. The paper specifies Δx in **centiMorgans**.
3. The scorer (`get_transitions`) computes `d = Δ/1e6 · effective_N · r` while the sampler (`get_transition_matrix_without_subject`) computes `d = Δ/1e6 · 4 · effective_N · r` — the two paths differ by a **factor of 4** that appears nowhere in the paper. Since Fig. 6 has no factor 4, the **sampler** is the deviation.
4. The non-anchor fill builds `anchor_pos_to_state` keyed on those row indices, then iterates `block_positions` (real base pairs) and tests `if pos in anchor_pos_to_state`. Comparing a base-pair coordinate against a row index, the guard **never fires**, so every site in the block — anchors included — falls through to the `np.argmin(np.abs(anchor_pos_array − pos))` fallback. With bp values on one side and row indices on the other, that argmin always selects the **last** anchor.

The arithmetic consequence: with row-index differences of order 1 to 10³, d lands around 10⁻¹⁰ to 10⁻⁷, so `exp(−d/86) ≈ 1`, so the switch probability q is on the order of 10⁻¹² per step. **The chain essentially never leaves its starting donor** — and even if it did, defect 4 would overwrite the whole block from a single donor anyway.

**Measured in this project (`docs/bugs.md`, `docs/memory.md` 2026-09-18).** Executing the shipped `HaplotypeHMM` class directly: total switch mass **5.8 × 10⁻¹²** at a 1-index gap and 5.8 × 10⁻⁶ at a gap of 10⁶. With the paper's own Fig. 6 constants the switch mass would be **0.13–0.99**; that is a gap of ten to eleven orders of magnitude. And on **300 of 300 randomly chosen HMM-path blocks of real chr21 for subject HG00438**, the sampled block was **identical to a single donor haplotype copied verbatim**. (Blocks were chosen at random from the HMM-path set; the earlier synthetic-panel result reproduced on real data.)

**What it teaches.** Everything Li–Stephens buys you comes from the *switching*. The mosaic is what makes the output a combination no individual possesses. Without switching, the released block is one other real person's real haplotype. That is not nothing — it does break the link to the target — but it is a categorically different object from a mosaic, and it means the released haplotype now carries a *different* individual's data verbatim, which is its own disclosure question. It is also a good lesson in why genetic distance is measured in centiMorgans rather than base pairs or row indices: recombination probability is a function of *genetic* map distance, and a row index is not even monotone in physical distance in any metric sense. The consequence for our own work is direct: "reuse PanMixer's cohort HMM" is not available as a plan. The log-space forward recursion is correct and O(K) per step and can be kept; the transition constants must be re-derived from Fig. 6 with a real genetic map.

### (c) The privacy scorer includes the target in its own reference panel

`forward_algorithm` uses `self.pangenome_haplotypes` — the concatenation of all 88 haplotypes, **including the target's own** — and `self.num_states` = 88, while normalising its transitions by the without-target count of 86. The sampler, by contrast, correctly uses `pangenome_haplotypes_without_subject` (86 states).

The consequence is a **ceiling**, and it needs to be stated carefully because the naive version of the statement is refuted by our own measurement. If the target's own haplotype is one of the states, then the model can always explain the observation by copying the target: pick that state (probability at least 1/2N = 1/88 under the uniform start) and match at every anchor (probability 1 − μ_mis per anchor, where **μ_mis = 10⁻⁴ is the HMM's per-anchor mismatch probability** — a fourth, unrelated use of the letter ε in the source, which is why I have renamed it). So

> p(h_bj) ≥ (1/88) · (1 − μ_mis)^n, where n is the number of anchors in the block,

and therefore

> ε_j = −log p(h_bj) ≲ log(88) + n·μ_mis = **4.4773 nats plus a small additive term**.

In words: **the block can never look much more surprising than log 88, no matter how rare it truly is, because the model always has the option of explaining the block as a copy of the target.** Note that this is a soft bound, not a hard cap: the additive per-anchor term means a measured value can legitimately sit slightly above 4.4773.

**Measured in this project** on **150 real HMM-path blocks for HG00438 on chr21** (`docs/bugs.md`; `docs/memory.md` 2026-09-18, Correction 3), by executing the shipped `HaplotypeHMM` class on real chr21 data: shipped ε_j spans **0.071 to 4.508**, median **1.242**, with only 9 of 150 within 0.05 nats of 4.4773. The maximum of **4.508 exceeds log(88) = 4.4773 by 0.031 nats**, which is exactly what the additive term predicts — 0.031 / 10⁻⁴ ≈ 300 anchors, an unremarkable size for a large block. So the measurement is *consistent with* the argument, not a counterexample to it. Re-scored with the target **removed** from the emission panel, the median moves only slightly, to **1.299**, but the **maximum rises to 78.1**. For the 9 ceiling-bound blocks the leave-one-out values have median 12.06 and maximum 78.12 — an understatement of up to **17.4×**.

An earlier version of this finding claimed ε_j was "pinned near log(2N) for essentially every block". That was **wrong** and is retracted; the retraction is recorded in `docs/memory.md` as Correction 3. The effect is a *ceiling*, not a *flattening*: typical blocks are barely affected (median ratio 1.03×) while the rarest are truncated hardest.

**What it teaches.** This is the most instructive defect in the set, because it is a *leave-one-out* error, and leave-one-out errors are endemic to privacy evaluation. The question "how surprising is this person's haplotype under a population model?" is only meaningful if the population model was not fitted on that person. Include them and you are instead asking "how surprising is this haplotype given that I already know someone has it", whose answer is "not very". Worse, the distortion is not uniform: it selectively flattens exactly the tail of the distribution that the whole mechanism exists to find. A knapsack that cannot see that block A is 17× more identifying than block B will not prioritise block A.

**What this does not establish.** Nothing about whether the *paper's* published numbers are affected. ε_j enters the LP only through a ranking, and we have not measured how much the selection changes when the ceiling is removed. The claim to make is the narrow one: the shipped scorer truncates rare blocks.

### (d) The "anchor" set is PanGenie, not 1000 Genomes

`hmm.py` loads `pangenome_to_thousand_g_alignments.pickle` into a variable named `self.snp_positions`, with the comment `#all alignments are snps`, and uses its keys as the anchor set. But `get_mappings.py` binds `thousand_g_alignments` to **`PG.vcf.gz`** — the PanGenie callset. So the variable named after 1000 Genomes is bound to PanGenie, while the paper defines anchors as **V₁₀₀₀G** — "Established variants from the Phase 3 phased 1000 Genomes dataset, used as trusted anchor points." The actual 1000G map is a different file, `pangenome_to_thousand_g_phased.pickle`, which `hmm.py` never opens.

There is an important bound on the blast radius here, and it corrects an inference this project itself got wrong once. PanGenie **is** GRCh38, so the anchor mapping is a same-build join and is **not** corrupted by defect (a). Consequently the **per-block routing decision** — HMM path versus allele-frequency path — is *not* degraded by the build mismatch. What defect (a) *does* corrupt is block **boundaries** (plink runs on the 1000G panel), the **attack database**, and the **AF/LD utility site sets**. Two different notions of "anchor" are in play; do not conflate them.

**Measured in this project (recorded in `docs/plan.md`'s comparison table):** **273,475 anchor variants** on our rebuilt chr21, with **10.4% of block-dictionary entries** reaching the ≥2-anchor threshold. Note that this count is not the same object as the 258,610 strict 1000G matches (a different callset) or the 384,720 PanGenie chr21 records (records in the source file, before intersecting with the pangenome), so the three figures are not expected to agree and comparing them is a category error.

**What it teaches.** Two things. The mundane one: **a variable name is not a type**, and a comment asserting a fact is not the fact. This one cost this project a wrong inference before it was caught. The more interesting one: "anchor" is a design choice with real consequences. Anchors determine which blocks are rich enough to model with an HMM and which fall back to independent allele-frequency sampling — and since the HMM path is the only path that respects LD, the anchor set decides where LD is preserved. Whether PanGenie or 1000 Genomes is the *better* anchor set is a genuinely open question and this project has not answered it. The defect is that the code does not do what the paper says.

### (e) The gap-score attack does not exclude the target from its own maximum

The published equation (7) is explicit: GapScore = L(g_target, g\*_target) − **max_{i ≠ target}** L(g_i, g\*_target). The code at `tools/downstream/privacy/gap_score.py` computes `highest_score_genotypes = np.max(scores_genotypes)` over the **entire** attack database, with no exclusion. The code is aware this can matter — `diploid_gap_score.py` contains a debug branch that checks whether `subject_name` appears in the attack database's sample list and prints its rank, and `tools/common/utils.py` does delete the self row, but only when the name is present in the panel's sample list, so which panel is loaded determines whether it happens at all — yet the exclusion is never applied to the score written to `gap_score_all.json`.

**What it teaches.** Whether this changes anything depends on whether the target is present in the panel being used as the attack database — and HPRC individuals largely do appear in the 1000 Genomes 30x release. If the target *is* in the database, the maximum can be the target's own entry, which mechanically drives the gap score negative and makes the attack look defeated when it was not. Since ε_private is *defined* as the threshold where the gap score goes negative for everyone, this is the quantity on which the entire empirical privacy claim is calibrated. The general lesson is the same as (c) from the other direction: **in privacy evaluation, who is and is not in the comparison set is not a detail, it is the definition.** Anyone attempting to reproduce the published 0.002 must exclude the target row before taking the maximum, and must pin which attack database was used — a denser panel is a strictly stronger attacker, so the choice moves the privacy axis.

A related operational hazard, recorded in `docs/bugs.md`: `apply_mask.sbatch` writes the masked attack database but is **not referenced** by `run_get_starting_data_pipeline.sh`, and `diploid_gap_score.py` prefers 30x-sample attack-database files that **no script in the repository produces**, while `constants.py` hardcodes `NUM_1000G_SUBJECTS = 3202` although the shipped pipeline downloads 2,504. Following the README end to end therefore either fails to find the attack database or silently falls back to a different, weaker one. This is the dangerous kind of defect: it still yields numbers.

---

## The pipeline at a glance

| Step | Input | What happens | Output | Why it matters |
|---|---|---|---|---|
| 1–3. Download | URLs (S3, Zenodo, EBI FTP) | Fetch HPRC v1.0 PGGB VCF, PanGenie callset, phased reference panel (shipped: Phase 3 / GRCh37 / 2,504 samples; this project substituted 30x / GRCh38 / 3,202 samples) | `pangenome.vcf.gz` (45 samples), `PG.vcf.gz` (368 samples), `1000g_phased.vcf.gz` | Supplies the data to protect, plus two external models of "what is common in humans" |
| 4–6. Clean & split | `pangenome.vcf.gz` | Drop chrX (sex-dependent ploidy), drop `chm13` (a human cell line used as a reference assembly, not a cohort donor), split per chromosome | `chr{N}/pangenome.vcf.gz`, 44 samples = 88 haplotypes; chr21 = 340,824 records | Establishes the uniform "everyone has exactly 2 haplotypes" invariant everything else assumes |
| 7. Allele counts | per-chr pangenome VCF | Count `len(ALT)+1` per record | `num_alleles.npy` (340,824,) | Bounds the allele-index space; used as a consistency assertion later |
| 8. LD blocks | reference panel VCF | `plink --make-bed` then `plink --blocks` | `blocks.blocks.det` — 14,137 blocks on our chr21 / 30x GRCh38 | Defines the unit of editing: variants inherited together must be edited together |
| 9. VCF → numpy | 3 per-chr VCFs | Two-pass parse into `int16[samples, sites, 2]`; `-1` = missing | `*.npy`, `*_positions.npy`, `*_subjects.npy` | Turns text into arrays; the `-1` convention keeps "missing" distinct from "reference" |
| 10. Mappings | 3 VCFs | Join on (POS, REF, ALT) and relaxed (POS, REF) | mapping pickles + index masks — 258,610 strict on our GRCh38 rebuild, 1,025 with the shipped GRCh37 panel | The join keys for frequencies, anchors and the attack; where the build mismatch bites |
| 11. SNP mask | per-chr pangenome VCF | Flag records that are single-ALT single-base SNPs | `biallelic_snp_mask.npy` — 274,458 = 80.5% on our chr21 | The attack's Hardy–Weinberg model is only valid on biallelic sites |
| 12. Block dict | plink intervals + positions | Assign each variant to a block; variants in no interval get their own | `blocks_dict.json` — 100,757 entries, of which 89,087 (88.4% **of entries**) hold one variant = 26.1% **of variants** | The index everything iterates over; most entries are trivial but most variants are not |
| 13. Allele freqs | pangenome/PG/panel arrays + mappings | 3-way priority (panel → PanGenie + graph → graph with padded denominator 2×(N_panel+N_pangenome) = 6,492 on our rebuild) | `allele_frequencies.npy` (340,824 × 90), rows sum to 1 | The population model; novel variants are scored as very rare, so they are protected first |
| 14. Privacy + utility + candidates | everything above, one subject | Self-information per block (per-site −log AF, or HMM forward for ≥1 anchor); one sampled replacement per block per haplotype (HMM only at ≥2 anchors); cost = Σ 1/support over changed sites | `pmi_blocks.npy` (k,2), `support_blocks.npy` (k,2), `utility_loss.npy`, `new_haplotypes.npy` | Converts biology into a value vector and a weight vector — the whole trade-off in two arrays |
| 15. Total utility | all subjects' utility vectors | Sum 2 × Σ utility_loss per subject, genome-wide | `utility_loss.json` | A reporting aggregate; the LP's actual **budget denominator** is recomputed per chromosome in Step 17 |
| 16. Experiment start | subjects × capacities | Build row table and directory skeleton | `data.csv`, `exp_N/data/chr*/row_id/` | `row_id` is the key for the whole run |
| 17. Optimizer | `support_blocks` (weights), `pmi_blocks` (values), capacity | LP relaxation of 0–1 knapsack via OR-Tools GLOP, then round | `xsol.npy` (k,2) binary. HG00438/chr21: cap 0.1 → 47,070 moves at exactly 10.000% (928.92/9,288.85); cap 0.5 → 177,865 moves at 33.1% (saturates) | The only place the trade-off is actually decided |
| 18. Stacker | `xsol`, candidates, original | Overwrite selected blocks; **leave unselected blocks verbatim** | `new_haplotypes.npy` (340,824, 2); 68,335 of 650,541 called alleles changed at cap 0.1 | Executes the decision — and fixes the guarantee structure: verbatim blocks ⇒ no worst-case pairwise bound |
| 19–20. Write & merge | stacked array + template VCF | Rewrite one sample column; bgzip; index; `bcftools concat` across chromosomes | `new_haplotypes.vcf.gz`, `merged.vcf.gz` | **The released artefact**: the cohort VCF with one column rewritten in place; graph topology untouched |
| 21. Gap score | release + attack DB + AF | Hardy–Weinberg genotype probabilities; rarity-weighted matching; self-score minus best other | `gap_score_all.json`; negative = attack defeated | Produces ε_private = 0.002 (published) — the entire empirical privacy claim |
| 22. Membership inference | release, original, cohort, panel | Count matching alleles/genotypes against each set | `MIA_privacy.npz` | Tests the weaker "is this person in the cohort at all" attack; new in the published version |
| 23–24. AF / LD loss | before and after cohort arrays | Wasserstein (= 2 × mean abs AF difference) and L1 change in r² over 5 kb windows | per-row JSON | Published: AF 0.006 all alleles / 0.006 SNPs / 0.005 rare SNPs (MAF<0.05) vs 0.027 removal baseline; LD 0.003 vs 0.017 |
| 25. Beagle | release + reference panel (targets excluded) | Imputation-based reconstruction with the pinned `beagle.27Feb25.75f.jar`; compare to original and to a most-frequent-allele baseline | refined VCF + accuracy stats | Tests whether obfuscation is reversible by off-the-shelf tools |
| 26. vg giraffe | release chr21 → graph indexes; FASTQ reads from 5 non-cohort individuals | `vg autoindex --workflow giraffe`, `vg giraffe`, `vg stats -a` | perfect / gapless / MAPQ-60 rates — published 77.82 / 95.61 / 77.01 vs baseline 77.83 / 95.62 / 77.01 | The external utility test: does the graph still reduce reference bias for outsiders? Not run in this project (inputs unavailable) |
| 27. Gather | all per-(chr, row) JSONs | Sum privacy/utility, average rates | aggregated tables | Privacy is budgeted per chromosome but claimed genome-wide |

---

### In plain words

PanMixer takes a pangenome graph containing 44 real people, picks one of them, and rewrites that person's data in the released file so they can no longer be matched to an outside database — while leaving the file still useful to everyone else. It gets there in two halves.

The slow half, run once and taking (per the repository's README, line 78) "about 1–2 days to run end-to-end", converts the problem into arithmetic. It downloads the graph and two external reference datasets; strips chrX and the CHM13 sample — a human cell line derived from a complete hydatidiform mole, used as a reference assembly rather than contributed by a cohort donor — down to 44 individuals and 88 haplotypes; uses plink on a large external panel to find LD blocks, which are runs of variants inherited as a unit; turns every VCF into an integer array where −1 means "missing", *not* "reference"; builds allele frequencies with a three-way fallback that deliberately scores unseen variants as very rare; and then, for every block, computes three things — how surprising the target's real block is under the population model (the privacy value ε_j), one sampled plausible replacement block, and how much allele-frequency distortion installing that replacement would cause (the utility cost η_j).

The fast half runs many times. An LP-relaxed knapsack picks the best set of replacements within a per-chromosome utility budget; a stacker installs exactly those and copies everything else through untouched; and one column of the cohort VCF is rewritten. A suite of downstream evaluations then measures what was bought and what was lost: a rarity-weighted linkage attack whose failure point defines ε_private = 0.002 (published), a membership-inference check, allele-frequency and LD distortion (published: Wasserstein 0.006 for all alleles, 0.006 for SNPs only, 0.005 for rare SNPs with MAF < 0.05, against a 0.027 removal baseline; LD 0.003 against 0.017), a Beagle imputation attack testing reversibility, and `vg giraffe` read mapping testing whether the graph still serves outsiders (published: 77.82% perfect / 95.61% gapless / 77.01% MAPQ 60, against a baseline of 77.83 / 95.62 / 77.01).

On the block structure, the number to keep straight is this: on our rebuilt chr21 there are 100,757 block-dictionary entries, of which 89,087 — **88.4% of entries** — hold a single variant, and those account for **26.1% of variants**. The 10.4% of entries with two or more anchors hold 71.4% of variants. Never say "88.4% of variants".

Reading the released code against the published methods turns up five discrepancies, all recorded with evidence in this project's `docs/bugs.md`. The LD blocks and variant joins are computed against a GRCh37 panel while the graph is GRCh38, giving 1,025 chr21 matches instead of 258,610 — a 252× gap. The Li–Stephens sampler's distances are row-index differences rather than centiMorgans, with a reciprocal N_e, a stray factor of 4, and a non-anchor fill that compares base pairs against row indices; the switch probability collapses to about 10⁻¹² per step, and on 300 of 300 real chr21 blocks for HG00438 each "mosaic" was one donor haplotype copied whole. The privacy scorer includes the target in its own reference panel, so the model can always explain the block as a copy of the target, bounding the score by about log(88) = 4.4773 nats plus a small per-anchor term from the 10⁻⁴ mismatch probability — which is why a measured maximum of 4.508 sits just above 4.4773 without contradicting the argument — and on 150 real blocks the leave-one-out maximum rises from 4.508 to 78.1, truncating the rarest and most identifying blocks by up to 17×. The HMM's anchor set is bound to the PanGenie callset rather than the 1000 Genomes set the paper defines, though since PanGenie is GRCh38 this one is *not* compounded by the build mismatch. And the gap-score attack takes its maximum over the whole database without the self-exclusion the published equation (7) requires.

These are defects in the released code relative to the published methods, established by reading and running it. We have not reproduced the paper's figures and make no claim that the published results are wrong. Each is pedagogically useful in its own right, and together they are a compact tour of the three things that most often go wrong in this field: **coordinate systems, genetic distance, and who is and is not in the comparison set.**

What this section sets up for Section 4 is a single structural fact, not a list of bugs. PanMixer's guarantee is an *empirical threshold against one measured attacker*, and the "leave unselected blocks verbatim" design makes a worst-case distributional guarantee impossible in principle rather than merely absent in practice. Closing that gap — moving from a measured threshold to a bound that holds against every attacker at once — is what our own mechanism is for. What carries over from PanMixer is its **data**, its **block structure**, the algebra of its **forward recursion**, and its **evaluation suite**; what does not carry over is its **sampler**, its **privacy score** and its **release rule**. A rough sense of proportion, offered as an **estimate with no line-by-line counting method** rather than a measurement: most of the preprocessing and evaluation code is reusable, and almost none of the mechanism is.