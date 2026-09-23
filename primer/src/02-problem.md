## Section 2 — The Problem Itself

### 2.0 What this section is for, and what it deliberately does not do

Section 1 gave you the lay of the land. This section does something narrower and harder: it states **exactly what the problem is**, in enough mechanical detail that you could judge whether *any* proposed solution actually solves it.

There is a deliberate discipline here. **No solution appears in this section.** Not PanMixer (the published obfuscation tool dissected in Section 3), not our own proposed mechanism (Section 4), not differential privacy as a design, not even a hint of "and the fix is…". The reason is that in privacy research the single most common failure mode is *evaluating a method against a problem statement that was written after the method*. If you read the solution first, you will unconsciously accept its framing of the problem, and you will not notice the attacks it happens not to defend against. So: problem first, in full, with the attacker given every advantage. Solutions in Section 3 and Section 4.

By the end of this section you should be able to answer, unprompted: *what is released, who sees it, what do they already know, what can they compute, what does success look like for them, and what does the data owner lose if we simply refuse to release anything at all?*

#### 2.0.1 How numbers are sourced in this section

This section mixes three kinds of fact, and they are marked differently every time:

| Kind of fact | How it is marked | Where it lives |
|---|---|---|
| Published by the PanMixer authors | "the published paper reports…", usually with a direct quotation | `paper/s41467-026-77591-0_reference.pdf` (Nature Communications, DOI 10.1038/s41467-026-77591-0). The superseded preprint is `archive_docs/Blindenbach2026_PanMixer.pdf`; some numbers changed between the two, so the published version is always the one quoted. |
| Measured by us | "measured in this project" | The dated ledger `docs/memory.md`, the artifact registry `docs/data.md`, the defect list `docs/bugs.md`, and the raw run logs under `logs/`. |
| Read out of the source code | "read from the code" | The pinned PanMixer checkout at `external/PanMixer`, commit `c182c38d5bc8bb6f00f4b0b101207c4a009ca045` — MIT-licensed, 8,661 lines of Python. |

Our own mechanism proposal is `paper/Private_Genome_Path_Release.pdf`. Nothing in this primer should ever cite a temporary scratch path; if you find one, it is a bug in the document.

Two facts about arithmetic and background: general results from the literature that we did not re-derive (for example the Homer et al. 2008 result in §2.4.2) are labelled *standard field background* and given an author–year attribution, so you know to go and read the original rather than trusting this text. Anything we could not verify is phrased as an open question rather than asserted.

---

### 2.1 Vocabulary used in this section

Every term below is used later in this section. Read this table once; each term is also re-explained at the point of use, so you do not need to memorise it.

#### 2.1.1 Genomes, variants, files

| Term | Meaning in one or two sentences |
|---|---|
| **Base pair (bp), kilobase (kb), megabase (Mb)** | A base pair is one rung of the DNA ladder — one letter of sequence. 1 kb = 1,000 bp; 1 Mb = 1,000,000 bp. |
| **Chromosome** | One long continuous DNA molecule. Humans have **23 chromosome *pairs*** in a normal cell, and **24 *distinct* chromosome sequences** exist in the species: 22 numbered autosomes plus X plus Y. The two are different counts and are constantly confused. The often-quoted figure of about **3.1 billion bp is the *haploid* genome** — one copy of each chromosome. A normal diploid body cell contains roughly twice that much DNA, about 6.2 billion bp across two copies. |
| **Autosome** | Any chromosome other than X or Y. Chromosome 21, the chromosome all of this project's measurements are made on, is an autosome and one of the smallest. |
| **Locus** (plural *loci*) | A specific position, or small region, in the genome. "Position 5,031,442 on chromosome 21" is a locus. |
| **Allele** | One of the alternative DNA sequences that can occur at a locus. At a simple locus there are two: the **reference allele** (what the standard reference sequence says, written `REF`, numbered 0) and an **alternate allele** (`ALT`, numbered 1). |
| **Variant** | A locus at which people differ, together with the list of alleles seen there. |
| **Biallelic / multi-allelic** | A variant with exactly one ALT allele — two alleles in total, numbered 0 and 1 — is **biallelic**. One with two or more ALT alleles (0, 1, 2, …) is **multi-allelic**. This distinction is load-bearing: the attacker's standard scoring model in §2.4.1 is defined for biallelic sites only. |
| **SNP** (single-nucleotide polymorphism) | A variant where exactly one DNA letter differs, e.g. some people have A and others have G at that position. Pronounced "snip". |
| **Indel** | A variant that is a short insertion or deletion of DNA (conventionally under 50 bp). |
| **Structural variant (SV)** | A large variant — big insertions, deletions, inversions, duplications — conventionally 50 bp or more, sometimes megabases. |
| **Reference genome** | One agreed-upon consensus DNA sequence used as a coordinate system for everyone else. **GRCh37** and **GRCh38** are two successive versions ("builds") of the human reference. Coordinates are *not* interchangeable between builds: the same physical site has different numbers in each. |
| **Assembly** | The process of reconstructing long contiguous stretches of a genome directly from sequencing fragments, rather than only aligning fragments to a reference — and also the resulting reconstructed sequence. An assembly can contain sequence that is simply absent from the reference; an alignment-based call set mostly cannot. |
| **Reference bias** | The systematic error that arises because the reference is one sequence: DNA fragments that differ from the reference align worse, so variation far from the reference is under-detected. This is the main scientific motivation for pangenomes. |
| **Haplotype** | The set of alleles carried on a *single* physical copy of a chromosome, in order along that copy. Humans are diploid: two copies of each autosome, so two haplotypes per autosome, one inherited from each parent. |
| **Genotype** | What a person carries at a locus, counting both copies but not saying which copy is which: e.g. "one reference allele and one alternate allele". Written 0/0, 0/1, 1/1. |
| **Dosage** | The number of alternate-allele copies a person carries at a site: 0, 1 or 2. It is just the genotype written as a number. |
| **Phased / unphased** | Data is **phased** if we know which allele sits on which of the two haplotypes ("A is on the copy from mum, G on the copy from dad"). It is **unphased** if we only know the unordered genotype. Phased data is strictly more informative — and strictly more identifying. |
| **VCF** (Variant Call Format) | The standard text file format for variants: one row per variant, one column per sample, each cell a genotype. A VCF is essentially a big table of "who carries what". |
| **Cohort** | The set of individuals whose genomes went into building something. |
| **Pangenome** | A representation of many genomes at once, rather than one reference sequence. |
| **Graph / node / edge / path** | In a pangenome graph, **nodes** hold short DNA sequences, **edges** say which sequences can follow which, and a **path** is a walk through the graph spelling out one continuous sequence. Each haplotype in the construction cohort corresponds to a path. A *new* person's genome can be mapped onto the graph, producing *their* path. |

#### 2.1.2 Population structure

| Term | Meaning in one or two sentences |
|---|---|
| **Allele frequency (AF)** | The fraction of haplotypes in some population that carry a given allele. |
| **Minor allele frequency (MAF)** | The frequency of the *less* common allele at a variant; it ranges from 0 to 0.5. "Common" usually means MAF above about 5%; "rare" below about 1%. The published paper's "rare SNP" stratum is defined as population MAF < 0.05. |
| **Allele-frequency spectrum (AF spectrum)** | The histogram of how many variants fall in each frequency band — how many are very rare, how many common. It is the summary statistic population genetics cares about most, which is why distorting it is the headline utility cost in §2.6. |
| **Private variant** | A variant carried by exactly one individual in the data set at hand. Maximally identifying by construction. |
| **Recombination** | During the formation of eggs and sperm, the two parental copies of a chromosome swap segments. This shuffles alleles between haplotypes across generations. |
| **centiMorgan (cM)** | The unit of genetic distance. Two loci are 1 cM apart if there is a 1% chance per generation that a recombination event separates them. Genetic distance is *not* proportional to physical distance in base pairs — recombination rate varies enormously along a chromosome. |
| **Linkage disequilibrium (LD)** | The statistical correlation between alleles at nearby loci, created by shared ancestry and limited by recombination. Nearby variants are usually correlated; far-apart variants usually are not. LD is the single most important structural fact in this whole problem. Quantitatively, for two biallelic loci, D = p(AB) − p(A)p(B) measures the departure from independence, and r² = D² / [p(A)p(a)p(B)p(b)] rescales it to lie between 0 and 1: r² = 1 means one locus predicts the other exactly, r² = 0 means it tells you nothing. |
| **LD block / haplotype block** | A stretch of the chromosome within which variants are strongly correlated with each other and weakly correlated with variants outside it. Measured in this project on our rebuilt chromosome 21, the standard tool PLINK's `--blocks` caller finds **14,137** such blocks (`docs/memory.md`, entry dated 2026-09-18). |
| **Identity-by-descent (IBD)** | Two people share a segment *identity-by-descent* when both inherited it unbroken from the same recent common ancestor. Because recombination chops inherited segments down a little each generation, close relatives share long segments and distant ones share short segments — which is why relatedness is detectable from segment length alone, and why consumer genealogy works. |
| **Imputation** | Statistically filling in alleles you did not observe, using a reference panel plus LD. It is routine, accurate, and freely available — which makes it an attacker's tool as much as a scientist's. |
| **Association study** | A study that tests, variant by variant, whether carrying a particular allele correlates with a trait or disease across a large sample. It is the main scientific consumer of allele frequencies and LD, and therefore the reason AF and LD fidelity count as "utility". |

#### 2.1.3 The two 1000 Genomes releases (do not conflate them)

A **panel** (reference panel) is a large published collection of phased haplotypes from many people, used as a statistical reference — for imputation, for allele frequencies, or, adversarially, as an attack database. The **1000 Genomes Project** panel is the canonical example, and there are two different releases in play throughout this primer. Every number derived from a panel below says which one produced it.

| Release | Build | Samples | Role here |
|---|---|---|---|
| **1000 Genomes Phase 3** | GRCh37 | 2,504 | What the *shipped* PanMixer pipeline downloads (read from the code and from the repository's own data-availability pointers). It is on the older reference build. |
| **1000 Genomes 30x GRCh38 ("high-coverage")** | GRCh38 | 3,202 | What **this project substituted**, because the pangenome is on GRCh38. Measured in this project: **1,002,753 chr21 records**, of which 1,002,752 are unique (POS, REF, ALT) triples (`docs/data.md`). "30x" means each base was covered by about 30 sequencing reads on average — deeper, and therefore more accurate, than the roughly 4x of the older Phase 3 release. |

Why this matters immediately: coordinates are not interchangeable between GRCh37 and GRCh38, so joining a GRCh38 pangenome against a GRCh37 panel silently produces almost no overlap. Measured in this project on chr21: **1,025 exact (POS, REF, ALT) matches against the Phase 3 panel versus 258,610 against the 30x GRCh38 panel** — a 252-fold gap (`docs/memory.md`, `docs/data.md`). That measurement is Section 3's business; it appears here only so that when you see "the panel" you always ask *which one*.

#### 2.1.4 Statistics and privacy vocabulary

This primer promises you only basic probability, so here is the inferential-statistics vocabulary the rest of the section needs, in one place.

| Term | Meaning in one or two sentences |
|---|---|
| **Null distribution** | What a statistic would look like if the thing you are testing for were *not* true — here, what the attacker's scores would look like if the target were *not* in the data. |
| **z-score (sigma)** | (observed value − mean) ÷ standard deviation: how many standard deviations from typical an observation sits. "13 sigma" means thirteen standard deviations. |
| **False-positive rate** | The probability of raising an alarm when the hypothesis you are testing for is false. |
| **Power** | The probability of correctly raising an alarm when it *is* true. |
| **Simple hypothesis** | A hypothesis that fully specifies a probability distribution, as opposed to a family with unknown parameters. |
| **Likelihood-ratio test (LRT)** | The standard optimal way to decide between two simple hypotheses: compute how likely the observed data is under each, take the ratio, compare it to a threshold. The **Neyman–Pearson lemma** says no other test achieves higher power at the same false-positive rate. Attackers use LRTs; that is why they are the right thing to bound. |
| **Total variation (TV) distance** | A measure of how far apart two probability distributions are: TV(P,Q) is the largest possible gap, over all events A, between P(A) and Q(A). ("Largest over all events" is technically a *supremum*, i.e. a least upper bound; on a finite set of possible outputs, which is our case, it is just a maximum.) Operationally: if you must guess which of two distributions a sample came from, with equal prior odds, your best possible accuracy is (1 + TV)/2. TV = 0 means literally indistinguishable; TV = 1 means trivially distinguishable. |
| **Differential privacy (DP)** | A formal *definition*, not an algorithm: a randomised procedure M is epsilon-differentially private if, for any two datasets differing in one record and any set S of outputs, P(M(D) in S) is at most e^epsilon times P(M(D′) in S). Its point is that it bounds *every* attacker, not just the one you tested. |
| **Re-identification** | Attaching a name (or an external record) to a de-identified piece of data. |
| **Membership inference** | Determining whether a *specific known person* is present in a dataset, without necessarily learning anything else about them. |
| **Quasi-identifier** | A field that is not an identifier by itself but becomes one in combination with others (postcode + birth date + sex is the textbook example). |

#### 2.1.5 A symbol-collision warning you must carry through the whole document

This was flagged in Section 1's symbol table and is repeated here because Sections 2, 3 and 4 all use the same Greek letters for different things:

| Symbol | In the PanMixer paper | In our proposal | Elsewhere |
|---|---|---|---|
| **eta (η)** | *Utility loss* — the cost of an obfuscation move | **eta_tau**, a *tilt strength* parameter | — |
| **epsilon (ε)** | *Privacy risk* — the per-block score PanMixer minimises; `eps_private` is its empirical threshold | **epsilon_tau**, a DP-style *privacy budget* | epsilon is also the standard symbol in the DP definition above, and PanMixer's code separately uses 1e-4 as an HMM *mismatch probability* |
| **p, q** | Allele frequencies (reference and alternate) | Input and output paths | Also stay/switch probabilities in Section 4's HMM |

Whenever you meet one of these letters, re-read the sentence that introduced it in that subsection. The collision is not our choice; it is what happens when two literatures meet.

---

### 2.2 The precise setting: what object is released, by whom, to whom

Privacy questions are meaningless until you pin down three things: **what is handed over**, **who hands it over**, and **what the recipient already knows**. Two settings look superficially similar and are in fact fundamentally different. Confusing them is the most common conceptual error in this area, so we separate them carefully.

#### 2.2.1 Setting (i) — the contributor setting

A consortium collects genomes from a cohort of donors. It assembles them into a pangenome graph and publishes the graph (and usually a VCF derived from it) as a community resource. Everyone in the cohort **contributed their genome to the thing being published**.

Concretely, in the data this project works with: the Human Pangenome Reference Consortium (HPRC) v1.0 graph, built with the PGGB pipeline (a graph-construction toolchain) and converted to a VCF with `vg deconstruct -a` (a pangenome-graph toolkit command that reads variants back out of a graph relative to a chosen backbone sequence) against a GRCh38 backbone. Measured in this project and registered in `docs/data.md`: that VCF holds **45 samples**, and on chromosome 21 it has **340,824 unique (POS, REF, ALT) variant records**.

One of the 45 samples is **CHM13**. CHM13 is *human*: it is a human cell line derived from a complete hydatidiform mole — an abnormal conception in which the resulting tissue carries essentially one haplotype rather than two, which is exactly why it is prized as a reference assembly (it has no heterozygous sites to confuse an assembler). It is used here as a **reference assembly, not as a cohort donor**, so there is no individual whose privacy it represents. The shipped pipeline drops CHM13 and drops chromosome X, leaving **44 individuals = 88 haplotypes**. Every "2N = 88" in this document comes from that.

Here the person at risk *is inside the published object*. Their haplotypes are literally paths in the graph and columns in the VCF. Any sanitisation must therefore modify the resource itself, and every modification degrades the resource for everyone who uses it. The tension is immediate and structural: the individual's privacy and the resource's fidelity are made of the same bytes.

#### 2.2.2 Setting (ii) — the external-target setting

A pangenome graph has already been built and published from a cohort. It is *fixed*. Now a **new individual, who was never part of that cohort**, sequences their genome and maps it onto the published graph. The result is *their* path through *someone else's* graph. The question is whether that path can be released — to a collaborator, an archive, a clinical service, a method-development benchmark — without exposing them.

Here the person at risk is **outside** the published resource. Nothing is written back into the graph. The graph, the cohort, the cohort's paths, the mapping algorithm: all of it stays exactly as it was. The only new object in the world is the released path.

This is the setting our project targets (the mechanism itself is Section 4). It matters that it is genuinely different:

- In setting (i), the cohort's data *is* the utility. Degrading it harms third parties — every scientist using the reference.
- In setting (ii), the cohort's data is untouched and remains perfectly useful. The only thing being traded off is how faithfully the released path represents *this one person*.
- In setting (i), "just remove the person" is a coherent (if destructive) option. In setting (ii), "just remove the person" means *do not release anything*, which is not a privacy mechanism — it is a refusal.

#### 2.2.3 What is public and what is private

| | Setting (i): contributor | Setting (ii): external target |
|---|---|---|
| The pangenome graph | Public — it *is* the release | Public — pre-existing, untouched |
| The construction cohort's identities and haplotypes | Public (and are the thing at risk) | Public (and not at risk) |
| Reference genome (GRCh38), coordinates | Public | Public |
| External panels (Phase 3, 2,504 samples, GRCh37; 30x, 3,202 samples, GRCh38) | Public | Public |
| The sanitisation algorithm and its source code | Public (assume so) | Public (assume so) |
| All parameters and settings | Public (assume so) | Public (assume so) |
| The random seed actually used | **Secret** | **Secret** |
| The target's true genome / true mapped path | **Secret — this is the only secret** | **Secret — this is the only secret** |
| The released (sanitised) object | Public | Public |

Note the asymmetry in the last three rows. Exactly **one** substantive thing is private: the target's genome, equivalently the target's path. Everything else — including, importantly, the algorithm — is assumed known to the attacker.

**The random seed deserves its own line, and it is not a pedantic point.** A "seed" is the number that initialises a pseudo-random number generator; the same seed makes a randomised program replay exactly the same "random" choices. A randomised mechanism run with a known seed is therefore a *deterministic function* of its input. Read from the code: the standalone driver `tools/panmixer/obfuscate.py` in `external/PanMixer` takes an explicit `--seed` argument, documented in its own usage text as "optional random seed for HMM/allele-frequency sampling", and passes it straight to `np.random.seed`. Measured in this project (`docs/memory.md`, entry dated 2026-09-18; raw output in `logs/pm_obf.log`): running it twice on chr21 for subject HG00438 with the same seed produced byte-identical `new_haplotypes.npy` and `xsol.npy`, and two different seeds produced different outputs.

The consequence: the seed plays the role a cryptographic key plays. If the seed leaks, an attacker can run the public algorithm on candidate inputs until the output matches, and every empirical privacy claim evaporates. Seeds must be generated from a secure source, never logged, never published in a supplementary table.

Keep two related but distinct statements apart, because they are easy to read as a contradiction when you reach Section 4:

- *Operationally*, seed secrecy is necessary for any mechanism whose output is a deterministic function of (input, seed). That is the claim just made.
- *Formally*, whether a mechanism's guarantee is **derived from** its randomness is a separate question. A method can sample candidates and still rest its published claim entirely on measured attack outcomes rather than on the distribution it samples from. Section 3 shows that this is the case for PanMixer, and Section 4 contrasts it with a mechanism whose guarantee is a statement about the output *distribution*. Both statements can be true at once: "the guarantee is not proved from the randomness" does not mean "the seed may be published".

#### 2.2.4 Why "the algorithm is public" is the right assumption (Kerckhoffs-style reasoning)

**Kerckhoffs's principle**, from 19th-century cryptography, says a system must remain secure even if everything about it except the key is public knowledge. The modern slogan is "no security through obscurity". The reasoning transfers to genomic privacy almost word for word, and for four concrete reasons:

1. **Methods get published.** That is the point of doing science. The method is in a paper; the code is on GitHub under a permissive licence — the PanMixer implementation this project studies is MIT-licensed, 8,661 lines of Python, pinned here at commit `c182c38d5bc8bb6f00f4b0b101207c4a009ca045`. Assuming the attacker has not read the paper is assuming your reviewers cannot read either.
2. **Obscurity is not testable.** A claim like "the attacker probably does not know we used LD blocks" cannot be verified, cannot be falsified, and cannot be composed with anything else. A claim like "even knowing everything, the attacker's success is bounded by x" can be.
3. **Genomes are not revocable.** If a password leaks you change it. If a genome leaks, it leaks for life — and it partially leaks for your siblings, children, and parents, who never consented to anything. A privacy failure here has no remediation path, so the assumptions must be conservative.
4. **The realistic adversary is a legitimate recipient.** The person most likely to attack a released genomic dataset is not an anonymous hacker; it is someone who downloaded it lawfully — a competing lab, a re-identification researcher, an insurer's analytics vendor, a forensic service. They have the paper, the code, the parameters, and institutional patience.

There is one more reason specific to this field, and it is the sharpest: **the attacker's auxiliary data keeps getting better, and you cannot recall your release.** Panels grow — the 1000 Genomes panel itself went from 2,504 samples at Phase 3 to 3,202 in the 30x release — imputation improves, and direct-to-consumer genealogy databases expand. A guarantee that depends on the attacker's current resources expires. A guarantee that holds against an attacker who knows everything does not.

Note, finally, that the published PanMixer discussion concedes precisely this axis: "It is possible that, with additional genomic information or knowledge of the parameters and obfuscation process itself, an attacker could perform stronger attacks than those PanMixer protects against." That is an honest statement by the method's own authors, and it is the gap this whole primer is about.

---

### 2.3 Why a haplotype path leaks identity

It is worth being very concrete about *why* a string of alleles identifies a person, because the intuition "it's just a list of A/C/G/T, surely lots of people share it" is exactly backwards.

#### 2.3.1 Rare alleles

Some variants are carried by a handful of people on Earth. Carrying one is nearly a name badge. If a released path contains an allele that appears in, say, 1 in 50,000 haplotypes, and an attacker's database of 3,202 people contains exactly one person carrying it, that single fact roughly identifies the individual — subject, always, to the attacker's database containing them at all. The rarest variants, **private variants** carried by one individual, are maximally identifying by construction.

#### 2.3.2 Combinations of *common* alleles: do the arithmetic

The more counterintuitive and far more important fact is that **you do not need any rare variant at all**. Common variants, combined, are more than sufficient. Here is the arithmetic, explicitly.

**Version A — counting combinations.** Take n variants, each biallelic (two alleles, 0 and 1), each with two roughly equally common alleles, and suppose for the moment they are statistically independent of each other. A haplotype specifies one allele at each site, so there are 2^n possible haplotypes over those n sites. In plain words: each site is one yes/no question, and n independent yes/no questions carve the population into 2^n boxes.

- n = 10 → 1,024 combinations
- n = 20 → 1,048,576 (about a million)
- n = 30 → 1,073,741,824 (about a billion)
- **n = 33 → 8,589,934,592 (about 8.6 billion)**

The human population is roughly 8 billion. So **33 independent 50/50 sites already generate more distinct haplotypes than there are people alive.** Thirty-three. Not thirty-three thousand.

**Version B — expected number of coincidental matches, which is the quantity that actually matters.** Uniqueness is really about how many *other people* share your pattern, and "more boxes than people" does not by itself mean every box holds at most one person. So compute the expectation directly.

At a site where both alleles have frequency 0.5, two randomly chosen haplotypes carry the same allele with probability 0.5² + 0.5² = 0.5. Over n independent such sites, they agree at *every* site with probability 0.5^n. There are about 8 billion people, hence about 16 billion haplotypes (two per person). The expected number of haplotypes on Earth matching yours at all n sites is therefore

    E[matches] = 1.6 × 10^10 × 0.5^n

*In plain words: "how many haplotypes exist" times "the chance any one of them agrees with you everywhere".* This drops below 1 when n is about **34** (since log₂(1.6 × 10^10) = 33.9).

Repeat with *less* balanced, still-common variants — say MAF = 0.2, so allele frequencies 0.2 and 0.8. Two random haplotypes agree at such a site with probability 0.2² + 0.8² = 0.68. We need 1.6 × 10^10 × 0.68^n < 1, i.e.

    n > ln(1.6 × 10^10) / ln(1/0.68) = 23.50 / 0.3857 ≈ 61 sites.

So: **34 to 61 ordinary, common, unremarkable variants make a haplotype expected-unique on the planet.**

Now compare that to the scale of real data. Measured in this project on our rebuilt chromosome 21: the pangenome VCF has **340,824 variant records**, of which **274,458 (80.5%) are biallelic SNPs** — that count is the number of `True` entries in the `biallelic_snp_mask.npy` artifact that the shipped preprocessing builds, and it is the site set the attacker's scoring model of §2.4.1 is actually restricted to in the code. Chromosome 21 is one of the smallest human chromosomes, and there are 24 distinct chromosome sequences in the species. The 1000 Genomes 30x GRCh38 panel has 1,002,753 records on that same chromosome. You have between four and five *orders of magnitude* more identifying material than you need.

**The honest caveat**, which you should always attach to this calculation: real variants are *not* independent — that is exactly what LD means — so the effective number of independent sites is far smaller than the raw count. But "far smaller" here means thousands, not dozens. The margin is so enormous that no realistic correction rescues anonymity. This is why the field treats a haplotype as an identifier rather than as an attribute.

**The consequence to internalise:** a haplotype is a **quasi-identifier** — a field that is not an identifier by itself but becomes one in combination. Unlike a postcode or a birth date, a haplotype is a quasi-identifier with hundreds of thousands of components, every one of which is also a sensitive attribute in its own right. There is no clean separation between "the identifying columns" and "the sensitive columns". Every column is both. The published paper reaches the same conclusion in its own words, after testing a redaction baseline: "the unique combination of common alleles within a haplotype remains a quasi-identifier".

#### 2.3.3 Relatives

You share about half your genome with a parent, sibling, or child; about an eighth with a first cousin. Those shared stretches are shared **identity-by-descent**: both people inherited the same segment unbroken from a recent common ancestor, so the segments are long, and long shared segments are easy to detect. This is the entire basis of genetic genealogy. Two consequences follow:

- **A released path leaks about people who never consented.** Your sibling's risk goes up when you publish.
- **Removing you from a dataset does not remove your information** if a relative remains. Privacy definitions that assume records are independent are therefore mis-specified for genomes, and this is not a technicality — it is why off-the-shelf privacy definitions imported from tabular data do not transfer cleanly (§2.7, point 3).

#### 2.3.4 The attacker holds a database

Throughout, assume the attacker has an external genotype database for many people: a public panel, a biobank they have legitimate access to, a consumer-genealogy dump, or a forensic collection. The published PanMixer evaluation models exactly this: the adversary "has access to the released pangenome graph and possesses auxiliary information in the form of an external genotype database (e.g., the 1000 Genomes dataset)", and their objective is "to link a specific, obfuscated haplotype path within the graph to a known individual in their external database, thereby exposing their identity."

This is the crucial structural point: **the attacker does not need to know what the released path means. They only need to compare it to things they already have.** Privacy here is not about hiding meaning; it is about preventing a *match*.

One concrete observation about how thin the line between "cohort" and "attack database" really is. Measured in this project by intersecting sample identifiers (the 44 pangenome subject names in `pangenome_subjects.npy` against the sample columns of the 30x GRCh38 panel VCF): **39 of the 44 pangenome individuals appear by sample ID in the 3,202-sample panel** — including HG00438, the subject all of this project's per-subject measurements use. This particular intersection was computed while revising this section and belongs in the dated ledger `docs/memory.md` alongside the rest of the chr21 rebuild; treat it as reproducible in one line rather than as an established ledger entry until it is written there. Its significance for evaluation is taken up in §2.4.1 and in Section 3: when the "attacker's external database" and the "cohort being protected" share most of their members, what a gap score means needs care.

---

### 2.4 The attacks, mechanically

Three families, with worked micro-examples. Note first the standard taxonomy of what an attacker might achieve:

- **Identity disclosure (re-identification):** "this anonymous path belongs to Erin."
- **Membership disclosure (membership inference):** "Erin is in this dataset" — harmful whenever dataset membership is itself sensitive (a cancer cohort, an HIV cohort, a psychiatric-genetics study).
- **Attribute disclosure (reconstruction):** "whatever you hid, Erin's genotype at locus j is 1/1."

A defence can stop one and leave the others wide open. Evaluate them separately, always.

#### 2.4.1 Linkage / re-identification, and the "gap score"

**Mechanism.** The attacker takes the released object and scores it against every individual in their external database, using a score that rewards agreement and rewards it *more* when the agreeing genotype is rare. Then they look at whether the true individual comes first, and by how much.

**The score, in symbols and in words.** The published paper defines a genotype linkage score L between database individual g_i and the released target g*_target as

    L(g_i, g*_target) = − Σ_{j ∈ S(g_i, g*_target)} log f_j            (paper, Eq. 5)

where S(g_i, g*_target) is the set of loci at which the two genotypes are *identical*, and f_j is the population frequency of that genotype at locus j.

*In plain words:* walk along every locus; at each locus where the candidate's genotype exactly equals the released genotype, add −log(how common that genotype is). A common shared genotype (f_j near 1) adds almost nothing; a rare shared genotype (f_j tiny) adds a lot. Mismatches add nothing at all — this is a *reward-only* score, not a distance. The quantity −log f is called the **surprisal** of an event: how many nats of "that's unexpected" the observation carries. Higher total L means a better, rarer, more surprising match.

A trap for re-implementers, and a correction to an earlier draft of this primer: some plain-text extractions of the published PDF drop the leading minus glyph from Equation (5), leaving what looks like "L = Σ log f_j". Read from the typeset PDF (`paper/s41467-026-77591-0_reference.pdf`, page 9), **the minus sign is present**, and the equation therefore agrees with the paper's prose description of "the negative sum of the log-frequencies of the shared genotypes, ensuring that shared rare variants contribute more to the linkage scores". There is no sign discrepancy in the paper. There is a discrepancy between the paper and any extraction pipeline that eats the glyph — which matters, because implementing the sign backwards ranks every candidate in reverse and silently produces a "privacy-preserving" result.

Two derived quantities:

- **Success criterion (paper, Eq. 6):** the target is "linked" if L(true individual, release) is greater than or equal to L(any other database individual, release).
- **Gap score (paper, Eq. 7):** GapScore = L(true individual, release) − max over all *other* individuals i ≠ target of L(i, release). Positive gap = successful re-identification; the magnitude measures the attacker's confidence margin.

*(An implementation detail read from the code, recorded here only so you recognise the definition when Section 3 dissects it: in `tools/downstream/privacy/gap_score.py` the routine computes `highest_score_genotypes = np.max(scores_genotypes)` — an unrestricted maximum over every entry of the attack database, with no i ≠ target exclusion. Whether that departs from Eq. 7 in practice depends on whether the target is in the attack database at all; per §2.3.4, 39 of the 44 pangenome individuals do appear in the substituted 30x panel by sample ID. Section 3 works out what follows.)*

**Worked micro-example.** Six loci. For each, f_j is the population frequency of whatever genotype the release carries there, and the per-locus score contribution is −ln(f_j), measured in nats:

| Locus | f_j | −ln f_j |
|---|---|---|
| j1 | 0.50 | 0.6931 |
| j2 | 0.40 | 0.9163 |
| j3 | 0.25 | 1.3863 |
| j4 | 0.10 | 2.3026 |
| j5 | 0.01 | 4.6052 |
| j6 | 0.002 | 6.2146 |

The release is (initially) Erin's true genome, unaltered. The attacker's database holds Alice, Bob, Carla, Dan — and Erin.

| Candidate | Loci where they match the release | Score L |
|---|---|---|
| **Erin (true)** | j1, j2, j3, j4, j5, j6 | **16.1181** |
| Alice | j1, j2, j3 | 2.9957 |
| Bob | j1, j3, j4 | 4.3820 |
| Carla | j1, j2, j4, j5 | 8.5172 |
| Dan | j2, j3 | 2.3026 |

GapScore = 16.1181 − 8.5172 = **+7.6009**. Re-identification succeeds, decisively. Notice *where* the margin comes from: the two rarest loci, j5 and j6, contribute 4.6052 + 6.2146 = 10.82 of Erin's 16.12. Rarity dominates — which is precisely why the obvious first fix is "delete the rare stuff", and precisely why that fix fails (§2.5(b)).

#### 2.4.2 Membership inference

**Mechanism.** Here the attacker *already has* the target's full genome. They want to know whether that person is represented in a released dataset. The published evaluation models this by scoring the target's known genome against **all** paths in the released graph with the same rarity-weighted metric, and asking whether one candidate stands out: "A linkage attack is considered successful when the score distribution contains a clear outlier: one candidate genome scores significantly higher than all others." Conversely, the defence works when "the true obfuscated genome does not emerge as the top-scoring outlier."

**Worked micro-example.** The released graph contains 88 haplotype paths. The attacker scores Erin's known genome against each. Suppose the 88 scores have mean 5.0 and standard deviation 2.0.

- **Case A:** the highest score is 31.0. Its z-score is (31 − 5)/2 = **13 standard deviations** above the bulk. No plausible null distribution — that is, no plausible picture of what these scores would look like if Erin were *not* in the data — produces a 13-sigma maximum among 88 draws by chance. Conclusion: Erin is in. Membership disclosed.
- **Case B:** the highest score is 9.2, i.e. z = 2.1, and the second-highest is 9.0. Is a 2.1-sigma maximum surprising? No — and here is the reasoning, because "the maximum of 88 draws" is exactly the kind of thing intuition gets wrong. If you take 88 independent draws from a bell-shaped (standard normal) distribution, the largest of them sits, on average, near the 1 − 1/88 = 98.9th percentile, which for a normal distribution is about z = 2.28. So a maximum around 2.1 sigma is *below* what you would expect from pure chance. The attacker learns nothing.

The lesson: membership inference is an **outlier-detection** problem, not a matching problem. What protects against it is not being wrong, it is being *unremarkable*.

**A second, cheaper form of membership inference, from summary statistics alone.** Suppose you publish no genotypes at all — only per-variant allele frequencies of the cohort. Can that leak? Yes, and this is *standard field background*, not something measured here: the landmark result is **Homer et al. (2008)**, which showed that aggregate allele frequencies from a pool can reveal whether a known individual contributed to it. The result is a decade-plus old, is the reason NIH and the Wellcome Trust withdrew open access to GWAS summary statistics in 2008, and you should read the original rather than take this paragraph's word for it. What follows is an order-of-magnitude reconstruction of *why* it works, in our specific numbers.

With N = 44 individuals, a cohort has 2N = 88 haplotypes. At each variant, the published cohort AF is (number of alternate alleles) ÷ 88, and the target contributes 0, 1, or 2 of those 88 — that is the target's **dosage** at the site. So if the target is in the cohort, the cohort AF is *nudged toward the target's dosage*. The size of the nudge at one site is roughly

    (target dosage − 2 × background AF) / (2N)

*In plain words: the target's own two alleles are 2 out of 2N, so one person can move the published frequency by at most about 1/N — here at most 2/88 ≈ 0.023.* Individually this is invisible, buried in sampling noise.

But the attacker does not look at one site. They correlate the cohort's deviations-from-background against the target's deviations-from-background *across all sites at once*. Two facts about summing:

- The systematic part — the nudge, which always points the same way if the target is in — adds up **in proportion to M**, the number of sites.
- The noise — which points in random directions — adds up only **in proportion to sqrt(M)**, because independent random errors partially cancel.

So the signal-to-noise ratio grows like M / sqrt(M) = sqrt(M), scaled by the per-site signal strength of about 1/(2N):

    z ≈ sqrt(M) / (2N)

Plug in chromosome 21 alone, with M = 274,458 biallelic SNPs (measured in this project, §2.3.2): sqrt(M) = 523.9, and 2N = 88. That is a heuristic z-score of about **6** — comfortably detectable — from summary statistics only, on one small chromosome, for a 44-person cohort.

Treat the constant in that calculation as approximate rather than exact. The true power depends on the AF spectrum (sites where the minor allele is nearly absent contribute almost nothing) and is reduced by LD (correlated sites do not count as independent evidence, so the effective M is smaller than 274,458). **Establishing the exact power for *this* data is an open question we have not answered, and it is worth doing properly.** But the **scaling** is robust and is the point: *small cohort × many sites = leakage*, and pangenome cohorts are very small (44 individuals) while site counts are very large (hundreds of thousands per chromosome, and 24 distinct chromosomes in the species). "We only released frequencies" is not automatically safe; it is a claim that requires a calculation.

#### 2.4.3 Reconstruction / imputation attacks

**Mechanism.** The attacker's goal is to undo whatever was hidden, using LD plus a reference panel. This is not an exotic capability — it is the ordinary, well-engineered, publicly available machinery of statistical genetics, run adversarially.

It is worth seeing what imputation *is*, mechanically, because it unifies two halves of this primer. An imputer models an unobserved haplotype as a **mosaic of the panel's haplotypes**: it walks along the chromosome, at each step either continuing to copy from the same panel haplotype or switching to another, with switch probability set by the genetic distance between sites. That is the Li–Stephens hidden Markov model of Section 1 — the very same machinery the block-builder in Section 3 uses. Imputation is that HMM run to *fill in* what you did not see; the mechanisms in Sections 3 and 4 are the same HMM run to *generate* something to release. Same model, opposite direction. If you understand one you understand the other.

The canonical attack pipeline, as specified in the published Methods: take the released path, extract the SNPs that overlap the attacker's panel, run **Beagle** (a standard, freely available genotype-imputation and genotype-refinement program) with default parameters, using the 1000 Genomes data *with the target individuals excluded* as the reference panel, and compare the refined genotypes to the originals. The repository's own tool-fetch script (`external_tools/get_tools.sh`, read from the code) pins the **27Feb25 build** of Beagle, `beagle.27Feb25.75f.jar`; that is the artifact you would actually be running, and it is the one to cite rather than any marketing version number.

The critical control is also specified in the Methods: compare against "a baseline reconstruction obtained by selecting the most frequent alleles from the 1000 Genomes dataset". That baseline matters enormously for interpretation — an imputer that merely reproduces the major allele everywhere has learned nothing about the individual, so "imputation accuracy" is only meaningful *relative to the major-allele guess*. The paper's own claim is framed exactly this way: it reports that the attacker's reconstruction accuracy is "no better than simply guessing the major allele at each site".

Note finally that a successful reconstruction attack does not have to be perfect. It only has to be good enough to push the linkage score of §2.4.1 back above zero.

#### 2.4.4 Why LD makes naive redaction fail — worked example

This is the mechanical heart of the difficulty, so we do it with numbers.

**Setup.** Two nearby SNPs, A and B, each biallelic (alleles A/a and B/b). Among 100 haplotypes in the population, the four possible two-locus haplotypes occur with these counts:

| Haplotype | Count |
|---|---|
| A B | 45 |
| A b | 5 |
| a B | 5 |
| a b | 45 |

**Frequencies.** p(A) = 0.5, p(a) = 0.5, p(B) = 0.5, p(b) = 0.5, p(AB) = 0.45.

**LD coefficient.** D = p(AB) − p(A)p(B) = 0.45 − 0.25 = **0.20**. *In plain words: how much more often A and B travel together than they would if they were independent.* If they were independent you would expect 25 AB haplotypes; you see 45.

**Correlation.** r² = D² ÷ [p(A)p(a)p(B)p(b)] = 0.04 ÷ (0.5 × 0.5 × 0.5 × 0.5) = 0.04 ÷ 0.0625 = **0.64**. *In plain words: D rescaled to sit between 0 and 1, so it can be compared across sites with different allele frequencies.*

**Now hide B and keep A.** The attacker computes the conditional probability

    P(B | A) = p(AB) / p(A) = 0.45 / 0.50 = **0.90**.

So redacting locus B leaves it recoverable with 90% accuracy, against a 50% baseline from guessing the major allele. You deleted a column and hid essentially nothing.

**Generalise, with our measured block structure.** This argument holds at every scale: if a variant sits inside a group of correlated variants, hiding it while releasing its neighbours is close to a no-op. So the operative question is how big those correlated groups are. Measured in this project on our rebuilt chromosome 21 (`docs/memory.md`, entry dated 2026-09-18; run log `logs/pm_chr21.log`):

| Quantity | Value | What it counts |
|---|---|---|
| LD blocks called by PLINK's `--blocks` | **14,137** | Genuine multi-variant LD blocks, called on the panel |
| Entries in the working block dictionary (`blocks_dict.json`) | **100,757** | The units the mechanism actually operates on — the blocks *plus* one entry for every variant that fell in no called block |
| Singleton entries (one variant each) | **89,087 = 88.4% of *entries*** | Variants belonging to no called LD block. These same 89,087 variants are **26.1% of the 340,824 variants** |
| Entries with ≥ 2 "anchor" variants (the HMM path) | **10.4% of *entries*** | Blocks big enough that the mosaic model of §2.4.3 is used. They hold **71.4% of all variants** |

**Say which denominator, every time.** "88.4%" is a percentage *of dictionary entries*, not of variants. Writing "88.4% of variants are singletons" is wrong by a factor of more than three and contradicts the 71.4% figure in the very next row. The correct pair of sentences is: *88.4% of entries are singletons, and those singletons hold 26.1% of variants; 10.4% of entries take the HMM path, and those hold 71.4% of variants.* (An "**anchor**" here is a variant that also appears in the external call set the block-builder consults — which call set, exactly, is a subtle point Section 3 takes up, because the code and the paper do not name the same one. The remaining roughly 1.2% of entries are multi-variant blocks with fewer than two anchors, which are routed down a simpler allele-frequency path; 26.1% + 71.4% = 97.5% of variants, leaving about 2.5% — roughly 8,400 variants — in that residual category. The arithmetic closes.)

In other words: **most *blocks* are trivial, but most *variants* live in substantial multi-variant blocks** where exactly the recovery argument above applies.

**Why this specific fact drives the whole field.** Because correlation defeats per-variant reasoning, any defence must reason about *correlated groups* — blocks, or whole paths — rather than about individual variants. That is a constraint on solutions, not a solution; Section 3 takes it up.

#### Step table — the attacker's pipeline

| Step | Input | What happens | Output | Why it matters |
|---|---|---|---|---|
| 1. Obtain release | Public archive / collaborator | Attacker downloads the released path, graph, or VCF | A set of alleles attributable to one anonymous individual | Everything downstream is free once this exists; releases cannot be recalled |
| 2. Obtain auxiliary data | Public panels, biobanks, genealogy databases | Attacker assembles genotypes for many named people (e.g. the 3,202 samples of the 30x GRCh38 panel, or the 2,504 of Phase 3) | An external database indexed by person | The attack is comparison-based; without a database there is no one to link *to* |
| 3. Harmonise coordinates | Release + database | Match variants by (position, REF, ALT); requires both sides on the same reference build | Set of loci comparable across both | A build mismatch silently destroys overlap — measured in this project, 1,025 exact chr21 matches against the GRCh37 Phase 3 panel versus 258,610 against the GRCh38 30x panel, a 252x gap (`docs/memory.md`, `docs/data.md`) |
| 4. Reconstruct what was hidden | Harmonised release + panel | Run imputation (Beagle, 27Feb25 build, with the targets excluded from the panel), exploiting LD; score it against a major-allele baseline | A partially restored version of the release | Redaction and masking are undone here; this is what makes naive fixes fail |
| 5. Score every candidate | Restored release + database | Compute the rarity-weighted linkage score L = −Σ log f_j over loci where genotypes match | One score per database individual | Rare matches dominate; unique combinations of common alleles also accumulate |
| 6. Rank and measure margin | Scores | Take the top-scoring candidate; compute GapScore = L(true) − max over others | Positive or negative gap | Positive gap = re-identified. The *margin* measures confidence, not just success |
| 7. Membership test | Target's known genome + all released paths | Score the known genome against every released path; test whether one is a clear outlier against the null | "In" / "not in" verdict | Membership alone can be the harm, entirely independent of re-identification |

---

### 2.5 Why the obvious fixes fail

Every one of these is the first idea a competent person has. Each fails for a specific, diagnosable reason. Knowing *which* reason is what lets you evaluate real proposals.

**(a) Remove the individual entirely.**
In setting (i) this is the honest upper bound on damage. It is used in the published evaluation precisely as the "Removed" baseline, where "target individuals were entirely excluded from the graph", described as representing "the maximum possible utility loss for a given individual". The published numbers make the comparison concrete: at the paper's privacy threshold eps_private = 0.002, the Wasserstein AF divergence (defined in §2.6) is **0.006 for all alleles, 0.006 for SNPs only, and 0.005 for rare SNPs (population MAF < 0.05)**, against **0.027 for the Removed baseline** — roughly five times larger. For pairwise LD, the L1 difference is 0.003 at eps_private versus 0.017 for removal.

In setting (ii) removal is not even a mechanism: removing the person means "publish nothing", which delivers zero utility by construction.

Beyond utility, removal also fails on privacy grounds: (1) **relatives remain**, and they carry much of the same information (§2.3.3); (2) **differencing attacks** — if version 1 of a resource contains the person and version 2 does not, subtracting the two reveals them exactly; (3) if cohort membership is published separately, as it typically is, removal does not conceal membership either. Removal is a baseline, not a defence.

**(b) Remove only rare or private variants.**
This is the intuitively appealing fix, since §2.4.1 showed that rare loci dominate the linkage score. It does not work, for three independent reasons.

*Reason 1 — the arithmetic of §2.3.2.* Combinations of common alleles are already sufficient to be unique on Earth (34 to 61 sites). Deleting rare variants leaves that intact.

*Reason 2 — measured by the paper's authors.* The published paper tested exactly this baseline: "a baseline strategy that only removes unique or private variants while leaving common variation intact". It reports that the strategy "failed to prevent re-identification, maintaining a positive Gap Score comparable to the original graph", concluding that "the unique combination of common alleles within a haplotype remains a quasi-identifier".

*Reason 3 — our own worked example, continued.* Return to the six loci of §2.4.1 and redact the two rarest, j5 and j6, replacing them with the common genotype at each. **State the assumption explicitly, because it changes the arithmetic:** we assume that the common genotype substituted at j5 and j6 is *not* the genotype Alice, Bob, Carla or Dan carries there, so no rival gains a new match — this is the assumption most favourable to the defender, and the point is that the defence fails even so.

Erin no longer matches at j5 or j6, so her score falls to 0.6931 + 0.9163 + 1.3863 + 2.3026 = **5.2983**. But Carla also matched at j5, by carrying the same rare genotype, so *her* score falls too, to 0.6931 + 0.9163 + 2.3026 = 3.9120. The new best rival is therefore Bob, unchanged at 4.3820. GapScore = 5.2983 − 4.3820 = **+0.9163 — still positive.** Re-identification still succeeds.

The general lesson is worth stating separately: redacting rare variants lowers *everyone's* score, including the impostors', because rare matches were the impostors' scores too. The *ranking* — which is all the attacker needs — is remarkably robust to it.

**(c) Add independent random noise, per variant.**
Flip each allele independently with probability q. This fails from both ends simultaneously.

*Privacy side.* Independent noise on correlated signal is the easiest thing in the world to remove — that is what a denoiser is. Recall §2.4.4: if B is recoverable from A at 90% accuracy, noising B independently just gets undone by the same conditional-probability calculation. Worse, the linkage score aggregates over hundreds of thousands of sites, and aggregation crushes independent noise.

Here is that second point as arithmetic. *The agreement rates used here are chosen to make the arithmetic legible; they are illustrative, not measured.* Over 10,000 sites, suppose the noised release agrees with the true individual at 90% of sites and with an unrelated individual at 60%. The expected counts are 9,000 and 6,000, a difference of 3,000. How big is the random wobble in such a count? For a count of n independent yes/no events each with probability p, the variance is n·p·(1−p); taking p = 0.6 gives 10,000 × 0.6 × 0.4 = 2,400, so the standard deviation is sqrt(2,400) ≈ 49. A gap of 3,000 against a wobble of 49 is roughly **60 standard deviations** of separation. Independent per-site noise at any survivable rate does not move the ranking.

*Utility side.* Independent noise destroys precisely the structure that makes the data valuable. LD *is* the correlation between sites; independent flips *reduce* every pairwise correlation toward zero. AF spectra are wrecked, because flipping common alleles at rate q manufactures spurious rare variants everywhere. And in a *graph*, per-variant noise can produce paths that are not valid traversals at all — sequences that no chromosome could carry. The published discussion makes exactly this point about naive DP-style approaches: "Standard DP mechanisms would require perturbing edges or introducing dummy paths, which breaks the long-range LD structure essential for downstream tasks such as variant calling and read mapping", and notes that the noise required to satisfy privacy bounds "destroys the structural integrity of the genome assemblies".

**(d) k-anonymity-style grouping.**
**k-anonymity** requires that every record be indistinguishable from at least k−1 others with respect to the quasi-identifiers — typically achieved by generalising them (postcode → region, age → decade). It is the workhorse of tabular de-identification and it fails structurally here:

1. **There is no quasi-identifier/sensitive-attribute split.** k-anonymity assumes you can name the identifying columns, generalise them, and leave the sensitive column alone. In genomics *every* variant is simultaneously identifying and sensitive. There is nothing left to protect once you have generalised the identifiers, because the identifiers *are* the payload.
2. **Generalisation at this dimensionality destroys everything.** To make a 274,458-site haplotype match k−1 others you must coarsen until essentially all individual-level signal is gone — which returns you to (a).
3. **k-anonymity does not bound inference.** Even if k people share a group, if they all carry the same allele at a locus, the attacker learns it for certain. (This is the classic homogeneity attack, and it is endemic in genomics because LD *makes* groups homogeneous.)
4. **It assumes bounded background knowledge**, which contradicts §2.2.4 — and background knowledge only grows.

**(e) Release only summary statistics.**
Two objections. *Privacy:* §2.4.2's scaling argument — a heuristic z-score around 6 on chromosome 21 alone for a 44-person cohort, in the same family as the standard Homer et al. (2008) result — says that per-site frequencies from small cohorts are not obviously safe against membership inference. The safe version of that sentence is that safety requires a power calculation rather than an assumption, and we have not done that calculation for this data. *Utility:* summary statistics are the wrong object entirely. You cannot align sequencing reads to a table of allele frequencies. You cannot genotype a new sample against it. You cannot recover haplotype structure, which is the thing a pangenome exists to provide. Releasing only summaries is not a privacy-preserving version of the release; it is a different, much less useful product.

**(f) A note on access control.** "Put it behind a data-access committee" is a real and widely used answer, and it is deliberately out of scope here — it changes *who* the recipient is, not what the released object reveals to that recipient. Every attack above is carried out by an authorised downloader. Access control and mechanism design are complements, not substitutes, and a primer about mechanisms should not pretend otherwise in either direction.

| Fix | Which attack still works | Which utility dies |
|---|---|---|
| Remove individual | Membership (via relatives, via differencing across versions, via separately published cohort lists) | Maximal: published Wasserstein AF divergence 0.027 for the "Removed" baseline vs 0.006 at eps_private; in setting (ii) there is no release at all |
| Remove rare variants only | Linkage — the paper measured a positive Gap Score "comparable to the original graph"; our worked example keeps gap = +0.9163 | Modest; but privacy gain is near zero, so the trade is bad at any price |
| Independent per-variant noise | Linkage (noise averages out over 10⁵ sites — about 60 sigma in the illustration above); reconstruction (LD denoises it) | LD structure, AF spectrum, graph/path validity |
| k-anonymity grouping | Attribute disclosure within homogeneous groups; any attacker with background knowledge | Everything, once generalisation is strong enough to reach k |
| Summary statistics only | Membership inference from AF shifts (scaling ≈ sqrt(M)/2N ≈ 6 on chr21) | Read mapping, genotyping, haplotype structure — the entire use case |

---

### 2.6 What "utility" even means here, and why it is contested

There is no single number called utility. There are at least four distinct notions in play, they are not equivalent, and they can be individually satisfied while the release is useless or unsafe.

**1. Allele-frequency (AF) fidelity.** Does the released resource preserve the distribution of allele frequencies — the AF spectrum? The measure used in the published work is the **Wasserstein divergence** between the AF distributions before and after. *In plain words:* the "earth-mover's distance" — the minimum amount of probability mass you would have to shift, and how far, to turn one distribution into the other. Small is good; the scale here runs 0 to 1.

One arithmetic caution that trips people up when they compare §2.6 with the paper's Methods: for alleles modelled as Bernoulli (0/1) random variables, the paper states that "the Wasserstein divergence is two times the mean absolute difference between AFs before and after obfuscation. The factor of two arises because each diploid individual carries two alleles per variant." So if you compute mean |p_j − q_j| yourself and compare it to a published divergence, you will be off by a factor of two and conclude one of the two is wrong. Neither is; the convention just carries the diploid factor.

Published calibration points (`paper/s41467-026-77591-0_reference.pdf`), all at the paper's privacy threshold eps_private = 0.002, with their labels attached because a number without its stratum is not a fact: **0.006 for all alleles, 0.006 for SNPs only, 0.005 for rare SNPs (population MAF < 0.05)**, against the **"Removed" baseline of 0.027**. Across the whole sweep of privacy levels the paper reports AF divergence staying below 0.008.

**2. LD fidelity.** Does the released resource preserve the correlation structure between nearby variants? Measured in the published work as the **L1 loss** — the sum of absolute differences — between LD matrices computed before and after, where LD is pairwise r² between SNPs within sliding 5 kb windows. Published values: 0.003 on average at eps_private, versus 0.017 for the complete-removal baseline. This matters because LD is what makes imputation, association testing, and haplotype inference work at all.

**3. Read-mapping performance.** Does the graph still do its job — align sequencing reads well, without reference bias? Measured as the percentage of reads that map **perfectly** (exactly, no mismatches), **gaplessly** (no insertions or deletions in the alignment), or at **MAPQ 60** (the aligner's top confidence score, meaning "essentially certain this is the right place"). Published values, again with their labels: **77.82% perfect, 95.61% gapless, 77.01% MAPQ 60**, against an unobfuscated baseline of **77.83 / 95.62 / 77.01**. This is the applied reason pangenomes exist, so it is arguably the most important metric — and it is the least mathematically tractable of the four.

**4. Fidelity to the target specifically.** Does the released object still tell you about *this person*? This is a completely different question from 1–3, and here is the observation that makes the whole evaluation problem subtle:

> **Cohort-level utility metrics can be satisfied essentially perfectly by an output that contains no information about the target whatsoever.**

Consider a "mechanism" that ignores the target entirely and outputs a randomly chosen haplotype from the public cohort. What happens to metrics 1, 2 and 3?

- Metric 1: the AF spectrum moves by at most one haplotype out of 88 — a tiny perturbation, of the same order as the published 0.006.
- Metric 2: LD is *perfectly* preserved, because the output is a real cohort haplotype and therefore carries exactly the real correlation structure.
- Metric 3: reads map beautifully, because the output is a genuine biological sequence that some human actually carries.

This fake mechanism scores well on every cohort-level metric, is perfectly private by construction (it never looked at the target), and is completely worthless.

Three consequences, and they should shape how you read every evaluation in this field:

- **Cohort-level metrics measure "did you damage the shared resource", not "did you keep the person".** They are necessary, never sufficient.
- **Any claim of a privacy–utility trade-off must state which utility.** A curve trading privacy against AF fidelity is measuring the wrong axis if the point of the release was target fidelity.
- **Target fidelity must be measured against a null.** "The release agrees with the target at 78% of sites" is meaningless until you know what a random cohort haplotype would score. If the null is 76%, you have released essentially nothing.

Utility is also *contested* in a second sense: different downstream users want different things. A population geneticist running association studies wants AF and LD fidelity. A read-mapping engineer wants sequence realism and graph topology. A clinician wants the target's specific alleles at specific pharmacogenomic loci to be exactly right, and would happily sacrifice AF fidelity to get it. There is no scalar that serves all three, so every method embeds a value judgement in its choice of metric. Read that choice as an argument, not a measurement.

---

### 2.7 The shape of the difficulty

Step back from the mechanics. Five structural facts explain essentially every difficulty above.

**1. Utility for the target and distinguishability of the target are the same quantity.**
If a released object helps anyone do *anything* about the target that they could not do by chance, then the distribution of that object depends on the target. If the distribution depends on the target, a likelihood-ratio test has power above chance to detect that dependence. There is no configuration in which the release is informative about the person and simultaneously statistically independent of them. This is not an engineering limitation to be designed around; it is a tautology. The only real question is *how much* dependence, measured how — and the only defensible answer is a quantitative one.

**2. A guarantee that holds against one tested attacker is weak evidence.**
An empirical evaluation says: *against this database, with this score function, at these loci, with this imputer, no attack succeeded.* It does not say *no attack succeeds.* Four specific gaps:

- The attacker's database changes. Panels grow — 2,504 samples at Phase 3, 3,202 in the 30x release — and consumer genealogy databases grow faster.
- The score function changes. The rarity-weighted score of §2.4.1 is one choice; a learned discriminator is another.
- The attacker may know the mechanism. The published work states this explicitly as a limitation: "It is possible that, with additional genomic information or knowledge of the parameters and obfuscation process itself, an attacker could perform stronger attacks than those PanMixer protects against."
- Empirical thresholds are not portable. The published work is equally explicit that eps_private "is not universal as it depends on the specific graph, external genotype database, and target individuals."

The published discussion names the resulting gap directly: "Providing formal, worst-case guarantees against membership inference remains a distinct and important direction for future work." That sentence is, in effect, the statement of the open problem this project exists to address — and it is a statement by the authors of the method, not a criticism from outside.

This is where **total variation distance** enters as a *standard of evidence* rather than as a technique. Saying "the output distributions induced by two different inputs are at most tau apart in TV" is a statement about *all* distinguishers at once, because TV is defined as the largest gap over *all* events (formally a supremum; on a finite output set, a maximum). It converts "we tried some attacks" into "no attack, present or future, can exceed accuracy (1 + tau)/2 at equal priors". Similarly, **differential privacy** is valuable not because of any particular noise-adding algorithm but because of the *shape* of its promise: worst-case over datasets, over outputs, and over adversaries, with composition rules telling you what happens when you release more than once.

**3. But importing worst-case definitions is itself hard, because genomic records are not independent.**
DP is defined over datasets that differ "in one record". For genomes that phrase is ambiguous in at least three ways:

- One person contributes two haplotypes and hundreds of thousands of correlated variant calls, not one record. What is "one record" — a variant? a haplotype? a person? a family?
- One person's data partially determines their relatives' data (§2.3.3), so removing a record does not remove the information.
- The utility the release must preserve *is* the correlation structure, which is precisely what generic noise destroys (§2.5(c)).

The published discussion makes the practical version of this point: the noise required to satisfy standard bounds "destroys the structural integrity of the genome assemblies". So the difficulty is not that formal guarantees are unavailable in principle — it is that the standard ones are mis-specified for correlated, high-dimensional, biologically constrained data, and adapting them is real work. That adaptation is Section 4's subject.

**4. High dimensionality plus correlation is the root cause of almost everything.**

- **High dimensionality** is why anonymity is impossible by default: 34 to 61 independent sites make you unique, and you have 274,458 biallelic SNPs on one small chromosome out of 24.
- **Correlation** is why redaction fails: neighbours predict what you hid (r² = 0.64 gave 90% recovery in §2.4.4).
- **Correlation** is why independent noise fails: the signal is correlated, the noise is not, and the difference is easy to exploit.
- **Correlation across people** — relatives — is why removal fails.
- **Correlation** is also why utility is fragile: the correlation structure *is* the utility, so any mechanism that degrades correlation to buy privacy degrades the product it is trying to protect.
- **Correlation** is finally why you cannot reason one variant at a time and must instead reason about blocks or whole paths — which is what makes the combinatorics hard and the algorithms non-obvious. And the block structure is not a convenient abstraction: measured in this project, 71.4% of chr21's variants sit in the 10.4% of dictionary entries that are genuine multi-anchor blocks.

**5. Two smaller but load-bearing points.**

- **Composition.** A bound for *one* release says nothing about ten. If the same person's path is released twice under independent randomness, an attacker can average the two releases, partially cancelling the randomness, and effective privacy degrades. Any guarantee must state its release budget, and a one-release-per-genome constraint is a real operational restriction, not a footnote.
- **The threshold is a policy question, not a mathematical one.** Suppose you can bound an attacker's accuracy at 55%. Is that acceptable? Nothing in the mathematics answers that. The mathematics produces a dial; where to set it is a decision about harm, consent, and context, and it belongs to people, not to the algorithm. Be suspicious of any paper — including ours — that presents a threshold as though it were derived.

---

### In plain words

We are trying to publish something derived from one person's genome without letting anyone work out whose genome it is, or whether that person is in a particular dataset, or what they were carrying at the parts we tried to hide. There are two different versions of the setting and they are easy to confuse: one where the person contributed to the published pangenome graph itself, and one where the graph is already public and the person is an outsider whose mapped path through it is the only private thing. In both, we assume the attacker has the graph, the cohort, the algorithm, the code and the parameters — everything except the target's genome and the random seed — because assuming otherwise is not testable, and because a leaked genome cannot be reissued the way a leaked password can. (The seed really is a key: we re-ran the shipped tool on chromosome 21 and the same seed gave byte-identical output while a different seed did not, which is logged in this project's dated ledger.)

A genome identifies you not through exotic rare mutations but through sheer combinatorics: about 34 balanced common variants, or about 61 slightly less balanced ones, already produce more combinations than there are people alive, and one small chromosome — chromosome 21, one of 24 distinct human chromosome sequences — gives you 274,458 biallelic SNPs among 340,824 variant records in this project's build. The attacks that exploit this are simple to run: score the release against every person in a public database with a rarity-weighted match score, which adds −log(how common that genotype is) at every locus where the candidate agrees, and see whether the true person wins and by how much — that margin is the gap score; or score a known genome against every released path and look for one score that sticks out from the other 87, which is membership inference and is an outlier-detection problem rather than a matching one; or run an ordinary imputation program, the same mosaic model used to build the mechanisms themselves, only run in reverse, to rebuild whatever was hidden from the correlated variants around it.

Every obvious defence fails for a specific, checkable reason. Deleting the person destroys their contribution — the published paper's own "Removed" baseline distorts allele frequencies about five times more than its method does, 0.027 against 0.006 — and still leaks through relatives and through differencing two versions of a resource. Deleting only the rare variants leaves the gap score positive, as both the published experiment and our own six-locus worked example show, because striking out rare matches lowers the impostors' scores as well as the true person's and leaves the ranking intact. Independent per-site noise is averaged away by a score that aggregates over hundreds of thousands of sites, while simultaneously wrecking exactly the correlation structure that made the data worth having and, in a graph, producing paths no chromosome could carry. k-anonymity has no quasi-identifier column to generalise, because in a genome every variant is both identifying and sensitive at once. And summary statistics alone both leak membership — the standard Homer et al. (2008) result, whose scaling argument gives a rough signal-to-noise ratio of about 6 for a 44-person cohort on chromosome 21 — and cannot be mapped to, genotyped against, or phased.

Utility is not one thing: allele-frequency fidelity, linkage-disequilibrium fidelity, read-mapping performance and fidelity to the target are four different goals, and the most important trap in the whole field is that the first three can be scored almost perfectly by an output that simply copies a random cohort haplotype and says nothing at all about the target. Cohort-level metrics tell you whether you damaged the shared resource; they do not tell you whether you kept the person.

The deep reason all of this is hard is that any output useful about a person is by definition statistically dependent on that person, so the question is never whether there is leakage but how much, measured against every possible attacker rather than the two or three we happened to try — which is why total variation distance and differential privacy matter here as *standards of evidence*, not as techniques. And the root of the practical difficulty is that genomic data is both very high-dimensional, which makes everyone unique, and highly correlated, which makes hiding any individual piece of it futile — so a defence has to act on whole correlated blocks or whole paths, not on variants one at a time. On our chromosome 21, 88.4% of the 100,757 block-dictionary entries are singletons holding just 26.1% of variants, while the 10.4% of entries that are real multi-anchor blocks hold 71.4% of them: most blocks are trivial, but most variants are not.