## Section 1 — The Landscape

This section builds the vocabulary and the mental models that everything later depends on. Read it slowly; nothing after it will make sense if a term here is still fuzzy.

**What this section assumes you already know.** That DNA is a long molecule spelled in four letters (A, C, G, T); that a gene is a stretch of DNA; that chromosomes are long DNA molecules; that a mutation is a change in DNA. On the computing side: arrays, hash tables, graphs, big-O notation, and probability at the level of conditional probability, Bayes' rule, and expectation.

**What this section does *not* assume, and therefore defines from scratch.** Anything about how DNA becomes computer-readable data (reads, alignment, assembly, coverage). Anything about population genetics (allele frequency, drift, selection, effective population size, linkage). Anything about the bioinformatics tool ecosystem. Anything about hidden Markov models. Anything about privacy, security, hypothesis testing, or information theory. If one of those words appears below without a definition attached, that is a bug in this document, not a gap in your background.

**The one-sentence version of the whole problem:** *a public genomic resource built out of real people's genomes — or used to describe a new person's genome — can leak enough information to identify those people, and we want to release something useful anyway.* Everything below is scaffolding for that sentence.

---

### 1.0 Notation, and a warning about symbol collisions

This document compares two pieces of work that were written independently and that **use the same Greek letters for different things**. If you do not know this in advance you will read Section 3 and Section 4 as contradicting each other when they are simply speaking different dialects. The collision is flagged here, at the front, and every later section re-states which dialect it is in.

The two bodies of work are:

- **PanMixer** — the published obfuscation tool this project dissects in Section 3. Published as an "Article in Press" in *Nature Communications*, DOI 10.1038/s41467-026-77591-0; the copy of record used here is `paper/s41467-026-77591-0_reference.pdf`. An earlier, superseded preprint is kept at `archive_docs/Blindenbach2026_PanMixer.pdf`; its numbers differ from the published ones and must not be quoted. Its source code is pinned in this repository as the checkout `external/PanMixer` at commit `c182c38d5bc8bb6f00f4b0b101207c4a009ca045`.
- **This project's own proposal** — *Private Genome Path Release against a Public Pangenome Graph*, at `paper/Private_Genome_Path_Release.pdf`, the mechanism taken apart in Section 4.

| Symbol | In PanMixer's paper / code | In this project's proposal | Elsewhere in this document |
|---|---|---|---|
| **epsilon (ε)** | ε_j = the **privacy risk** of block *j*, defined as −log p(h_bj); ε_private = 0.002 is the paper's chosen operating threshold | ε_τ = 2·arctanh(τ), a **differential-privacy-style budget** in log units | ε is also the standard symbol for the DP **privacy budget** in §1.3.6; and in Section 3's defect discussion ε names the HMM's **per-anchor mismatch probability** (1e-4), which is none of the above |
| **eta (η)** | η_j = the **utility loss** caused by obfuscation move *j* | η_τ = arctanh(τ), the **tilt strength** — how hard the mechanism is allowed to lean toward the target | — |
| **p, q** | allele frequencies at a biallelic site (p = reference-allele frequency, q = 1 − p) | **paths**: p is the private input path, q another admissible input path | in Section 4's transition model p and q are the **stay** and **switch** probabilities of the HMM |
| **r** | r = 1.26, a **recombination-rate constant** in the paper's Fig. 6 transition formula | — | r² (§1.1.4) is the **squared correlation** between two loci — a completely unrelated quantity that merely shares a letter |
| **N** | N = number of individuals in the cohort (44 here), so 2N = number of haplotypes (88) | same | — |
| **T, K** | — | T = number of positions along the chromosome (sites or blocks); K = number of hidden states (haplotypes) | same, throughout |

Two rules follow. First, **never carry a symbol across the boundary between the two works** without re-reading its definition. Second, whenever this document quotes a number, it also says which work and which quantity it belongs to.

**Where the numbers in this document come from.** Three different kinds of number appear, and they are always labelled:

1. **Published results** — quoted from `paper/s41467-026-77591-0_reference.pdf`, always with the label the paper attaches (which metric, which baseline, which threshold).
2. **Measured in this project** — produced by running code here, on data here. Every such number is recorded in the repository's permanent documentation: `docs/memory.md` is the dated ledger of what was run and what came out, `docs/data.md` is the artifact registry (what each file is, its size and checksum), `docs/bugs.md` records defects with their evidence, and the raw run output lives under `logs/`. This document cites those locations, never a temporary scratch path.
3. **Standard field background** — textbook or landmark results, given with an author–year attribution and explicitly marked as background rather than as something measured here.

---

### 1.1 The biology you need

#### 1.1.1 Genome, chromosome, locus, variant, allele

A **genome** is the complete DNA sequence of an organism. DNA is double-stranded, and the two strands pair letter-for-letter (A with T, C with G), so one rung of the ladder is called a **base pair**, abbreviated **bp**. Two size units recur constantly and are worth fixing now: a **kilobase (kb)** is 1,000 bp and a **megabase (Mb)** is 1,000,000 bp.

Human genome size has to be stated carefully because two different figures are both correct:

- The **haploid** genome — one copy of each chromosome — is about **3.1 billion bp (3.1 Gb)**. This is the number usually meant by "the size of the human genome," and it is the size of a reference sequence.
- A normal human cell is **diploid**: it carries *two* copies of each non-sex chromosome. Its total DNA content is therefore roughly **twice** that, about 6.2 billion bp.

The genome is not one molecule. It is packaged into **chromosomes**. The count also has two correct forms that are easy to confuse:

- A human karyotype contains **23 pairs** of chromosomes.
- Those pairs are drawn from **24 distinct chromosome sequences**: 22 numbered **autosomes** (chromosome 1 through chromosome 22, present as a pair in everyone) plus the two sex chromosomes **X** and **Y**. A reference assembly must therefore contain 24 distinct sequences, even though any one person carries at most 23 distinct ones.

Diploidy is the single most important structural fact for this project. When we say "a person's chr21," we always mean *two* physical chr21 molecules — one inherited from the mother, one from the father — and they are not identical to each other.

A **locus** (plural **loci**) is a location in the genome: a coordinate, or a short coordinate range. An **allele** is one of the alternative DNA sequences that can occur at a locus. A **variant** is a locus at which people are observed to differ, together with the list of alleles seen there.

Concretely: at some position on chr21, most people have a `C` and some people have a `T`. That position is a variant; `C` and `T` are its two alleles. By convention one allele is designated the **reference allele** (**REF**) and the others are **alternate alleles** (**ALT**), and they are numbered: allele `0` is REF, allele `1` is the first ALT, allele `2` the second, and so on. This numbering is the reason a cohort of genomes can eventually be stored as a matrix of small integers.

Two words for the shape of that allele list, both load-bearing later:

- A variant with exactly one ALT allele — two alleles total, `0` and `1` — is **biallelic**.
- A variant with two or more ALT alleles is **multi-allelic**.

The distinction matters because several standard pieces of machinery, including the genotype model used by the linkage attack in Section 3, are only defined for biallelic sites.

Variants also come in flavours distinguished by size and shape:

| Name | Abbreviation | What it is | Typical size |
|---|---|---|---|
| Single nucleotide polymorphism | **SNP** | One base substituted for another (C→T) | 1 bp |
| Insertion / deletion | **indel** | A few bases inserted into, or deleted from, the sequence | 1–49 bp (by convention) |
| Structural variant | **SV** | A large rearrangement: deletion, duplication, inversion, large insertion, translocation | ≥50 bp, up to megabases |

SNPs are by far the most numerous and the easiest to handle; a single individual differs from the reference at millions of SNP positions. SVs are far rarer per person but touch far more total base pairs, and historically they were the hardest to detect — which matters, because the pangenome (§1.1.6) exists in large part to represent them properly.

*Measured in this project:* on our rebuilt chr21 there are **340,824 variant records**. The chr21 artifact set is registered in `docs/data.md` and the preprocessing run is logged in `logs/pm_chr21.log`.

**274,458 of those records are biallelic SNPs — 80.5%.** "Biallelic SNP" means a single-base `REF` and exactly one single-base `ALT`, and this figure is the count of `True` entries in the biallelic-SNP mask the shipped preprocessing builds (`build_biallelic_snp_mask`): *measured in this project*. Place it among its neighbouring counts, because they are close together and easy to confuse. **65,176** records are indels or structural variants (measured in this project), so **275,648** records are single-base ones (derived: 340,824 − 65,176); **1,190** of those are *multi-allelic* SNP records — a single-base `REF` with two or more single-base `ALT`s, which the mask excludes — leaving **274,458** (derived). One number is a subset of the other; they do not have to be equal. Two further counts on this same chromosome are **different filters** and must never be quoted as the mask count: **288,020 records (84.5%)** is what `bcftools view -v snps` keeps, a looser filter that retains mixed multi-allelic records, and **238,052 records (69.8%)** is what the gap-score attack keeps after applying the mask *and* requiring a match in the 1000 Genomes panel (both measured in this project). Three filters, three counts; only **274,458** is the mask.

#### 1.1.2 Genotype versus haplotype, and what "phased" means

Because you carry two copies of each autosome, at any variant locus you carry *two* alleles.

- A **genotype** is the *unordered pair* of alleles you carry at a locus. Written `0/0` (two reference copies, "homozygous reference"), `0/1` ("heterozygous"), `1/1` ("homozygous alternate"). The slash `/` means *we do not know which physical chromosome copy each allele sits on*.
- A **haplotype** is the *ordered sequence of alleles along one single physical chromosome copy*. If we know that the `1` at locus A and the `1` at locus B sit on the *same* copy, we have haplotype information. Written with a pipe: `0|1` means "reference allele on copy 1, alternate allele on copy 2," and the copy numbering is consistent across every locus on that chromosome.

**Phasing** is the process of determining which allele goes on which copy — converting genotypes into haplotypes. The name comes from saying that the two chromosome copies are "in phase."

*Why phasing is hard.* If a person is `0/1` at locus A and `0/1` at locus B, there are two possibilities: the two `1`s are together on one copy (haplotypes `1 1` and `0 0`), or they are on opposite copies (haplotypes `1 0` and `0 1`). Standard short-read sequencing (§1.1.5) reads only a few hundred bases at a time, so if A and B are 10,000 bases apart, no single read spans both and the raw data simply does not say which arrangement is true. Phase is recovered by one of three routes:

1. **Long-read or assembly-based sequencing**, which spans many variants in one contiguous piece of evidence. **Assembly** is the process of reconstructing long contiguous stretches of a genome directly from sequencing fragments, rather than only aligning fragments onto a pre-existing reference; the reconstructed sequence is itself called *an assembly*. This distinction matters throughout, because the pangenome resource used here is built from assemblies, not from read alignments.
2. **Trio sequencing** — sequencing a mother, father and child and reasoning about which allele must have come from which parent.
3. **Statistical phasing** — guessing the phase by asking which arrangement most resembles haplotypes already observed in a large reference population. We formalise exactly this idea as a hidden Markov model in §1.2.

*Why phasing matters here.* Haplotypes are much more identifying than genotypes. A genotype says "you have one copy each of these two rare things." A haplotype says "you have these two rare things *side by side on the same molecule*," which is a far more specific — and therefore far more identifying — statement. The pangenome resource at the centre of this project is explicitly a collection of **haplotypes**.

*Measured in this project:* the pangenome variant file used here carries **45 sample paths**. One of them is **CHM13**, and it needs a correct description because an earlier draft of this primer got it wrong. CHM13 **is human**. It is a human cell line derived from a *complete hydatidiform mole* — an abnormal human conception in which all chromosomes come from the father, so the resulting cell line is effectively derived from a single haploid genome and is homozygous nearly everywhere. That homozygosity is what makes it an outstanding substrate for building a reference assembly, and it is why CHM13 is used as a **reference assembly** rather than treated as a cohort donor whose privacy is at stake. The pipeline drops CHM13 and drops the X chromosome, leaving **44 individuals = 88 haplotypes**. (Recorded in `docs/data.md`, entry for `external/PanMixer/starting_data/pangenome.vcf.gz`, and in `docs/memory.md`.)

#### 1.1.3 Allele frequency, the frequency spectrum, and why rare things are dangerous

First, one term used constantly from here on. A **panel** (or **reference panel**) is a large, published collection of phased haplotypes from many people, used as a statistical reference — the thing you compare a new genome against. Panels are how statistical phasing, imputation, and most of the privacy machinery in this document get their notion of "what human haplotypes normally look like."

The **allele frequency** (**AF**) of an allele is the fraction of chromosome copies in a population that carry it. If 88 haplotypes are examined and 12 carry the alternate allele, its AF is 12/88 ≈ 0.136. The **minor allele frequency** (**MAF**) is the frequency of the *less* common allele at that site, so MAF is always ≤ 0.5. Loose but standard conventions: **common** means MAF > 5%, **rare** means MAF < 1%.

The **allele-frequency spectrum** (**AF spectrum**) is the histogram of how many variants fall into each frequency band — how many are very rare, how many are intermediate, how many are common. It is the single summary statistic population genetics cares about most, which is why distorting it is the headline utility cost of any privacy mechanism that perturbs genotypes. The main scientific consumer of allele frequencies is the **association study**: a test, run variant by variant across a large sample, of whether carrying one allele correlates with a trait or disease.

Why do some variants end up common and others rare? Three forces, each defined here because none of them is basic biology:

1. **Age.** A mutation that arose 100,000 years ago has had time to spread to many descendants. One that arose three generations ago exists in a handful of people.
2. **Genetic drift and population history.** **Drift** is random change in an allele's frequency arising purely from the finite, chance sampling of who actually reproduces — with no selective advantage involved at all. A **bottleneck** is a sharp temporary reduction in population size, which amplifies drift because a small number of individuals supply the next generation's whole gene pool. Expansions do the opposite, preserving newly arisen rare variants that would otherwise be lost.
3. **Selection.** Differential survival or reproduction that systematically pushes an allele's frequency up or down. Alleles that harm survival or fertility are held down; a few advantageous ones are pushed up.

The practical consequence for privacy: the human population has grown explosively in recent history, so most mutations are recent, so the **overwhelming majority of variants are rare**. And *rare* is the same word as *informative*. If I learn you carry an allele present in 50% of people, I have learned almost nothing about *which* person you are. If I learn you carry an allele present in 0.01% of people, I have narrowed the world by a factor of 10,000.

Formally, observing an allele of frequency *f* carries about **−log(f)** units of information, a quantity called **self-information** (or "surprisal"): the rarer the allele, the larger the number. In plain words, −log(f) is *how surprised you should be to see this*, and surprise is exactly what narrows down who someone is. **The units depend on the base of the logarithm: base 2 gives bits, natural log gives nats** (1 nat ≈ 1.443 bits). §1.3.1 below reasons in bits because bits count yes/no distinctions cleanly; from Section 3 onward this document uses **nats**, because PanMixer's code uses the natural logarithm. Whenever a number like "4.4773" appears later, it is in nats.

Hold on to the −log(f) quantity: the privacy score inside PanMixer (the published tool dissected in Section 3) is built out of exactly it.

#### 1.1.4 Recombination, linkage, centiMorgans, and LD blocks

When you make an egg or sperm cell (**meiosis**), your two copies of each chromosome physically line up and swap segments — a **crossover**, or **recombination** event. So the chr21 you pass to your child is not one of your two chr21s; it is a **mosaic**, made by copying from your maternal chr21 for a while, then switching to your paternal chr21, and so on, typically with only one or two switches per chromosome per generation.

Two consequences follow immediately:

- **Nearby variants travel together.** Two loci that are physically close are very unlikely to be separated by a crossover in any one generation, so the allele combination you inherited tends to stay intact across many generations. Two loci far apart get shuffled essentially independently.
- **Every chromosome is a mosaic of ancestral chromosomes.** Go back enough generations and any modern chromosome is a patchwork of segments, each inherited intact from some ancestral haplotype. This is the biological fact that the Li–Stephens model (§1.2.3) turns into an algorithm.

The natural unit for "how likely are these two loci to be separated" is **genetic distance**, measured in **centiMorgans (cM)**. **One centiMorgan is the distance over which one crossover is expected in 1% of meioses** — i.e. a 1% chance per generation that the two loci get separated. Genetic distance is *not* the same thing as physical distance in base pairs: in humans it averages roughly 1 cM per megabase, but the recombination rate varies by more than an order of magnitude along a chromosome, with narrow "hotspots" where crossovers concentrate and long "cold" stretches where they almost never happen. Any model that reasons about how far a haplotype segment extends must use cM, not bp — and certainly not "number of rows in a file." (Flag this now; it becomes a concrete, measured defect in Section 3.)

**Linkage disequilibrium** (**LD**) is the statistical shadow this leaves in data. Two loci are *in linkage equilibrium* if the allele you carry at one tells you nothing about the allele you carry at the other — i.e. the joint frequency factorises, P(A and B) = P(A)·P(B). They are *in linkage disequilibrium* when it does not factorise. The raw measure is

> **D = P(A and B) − P(A)·P(B)**
>
> *In words: how much more often the two alleles occur together than they would if they were independent. D = 0 means no association.*

D depends on the frequencies themselves, which makes raw values incomparable across sites, so LD is normally reported as the squared correlation

> **r² = D² / [ p(A)·p(a)·p(B)·p(b) ]**, where p(a) = 1 − p(A) and p(b) = 1 − p(B)
>
> *In words: treat each locus as a 0/1 variable ("do you carry this allele?") and take the squared correlation between the two. r² = 1 means knowing your allele at one locus tells you your allele at the other exactly; r² = 0 means it tells you nothing.*

Note the letter collision flagged in §1.0: this **r²** has nothing to do with the recombination-rate constant *r* = 1.26 that appears in PanMixer's transition formula in Section 3.

Nearby variants are typically in strong LD because recombination has not had time to break their association. Because LD decays with distance but in a lumpy way, the genome can be carved into **LD blocks**: contiguous runs of variants strongly correlated within the block and much less correlated across block boundaries. Blocks are a modelling convenience — there is no physical wall in the DNA — but they are enormously useful, because they let you treat each block as a semi-independent unit.

> **Caution** — "semi-independent unit" is a statement about **correlation**, not about **geometry**, and the two come apart on real graph data. A record is assigned to a block by its **`POS` alone**, while the sequence that record describes occupies the span `[POS, POS + len(REF) − 1]`, so a long record's span crosses block boundaries. *Measured in this project* on all **340,824** chr21 records: **31,487 records (9.2%)** have their `POS` inside an earlier record's span, **15,540 of 100,757 block-dictionary entries (15.4%)** contain at least one such overlapped record, **78,937 variants (23.2%)** live in an affected entry, and **1,229 of 10,502 (11.7%)** of the entries routed to the HMM are affected. §1.1.8(c) gives the full table and the cross-check behind it. So blocks are **not** disjoint in sequence, and the verdict recorded in `docs/memory.md` is that this is structural rather than an edge case. Both the additive per-block accounting of §A.6 and our own per-block product argument in §4.16.3 assume disjointness; §4.16.8's open question 8 is where that assumption is paid for.

*Measured in this project, on our rebuilt chr21.* Blocks are computed with **PLINK**'s `--blocks` routine. (PLINK is a long-established open-source genetics toolkit; `--blocks` is its confidence-interval LD-block caller, run on an external reference panel rather than on the pangenome itself.) That produces **14,137 LD blocks**. The working **block dictionary** consumed downstream, however, has **100,757 entries**, because every variant that falls inside no called LD block becomes its own one-variant entry. Stating this correctly requires naming the denominator every single time, because there are two very different denominators in play:

| Quantity | Count | Denominator | Percentage |
|---|---|---|---|
| Block-dictionary entries | 100,757 | — | — |
| **Singleton entries** (entries holding exactly one variant) | **89,087** | of 100,757 **entries** | **88.4% of entries** |
| Variants sitting in those singleton entries | 89,087 | of 340,824 **variants** | **26.1% of variants** |
| Entries with **≥ 2 anchors** (the ones routed to the HMM path) | — | of 100,757 **entries** | **10.4% of entries** |
| Variants sitting in those multi-anchor entries | — | of 340,824 **variants** | **71.4% of variants** |

**Never write "88.4% of variants."** 88.4% is a fraction of *entries*. The same singletons are only 26.1% of *variants*. Writing it the other way makes the next line ("10.4% of entries hold 71.4% of variants") arithmetically impossible, which is how the error was caught.

One term used in that table needs defining at first use. An **anchor** is a variant in the pangenome that also appears, matched exactly, in a chosen external callset — a variant for which external population information is available. (In the shipped code the anchor set is the overlap with the PanGenie callset; the published paper instead defines anchors against the 1000 Genomes panel, a discrepancy Section 3 takes apart. On our rebuilt chr21 there are **273,475** anchors.) Blocks with at least two anchors have enough external cross-variant information to be modelled with a haplotype HMM; blocks with fewer are handled by a simpler allele-frequency route.

That skew — a small number of fat blocks holding most of the data, surrounded by a very long tail of singletons — shapes both tools discussed later. *(Block counts, singleton fractions and the 10.4% / 71.4% split are recorded in `docs/memory.md` under the dated entry 2026-09-18 "A-DAT-grch38-repin", with the run logged in `logs/pm_chr21.log`; the anchor count is in `docs/plan.md`.)*

**LD is also the reason that hiding a variant is not the same as hiding its information.** If variant X is in strong LD with variants Y and Z, then deleting X from a release does not remove it: an attacker predicts X from Y and Z. Doing this for scientific purposes is called **imputation**; doing it adversarially is called a **reconstruction attack**. It is the same algorithm either way — and, as §1.2.3 will make clear, it is usually the *same HMM* run in the other direction.

#### 1.1.5 How DNA becomes data: reads, alignment, coverage, assembly

Everything in §1.1.6 and much of Section 3 depends on the mechanics of turning a physical sample into a file, so here they are.

- A **read** is a short DNA fragment whose letters have been determined by a sequencing machine. **Short-read** sequencing, the industry standard, produces reads of roughly 100–300 bp. **Long-read** technologies produce reads of tens of thousands of bp, at higher per-base error but spanning far more variants at once.
- **FASTQ** is the standard plain-text format holding raw reads plus a per-base quality score.
- **Alignment** (also called **mapping**) is finding the location in a reference sequence that a read most closely matches, tolerating a few mismatches and small gaps. An aligner also reports **mapping quality** (**MAPQ**) — its confidence that it put the read in the right place; MAPQ 60 is the usual "as confident as it gets" value.
- **Coverage depth** is how many reads overlap a given base on average. "**30x**" means each base was covered by about 30 reads on average — deep enough for confident genotype calls. This label matters here: one of the two 1000 Genomes releases used in this document is a 30x release, and the older one is roughly 4x, so the newer one is substantially more accurate as well as being on a different coordinate system.
- **Assembly**, as defined in §1.1.2, is reconstructing long contiguous stretches of the genome from the fragments themselves, without leaning on a reference. A **contig** is one such named contiguous sequence; for a finished human reference a contig is usually one whole chromosome. The contig *name* is the key that every coordinate is relative to, which is why a naming-convention mismatch between two files silently breaks joins.

The reason to separate alignment from assembly: alignment can only ever report what the reference already contains a place for, while assembly can produce sequence the reference has never seen. That asymmetry is the whole of the next subsection.

#### 1.1.6 Reference genomes, reference bias, and why pangenomes exist

A **reference genome** is a single, agreed-upon linear DNA sequence that everyone uses as a shared coordinate system. "chr21 position 14,000,000" only means something relative to a named reference **build**. Two builds appear in this project: **GRCh37** (also called hg19) and **GRCh38**. They are different sequences with *different coordinate systems*, so the same biological site carries different position numbers in each, and the same position number denotes different sites. Converting annotations between builds is called **liftover**, and it is lossy and error-prone. Mixing builds silently is one of the classic ways to produce a pipeline that runs cleanly and computes nonsense. (Flag this too; §1.2.1 gives the measured consequence.)

The deeper problem is **reference bias**. The reference is *one* sequence, assembled largely from a small number of donors of predominantly European ancestry. When you sequence a new person, you align their short reads to this reference. A read carrying the reference allele matches perfectly and aligns easily. A read carrying a divergent allele — especially a large insertion that simply *is not present in the reference at all*, so there is nowhere for it to align — mismatches, aligns poorly, aligns to the wrong place, or fails to align entirely. The result is a systematic **under-detection** of non-reference variation, worst for exactly those individuals whose ancestry is most distant from the reference's donor composition. The instrument used to measure human diversity is itself biased against the most diverse samples.

The fix is to stop using one sequence. A **pangenome** is a representation of *many* genomes at once. The natural data structure is a **pangenome graph**:

- **Nodes** are DNA sequence segments (from a single base to many kilobases).
- **Edges** connect nodes observed adjacent in at least one real genome.
- A **path** is a walk through the graph; reading off the node sequences along a path spells out a contiguous DNA sequence.
- Each contributing **haplotype** is stored as a path through the graph. Where all haplotypes agree they share one node; where they differ, the graph opens into a **bubble** — parallel alternative node sequences that diverge and then reconverge. **A bubble is a variant.** Bubbles can be **nested**: a small SNP bubble sitting inside one branch of a large insertion bubble.

Two operations matter for us:

- **Mapping** (or **projecting**) a new genome onto the graph: aligning that individual's reads or assembled sequence to the graph and determining which path best explains them. In a diploid setting the output is a *pair* of paths.
- **Deconstructing** the graph: walking its bubbles and writing each one out as a classical variant record relative to a chosen linear **backbone** path, converting a graph back into a familiar variant file.

The specific resource here is the Human Pangenome Reference Consortium (**HPRC**) v1.0 graph, built with **PGGB** (the Pangenome Graph Builder), and deconstructed to a variant file with `vg deconstruct -a` against the **GRCh38** backbone. (`vg`, short for "variation graph," is the standard open-source toolkit for building, indexing and aligning against pangenome graphs; `vg deconstruct` is its bubble-to-VCF converter, and the `-a` flag tells it to retain nested variants rather than reporting only top-level bubbles.) The resulting file contains 45 sample paths; as described in §1.1.2, the pipeline drops CHM13 and chrX, leaving **44 individuals / 88 haplotypes**, with **340,824 variant records on chr21**. *(Registered in `docs/data.md`.)*

#### 1.1.7 VCF: what it holds, and why a graph can be written as one

**VCF** stands for **Variant Call Format**. It is a plain-text, tab-delimited table, usually compressed with `bgzip` (a block-wise gzip that stays randomly seekable) and indexed with `tabix` (which builds a coordinate index so you can jump straight to a genomic range without reading the whole file). It contains:

- header lines beginning with `##`, declaring the reference build, the contig names and lengths, and the meaning of every field;
- one line per **variant site**, with fixed columns `CHROM` (contig/chromosome name), `POS` (1-based position on the reference), `ID`, `REF` (the reference allele sequence), `ALT` (a comma-separated list of alternate allele sequences), then `QUAL`, `FILTER` and `INFO` — three quality-control and annotation fields that this project does not use and that you can ignore throughout;
- then **one extra column per sample**, holding that sample's call — typically the `GT` (genotype) field: `0|1`, `1|1`, `0/0`, and so on, where the integers index into the list `REF, ALT[0], ALT[1], …`.

The reason a graph can be flattened into a VCF is the bubble correspondence: **one bubble becomes one VCF record.** The backbone path through the bubble supplies `REF`, the alternative traversals supply `ALT`, and *which traversal each haplotype takes* becomes that haplotype's allele index. So the VCF is a lossy but faithful-enough tabular encoding of "which branch did each haplotype take at each branch point." Conversely, given the graph plus a set of VCF rows, a haplotype's **path** and its **allele vector** are two views of the same object. This equivalence is used constantly in this project; it is what lets algorithms designed for VCF-shaped data operate on a graph.

**The variant files this project touches on chr21, and which panel is which.** This table needs one warning attached, because two *different* 1000 Genomes releases appear in this work and conflating them makes Section 3 unreadable:

| File | Build | Records on chr21 | Samples | Haplotypes | Role |
|---|---|---|---|---|---|
| HPRC v1.0 PGGB graph, deconstructed | GRCh38 | 340,824 | 44 (after dropping CHM13) | 88 | the pangenome cohort being protected |
| PanGenie callset (Zenodo 7669083) | GRCh38 | 384,720 | — | — | second source of allele frequencies, and the anchor set the shipped code actually uses |
| **1000 Genomes Phase 3** | **GRCh37** | — | **2,504** | 5,008 | **the panel the shipped pipeline downloads** |
| **1000 Genomes 30x high-coverage** | **GRCh38** | **1,002,753** | **3,202** | 6,404 | **the panel THIS PROJECT substituted** |

A **callset** is simply the set of variants that a particular method reported from a particular dataset. **PanGenie** is a genotyping method; its published callset is used here as a source of allele frequencies, especially for structural variants that SNP-oriented panels do not carry.

The two 1000 Genomes rows are the important ones. **The shipped PanMixer pipeline downloads 1000 Genomes Phase 3 — 2,504 samples, on GRCh37.** The pangenome it is matched against is on GRCh38. This project substituted the **30x GRCh38 panel — 3,202 samples, 1,002,753 chr21 records** — which is what all "our rebuilt chr21" numbers in this document are computed on. Whenever a number derived from a panel is quoted anywhere in this document, the panel that produced it is named. §1.2.1 gives the measured size of the difference.

#### 1.1.8 The correspondence in full: indels, nested bubbles, missing values, cycles, and the worst site on chr21

§1.1.7 stated the correspondence between a graph and a VCF in one sentence — *one bubble becomes one VCF record* — and then moved on to the file inventory. That sentence is true and it is not enough. Every hard question in this project turns out to be a question about the places where the correspondence is *awkward*: where a variant is not a single letter, where one variant sits inside another, where a haplotype has no value at all, where a path loops back on itself, and where a single VCF row is 28 megabytes long. This subsection works through all five, with the measured chr21 numbers attached, because Sections 3 and 4 both assume you have them.

The organising claim, stated once now and justified over the rest of the subsection:

> **Key idea** — The **graph is the primary object** and the **VCF is a projection of it**. The graph holds nodes, edges and paths; a path is a walk, and a walk is the ground truth about what one haplotype's DNA actually is. A VCF record is a *report about a region of that graph*, written in the coordinate system of one chosen backbone path. Reports are lossy. Every oddity below is a place where the report cannot carry something the walk knows — and the honest way to read a VCF built from a graph is as a flattened summary, not as the data.

##### (a) An insertion and a deletion, in a VCF

A SNP is easy to write in a VCF because it is a one-for-one substitution: `REF` is one base, `ALT` is one base. An insertion or a deletion is not one-for-one — one of the two sides is *shorter* — and the VCF format has no way to write an empty string in `REF` or `ALT`. The format's solution is the **padding base**.

> **Definition** — **Padding base (also "anchor base", a different use of the word "anchor" from §1.1.4's).** A shared base immediately to the left of the event, included at the start of **both** `REF` and `ALT` so that neither string is empty. `POS` then points at the padding base, **not** at the first inserted or deleted base. This is the standard VCF convention, not a choice anyone in this project made.

> **Definition** — **Left alignment (also "left normalisation").** In a repetitive stretch the same event can be written at several coordinates: deleting one `A` from `AAAA` gives `AAA` no matter which `A` you delete. The convention is to shift the record as far **left** (toward lower `POS`) as it can go while still describing the same alleles, and to trim any bases that are identical at the right-hand end. Two files describing the same biological event will only join on `(POS, REF, ALT)` if both were normalised the same way — which is part of why §1.2.1 calls that key fragile. (§1.2.1's *measured* 0.30%-versus-75.88% match rates are attributed in `docs/bugs.md` to a reference-build mismatch, not to normalisation; no measurement in `docs/` apportions any of that gap to normalisation, so none is claimed here.)

> **Worked example** — *illustrative only.* The letters and the coordinate in this box are **invented** to show the format convention; nothing in this box is a measurement. Every measured number elsewhere in §1.1.8 is labelled as measured.
>
> Suppose the backbone reads `… T C A G T …` with the `C` at position 1,000.
>
> | | how it is written | reading |
> |---|---|---|
> | **SNP**, `C`→`T` at 1,000 | `POS=1000  REF=C  ALT=T` | one base for one base; no padding needed |
> | **Insertion** of `GG` after the `C` at 1,000 | `POS=1000  REF=C  ALT=CGG` | `REF` is the padding base alone; `ALT` is the padding base plus the inserted sequence. `len(ALT) − len(REF) = +2` |
> | **Deletion** of the `AG` that follows the `C` at 1,000 | `POS=1000  REF=CAG  ALT=C` | `REF` is the padding base plus the deleted sequence; `ALT` is the padding base alone. `len(ALT) − len(REF) = −2` |
>
> Three consequences worth holding on to. First, **the sign of `len(ALT) − len(REF)` is what makes a record an insertion or a deletion**, not anything in the `INFO` column. Second, **`POS` is the padding base**, so the sequence a record actually describes occupies `[POS, POS + len(REF) − 1]` — a span, not a point. That span is the reason "which block is this variant in?" is a harder question than it looks (§1.1.4, and Step 12 of Section 3). Third, in the allele-index view of §1.1.7 **all three rows above look identical**: a haplotype carrying the insertion and a haplotype carrying the SNP both just have the integer `1` in their column. The integer says *which branch*, never *how big the branch was*. Section 3's §A.5 costs an edit at one over the number of called haplotypes **per record**, with no term in base pairs; the representation gives it nothing to weight by, because the allele *lengths* are not in the matrix at all (see §1.1.8(f), which lists what the converter keeps and what it drops).

*Measured in this project, on our rebuilt chr21:* of **340,824** variant records, **65,176** are non-SNP — indels and structural variants with a length change of the kind just described. That is **19.1%** of records (derived: 65,176 / 340,824), so **roughly one record in five** is an event whose `REF` and `ALT` differ in length, which is far too many to treat as an afterthought. (The 340,824 record count is registered in `docs/data.md`; the 65,176 non-SNP count and the 274,458 biallelic-SNP count are both in `docs/memory.md`. Do keep in mind that "the SNPs" is still not one well-defined set downstream: two *other* filters on this chromosome keep **288,020 records (84.5%)** and **238,052 records (69.8%)** respectively, as §1.1.1 sets out.)

##### (b) The same two events in the graph, and why the graph is the primary object

In the graph an insertion and a deletion are **the same shape**: a bubble in which one branch carries sequence that the other branch lacks. There is no padding base, because there is nothing to pad — a branch is a list of nodes, and a branch may simply have fewer nodes than its sibling. There is no left-alignment question either, because a node either is or is not on the walk.

```
        (a) An insertion, as a bubble
                                  ┌──[ GG ]──┐
        ──[ …TC ]───────────────>─┤          ├─>───────────[ AGT… ]──
                                  └────ε─────┘          (ε = no node: the
                                                         branch is a direct edge)
        Haplotype 1 walks  …TC → GG → AGT…      (carries the insertion, allele 1)
        Haplotype 2 walks  …TC ──────> AGT…     (does not, allele 0)


        (b) A deletion, as a bubble  — structurally identical, only the
            backbone's side of the fork has changed
                                  ┌──[ AG ]──┐
        ──[ …TC ]───────────────>─┤          ├─>───────────[ T… ]──
                                  └────ε─────┘
        Haplotype 1 walks  …TC → AG → T…        (has the sequence, allele 0 = REF)
        Haplotype 2 walks  …TC ──────> T…       (deleted it,      allele 1 = ALT)
```

Note what just happened in panel (b): whether the bubble is called "an insertion" or "a deletion" depends **entirely on which branch the chosen backbone path takes**. The graph is symmetric; the VCF is not, because a VCF must nominate a `REF`. That asymmetry is exactly the reference bias of §1.1.6 re-appearing in the file format: the backbone's branch is the one that gets to be "normal".

> **Key idea** — This is why the graph is primary and the VCF is a projection. Deconstruction (§1.1.6) has to choose a backbone, choose a coordinate for each bubble, choose which branch is `REF`, linearise every branch into a flat string, and then throw away everything that does not fit those choices. The walk existed first. Given the graph *and* the VCF rows, a haplotype's path and its allele vector are two views of one object — but given only the VCF, several things the walk knew are gone. The rest of §1.1.8 is those things.

##### (c) Nested variants: a bubble inside a bubble, and the `LV` / `PS` / `AT` tags

§1.1.6 said in one clause that bubbles can be nested. Here is what that means and how much of the chromosome it touches.

> **Definition** — A **snarl** is the graph generalisation of a bubble: a region of the graph delimited by a pair of boundary node-sides, such that every path that enters the region through one boundary leaves through the other. A **bubble** is the simple, acyclic special case — two or more parallel branches that diverge and reconverge. Snarls **nest**: a snarl can lie wholly inside one branch of a larger snarl, and that containment relation forms a tree, the **snarl tree**. The tree's root-most snarls are *top level*; the ones inside them are their *children*.

> **Definition** — A **nested variant** (equivalently a **child variant**) is a variant that sits inside another variant's allele — a bubble within a bubble. `vg deconstruct -a` (§1.1.6) emits **every level** of the snarl tree rather than only the top one, and tags each record in the `INFO` column with:
>
> - **`LV`** — the **level** in the snarl tree. `LV=0` is top level; `LV=1` is a child of a top-level bubble; `LV=2` a grandchild, and so on.
> - **`PS`** — the **parent snarl**: the identifier of the bubble this record sits inside. Absent for `LV=0` records.
> - **`AT`** — the **allele traversals**: for each of the record's own alleles, the list of graph nodes (with orientation) that the corresponding branch walks through. This is the one field that still speaks the graph's own language.

```
        A top-level bubble with a child bubble inside one of its branches

                            ┌─────── ALT branch: [ AAG ] ── [ C ]* ── [ TT ] ───┐
        ──[ backbone… ]──>──┤                                                   ├──>──[ …backbone ]──
                            └─────── REF branch: [ AG ] ────────────────────────┘

                                                        * this node is itself a bubble:
                                                              ┌─[ C ]─┐
                                                          ──>─┤       ├─>──
                                                              └─[ G ]─┘

        PARENT record: LV=0, the whole fork.  Its alleles are the two long branches.
        CHILD  record: LV=1, PS = the parent's id, the little C/G fork.
                       It exists ONLY on haplotypes that took the ALT branch above.
```

*Measured in this project, on our rebuilt chr21* (all from `docs/memory.md` and `docs/terms.md`):

| Metric | Value |
|---|---|
| records at `LV=0` (top level) | **90.8%** |
| records at `LV=1` | **8.2%** |
| records at `LV=2` | **0.9%** |
| records at `LV=3` | **0.1%** |
| records at `LV=4` | **10 records** |
| total records at `LV ≥ 1` (nested) | **31,489** |

And measured a second, independent way — by geometry rather than by tag, taking each record's span as `[POS, POS + len(REF) − 1]` as in §1.1.8(a):

| Metric | Value |
|---|---|
| records whose `POS` lies inside an **earlier record's span** | **31,487 = 9.2%** of 340,824 records |
| …inside a record of `REF` length ≥ 10 kb | **19,211 = 5.6%** |
| …inside a record of `REF` length ≥ 100 kb | **4,608 = 1.4%** |
| block-dictionary entries containing ≥ 1 such overlapped record | **15,540 of 100,757 = 15.4%** |
| variants living in an affected block | **78,937 = 23.2%** of 340,824 variants |
| HMM-path blocks (≥ 2 anchors) affected | **1,229 of 10,502 = 11.7%** |

The two counts — **31,489** records tagged `LV ≥ 1` and **31,487** records whose `POS` falls inside an earlier span — agree to within two records, which is the cross-check that the tag and the geometry are describing the same phenomenon (both measured in this project). The verdict recorded in `docs/memory.md` is worth quoting exactly, because it is the reason this is a Section 1 topic rather than a Section 3 footnote: *"not documentable as an edge case. Roughly a tenth of records and a quarter of variants-by-block are involved, so block disjointness … fails broadly, not locally."*

> **Caution** — A nested variant and its parent are **both rows** in the VCF and therefore **both rows** in the matrix of §1.2.1. The same stretch of DNA is described twice, at two levels of resolution. Any per-row sum — a count of variants, a total privacy score, a total edit cost — therefore counts the nested region more than once. Nothing in the tooling warns you about this, because (see (e)) the tags that would tell you are discarded on the way into the matrix. And the follow-on question — what the released file *literally contains*, and therefore which allele a reconstruction should believe, when a child record's block is rewritten and its parent's block is not — is not answered here. It is answered in Section 3's walkthrough of the stacker and the VCF writer, which is where the released bytes are actually assembled.

##### (d) One parent and one child, written out as VCF rows

Here is a real pair from our chr21, chosen because it is small enough to read. *All values measured in this project* (`docs/memory.md`, 2026-09-25 missingness entry).

The parent is **row 37624**, at **`POS` 13,205,133**, with a `REF` of **13 bp** and **five** `ALT` alleles of **15, 1, 17, 11 and 13 bp**. The child is **row 37625**, at **`POS` 13,205,136** — three bases inside the parent's span.

Schematically — with `REF`/`ALT` sequences abbreviated by their lengths, because the actual strings are not what matters here, and with the sample columns collapsed to the one sample whose call is interesting. Which of the 44 samples that is, which of its two strands carries allele 2, and which allele the *other* strand carries are none of them recorded in `docs/`. The schematic therefore writes the column as *the sample carrying parent allele 2*, puts allele 2 on the right-hand strand, and shows the companion strand as `0` — **those three choices are illustrative, not measured.** What *is* measured is only that exactly one of the 88 haplotypes carries parent allele 2 and has no child call:

```
  #CHROM  POS         ID  REF      ALT                                INFO                   <the sample carrying allele 2>  …
  chr21   13205133    .   <13bp>   <15bp>,<1bp>,<17bp>,<11bp>,<13bp>  LV=0;AT=…              0|2                             …
  chr21   13205136    .   <…>      <…>                                LV=1;PS=<parent>;AT=…  0|.                             …
```

Read the two `GT` fields together. On the parent row this haplotype pair is `0|2`: one strand takes the reference branch, the other takes **allele 2 — the 1 bp allele**. A 13 bp `REF` replaced by a 1 bp `ALT` is, by (a)'s sign rule, a **12 bp deletion**: on that chromosome the region the child bubble lives in has been deleted. So on the child row that strand has nothing to report, and the VCF writes `.`.

*Measured:* across the cohort, **87 of the 88 haplotypes** carry parent allele **0, 1, 3, 4 or 5**, and **all 87 have a child allele**. The **one** haplotype carrying parent allele **2** has `.` at the child. Nothing is broken, nothing is missing in the ordinary sense, and no sequencing failed. The child variant simply **does not exist on that chromosome**.

> **Caution** — **VCF has no way to say "not applicable."** The format offers one symbol, `.`, for every kind of absence: not applicable, not sequenced, not callable, not decided. The information that these are different things is present in the graph and in the `INFO` tags, and it is gone the moment you read only `GT`.

##### (e) Why a haplotype can have no value at a variant — all five causes

This is the passage to read twice, because it is the one a design decision rests on.

> **Question** — *Isn't the whole point of a variant that every haplotype has some value there? How can a donor fail to traverse a parent bubble at all?*

The short answer has two halves, and the first half is reassuring.

**For top-level variants it cannot happen.** Every haplotype in this cohort is stored as a path that runs along the **GRCh38 backbone** (§1.1.6). A top-level bubble is, by definition, a fork *on* that backbone, so every path reaches it and every path takes one of its branches. A `LV=0` variant therefore has an allele for everybody, in principle. *Measured in this project:* top-level records nevertheless show **2.25%** (SNPs) and **3.11%** (non-SNPs) mean missingness, and `docs/memory.md` attributes that residue to **assembly gaps and conflicts, not structure** — causes 3 and 4 below. (That attribution is stated in the ledger without being split between the two; the split is **not measured**, and no number for it exists.)

**For nested variants it can, and routinely does** — for exactly the reason (d) shows. A child bubble sits inside **one particular parent allele**. A haplotype that took a different parent allele never visits that part of the graph.

*Measured in this project,* on the class where it bites hardest:

| Metric | top-level SNP | top-level non-SNP | nested (`LV ≥ 1`) |
|---|---|---|---|
| records | **253,000** | **56,335** | **31,489** |
| mean haplotypes uncalled | **2.25%** | **3.11%** | **22.53%** |
| records > 50% uncalled | **421** | **722** | **3,515** |

Nested variants are roughly **ten times** more often uncalled than top-level ones (derived from the row above). Keep that ratio: it is the real justification for the "anchor" restriction you met in §1.1.4 and will meet again in Section 3.

Now the full enumeration. `docs/memory.md` puts it as "`-1` means **at least five** different things" — *at least*, because the list is an inventory of causes found, not a proof that no sixth exists. They are genuinely different things that a model should arguably treat differently.

1. **Not applicable — the DNA does not exist on that chromosome.** The case in (d): the haplotype took a parent allele that does not contain this child bubble. Nothing failed; there is no allele to report because there is no sequence to report it about. *Measured in this project:* over the **5,735** nested records that have any missing child call, **66.1%** of the missing child calls — among haplotypes that *do* have a parent call — are **perfectly determined by which parent allele the haplotype carries**. "Perfectly determined" is the strong claim here: knowing the parent allele tells you with certainty that the child call will be absent.

2. **Inherited from a missing parent.** If the parent call is itself absent, the child's status cannot be worked out at all — you do not know which branch the haplotype took, so you do not know whether the child bubble was on its route. *Measured in this project:* **12.1%** of that same nested missingness (same denominator, 5,735 nested records).

3. **An assembly gap.** The donor's assembly simply does not cover the region: a contig break, not a biological statement. Because a contig break is a *contiguous stretch* of the chromosome, this cause leaves a distinctive signature — a long **run** of consecutive missing positions in one sample. *Measured in this project,* run lengths of consecutive missing positions on strand 0 across all 44 subjects: **28,052 runs**, median **1**, mean **20.0**, maximum **23,478**. **54.4%** of runs are **isolated single positions** — that is the structural, not-applicable shape of cause 1. At the other end, **349 runs of length ≥ 100** cover **448,283** matrix entries; those are contig breaks.

4. **A `CONFLICT` record.** The sample had **more than one graph path** through the region, and they disagreed, so `vg` declined to pick and wrote the genotype as `.|.`. *Measured in this project:* **1,855** chr21 records (**0.54%**) carry a `CONFLICT` tag, and at such a record a mean of **29.4 of 88** haplotypes are uncalled, against **3.6 of 88** at non-`CONFLICT` sites. The dominant cause of `CONFLICT` itself is **not established** — assembly fragmentation is the most plausible candidate recorded in `docs/memory.md`, but it *"cannot be determined from the VCF alone and has not been verified."* Do not repeat it as fact. What *is* measured is that it is not mainly about repeats: see (g).

5. **A haploid genotype, turned into a missing entry by the converter.** This cause is not in the data at all — it is manufactured downstream. Where the VCF gives **one** allele rather than two (a `GT` of `0` rather than `0|0`), PanMixer's `VCFtoNP.py` writes the allele into strand 0 and **`-1` into strand 1**. *Measured in this project,* over **150,000** chr21 records: **3.135%** of sample `GT` fields are **haploid**, and **5.581%** of the two allele slots are an **explicit `.`**. So a measurable share of the matrix's missingness is a **conversion decision**, not a property of the genome.

*Measured in this project,* the totals those five add up to, in the matrix of §1.2.1: **1,278,653** entries missing on chr21 = **4.26%** of all entries, spread over **52,926** sites = **15.5%** of sites having at least one missing entry.

> **Open question** — **There is no measurement apportioning those 1,278,653 missing entries across the five causes.** The 66.1% and 12.1% figures are shares of *nested* missingness with denominator 5,735 nested records; 3.135% is a share of *`GT` fields*; 349 runs covering 448,283 entries is an absolute count on one strand. The denominators are all different and none of them is "all missing entries." A single apportionment table would be the right thing to have before choosing a policy, and it does not exist in `docs/`. Do not assemble one by arithmetic from the figures above; they are not commensurable.

Why this matters enough to be the centre of the subsection: everything downstream treats the five alike. The edit cost of §A.5 is one over the count of **non-missing** haplotypes, so all five causes inflate it equally. The privacy score of §A.4 is set to **zero** wherever the target's own entry is missing, for all five causes equally. And the choice of what a model should *do* about a missing donor entry — Section 4 calls it the **missing-data policy** — is a live, undecided question precisely because the right answer differs by cause: "not applicable" argues that the donor should be dropped from consideration at that position, whereas an **assembly gap** means the donor almost certainly *has* an allele that we merely do not know, and dropping them penalises poorly-assembled donors for reasons that have nothing to do with genetics.

> **Caution** — The scale of that decision, *measured in this project:* mean missingness at the chain positions PanMixer's HMM actually steps over — anchors inside blocks holding two or more anchors — is **0.26%**; across **every** variant, which is what v1 makes a chain position, it is **4.26%**, with **30,800** records over 25% missing. Worst is the class v1 adds *inside* blocks PanMixer already chains — the non-anchor records of multi-variant blocks — at **11.31%**, with **7,555** over 25% missing (all measured in this project; the per-class table is in `docs/memory.md`, 2026-09-26). **Each of those means is over a different set of positions; none of them is quotable without naming its set.** A policy that is harmless at 0.26% is not automatically harmless at 4.26%, still less at 11.31%. This is recorded as an **open decision** and is deliberately not resolved anywhere in this document.

##### (f) The five causes at a glance

Rows are the causes; the columns are what each one means, whether the VCF can still tell you which one you are looking at, and whether that distinguisher survives the conversion into the matrix. *(Measured shares as cited in (e); the field and survival columns are derived from the converter's source, `external/PanMixer/tools/common/VCFtoNP.py`.)*

| Cause | What it means biologically | Measured share | Distinguishable in the VCF? From which field | Survives conversion to the matrix? |
|---|---|---|---|---|
| **1. Not applicable** | the child bubble is not on this haplotype's route; **the DNA does not exist** on that chromosome | **66.1%** of nested missingness (of 5,735 nested records with any missing child call) | **Yes** — `LV` and `PS` identify the record as a child and name its parent; the parent row's own `GT` then says which branch this haplotype took | **No** — the whole `INFO` column is discarded |
| **2. Inherited** | the parent call was itself absent, so the child's status is unknowable | **12.1%** of the same nested missingness | **Yes** — same two fields, plus the parent row's `GT` being `.` | **No** — same reason |
| **3. Assembly gap** | the donor's assembly does not cover the region; a contig break, not biology | **349 runs of ≥ 100** consecutive missing positions covering **448,283** entries (strand 0, 44 subjects) | **Partly** — by **run length**, not by a tag: a long consecutive run in one sample is the signature | **Not as such** — run length is in principle visible in the matrix, but nothing computes it and no cause label is stored |
| **4. `CONFLICT`** | the sample had multiple graph paths through the region that disagreed, so `vg` wrote `.` | **1,855** records (**0.54%** of chr21); mean **29.4 of 88** haplotypes uncalled there | **Yes** — the `CONFLICT` `INFO` tag, which also **names the samples** | **No** — the whole `INFO` column is discarded |
| **5. Haploid genotype** | nothing is absent biologically; the VCF gave one allele and the converter wrote `-1` into strand 1 | **3.135%** of sample `GT` fields haploid; **5.581%** of allele slots an explicit `.` (of 150,000 chr21 records) | **Yes** — `GT` **arity**: whether the field has a separator at all | **No** — the converter resolves the arity and keeps only the result |

*The converter keeps exactly two things per sample per record:* `fields[1]` (the `POS`) and `fields[9+i].split(':')[0]` (the `GT` subfield). `REF`, `ALT`, the record `ID` and the **entire** `INFO` column — `LV`, `PS`, `AT`, `CONFLICT` — are dropped. (Read from the released converter source, as recorded in `docs/bugs.md`; this is a reading of code, not a measurement.) So by the time the data reaches the matrix, **all five causes look identical: the integer `-1`.**

> **Key idea** — Every one of the five is recoverable at the moment of conversion, and none is recoverable afterwards. `docs/plan.md` records the resulting recommendation in five words: *"Cheap now, unrecoverable later."* The concrete form is to emit an auxiliary **reason array** alongside the allele matrix during conversion, one label per missing entry, so that a policy can later be chosen per cause instead of once for all five.

##### (g) Cycles: a walk that leaves a node, returns to it, and leaves again

> **Question** — *Can a pangenome graph contain loops? Isn't it a directed graph?*

It is a directed graph, and it can still contain loops. Those are not contradictory: "directed" constrains which way you may traverse an edge, not whether the edge set contains a cycle. A directed graph with no cycles is a **directed acyclic graph (DAG)**; a directed graph that is *not* a DAG contains at least one cycle. A walk may then enter a node, continue forward, and arrive back at a node it has already used — and then leave it again along a different edge than last time. `docs/` makes no general claim about whether **PGGB** (§1.1.6) graphs are DAGs, so none is made here. What it does record is the specific case, *measured in this project:* our HPRC v1.0 chr21 graph contains at least one traversal that revisits nodes, so **this** graph is not a DAG. That traversal is the subject of the rest of (g).

Nothing about a path "remembers" that it has been somewhere before. A path is just a **sequence of steps**, and a node may appear at more than one index in that sequence. The walk `n1 → n2 → n3 → n2 → n4` is perfectly well formed; it visits `n2` twice and spells out `n2`'s sequence twice when you read the DNA off it.

**This is how a tandem duplication is encoded.** A tandem duplication is a stretch of DNA present twice in a row. In a linear reference you would have to write the second copy out as inserted sequence. In the graph you do not write it at all: the walk simply goes around the region again. The sequence is stored **once**, in the nodes; the *repetition* lives in the walk, not in the storage. That is the compactness a graph buys you.

Here is the real one, *measured in this project* (`docs/memory.md`, 2026-09-22): at the chr21 site examined in (h), **allele 58** of the top-level bubble walks **10,008 node steps** over only **5,014 distinct nodes**. **4,994** of those nodes appear **exactly twice** and **20** appear **once**. At around step **5,004** the walk reaches node `>102277684`, which sits beside the bubble's *closing* boundary, and jumps back to `>102270112`, which sits beside its *opening* boundary — and then walks the whole region a second time. That is a **whole-region tandem duplication**, read directly off the traversal.

> **Question** — *Then why can a duplication not just be written as an insertion variant, with the second copy as inserted sequence?*

Two reasons, and the second is the deep one.

1. **The copies are not identical.** Of allele 58's nodes, 4,994 appear twice but **20 appear once** — which means there are roughly **10 interior positions where the second copy diverges from the first** (measured in this project). So the region is not "this sequence, times two". It is "this sequence, then a slightly different version of this sequence". You cannot collapse it into a clean copy-number count, because there is no single thing being counted.
2. **VCF positions are reference coordinates, and inserted sequence has none.** Suppose you did write the duplication as one big insertion. Now you want to record that the second copy differs from the first at those ~10 interior positions. There is **nowhere to put those records**: a VCF row must carry a `POS`, `POS` is a coordinate on the backbone, and the second copy is not on the backbone. The format offers no coordinate for "inside the inserted sequence." **This is precisely the limitation pangenome graphs exist to remove** — in the graph, an interior variant of the duplicated region is just another bubble on the walk, and it needs no backbone coordinate at all.

That second point also explains something that otherwise looks like a bug: **why the VCF line for this site is 28 megabytes long.** The graph stores the repeat compactly, duplicating no sequence. The ~28 MB is purely a **VCF export artefact** — a VCF allele is a **flat string**, so a cyclic traversal has to be **linearised**, writing every base of both copies out in full, for every allele. *Nothing about the graph is 28 MB.*

*Measured in this project,* how common this is, and one important limit on that measurement:

| Metric | Value |
|---|---|
| records whose own `AT` field contains a **revisiting** traversal | **37 = 0.011%** of 340,824 |
| records carrying a `CONFLICT` tag | **1,855 = 0.54%** |
| records with **both** | **11** |
| revisiting but **not** flagged `CONFLICT` | **26** |
| flagged `CONFLICT` but **no** revisiting traversal | **1,844** |
| traversals at the (h) site using **reverse** orientation (an inversion) | **0 — none** |

> **Caution** — That 37-record count is a **floor, not an estimate**, and the reason is worth understanding. A record's `AT` lists the traversals of **its own** bubble. It therefore cannot reveal that a sample passed through **this** bubble twice — that fact is only visible at the **enclosing parent**. So the cross-tabulation **understates** the link for child records sitting inside a duplicated region. Concretely: the 257 kb parent at (h) **has** a revisiting allele and is **not** flagged `CONFLICT`, because at the parent level the sample `HG01891` has exactly **one** allele — number 58 — whose walk happens to loop. The conflict only materialises at the **child** records, where that single parent allele corresponds to **two passes that may disagree**, and there is no way to say in the child's `GT` which pass you mean. *Measured:* `CONFLICT=HG01891` duly appears at **`POS` 14,577,263**. The mechanism is real; **1,844 conflicts with no revisit** shows it is not the main source of conflicts genome-wide.

##### (h) A site where every haplotype differs: `grch38#chr21:14569980`

> **Question** — *Suppose a single position is different for every single haplotype. What would that look like in the VCF versus in the graph?*

It looks like this, and it is real. *Every value in the table is measured in this project* (`docs/memory.md`, 2026-09-22; `docs/terms.md` entry `hypervariable_site`).

The site is **`grch38#chr21:14569980`**, graph bubble **`>102270111>102277685`**, **`LV=0` — top level**, VCF **row 51791**.

| Property | Value |
|---|---|
| alleles **declared** on the record | **90** (89 `ALT`) |
| distinct alleles actually **carried by the cohort** | **88, across 88 haplotypes — every haplotype unique** |
| `REF` length | **257,439 bp** |
| `ALT` lengths | **257,124 – 514,768 bp** |
| size of the **single VCF line** | **27,983,949 bytes (~28 MB)** |
| `AC` for allele 1 | **0** — an orphan left behind by `view -s ^chm13` **without** `--trim-alt-alleles`: the sample was removed, its allele was not |
| **nested child records** (`LV ≥ 1`, `PS` = this bubble) | **2,464** |
| total matrix rows for this one region | **2,485 = 0.7%** of chr21's 340,824 rows |
| distinct graph nodes inside the bubble | **7,575**; traversals span **4,971 – 10,008** nodes |
| allele 58, structurally | a **whole-region tandem duplication** — 10,008 steps over 5,014 distinct nodes, ≈ 4,994 revisits (see (g)) |
| inversions | **none** — 0 traversals use reverse orientation |
| block-dictionary entries the region **spans** | **370** overlapping entries |

*Measured in this project,* the class is small but not unique: **7** chr21 sites have **≥ 88** declared alleles, **155** have **≥ 44**, and **618** have **≥ 20**.

**Now the question that this site is really about.** How can **one** VCF record describe a region in which thousands of interior differences exist? Look at the last two rows of the table against the row above them. There are **2,464 nested child records** here, so the nested decomposition **did** happen — the interior differences *are* broken out into their own rows, one bubble at a time, and they are exactly the `LV ≥ 1` records of (c). The top-level record is emitted **in addition to** all 2,464 of them, and there are two reasons it has to be:

1. **Some alleles cannot be expressed as a combination of the children.** Allele 58 traverses the region **twice** (see (g)). No set of independent per-site substitutions can encode "this whole region is duplicated" — a duplication is a statement about the *shape of the walk*, and the child records only describe *which branch at each fork*. If you kept only the children you would lose allele 58 entirely.
2. **VCF has no syntax for "this allele = this combination of child variants."** A VCF allele is a flat string and nothing else. There is no way to write "allele 58 is: take the parent's REF branch, go round twice, and at such-and-such a child bubble take allele 2." So a top-level traversal has to be **spelled out in full**, base by base — all 257,124 to 514,768 of them, ninety times over. **Hence 28 MB.**

So the answer to "one record or thousands?" is **both, simultaneously, by design**: one top-level record that can express walk-shape but costs 28 MB to write, plus 2,464 child records that are cheap but cannot express walk-shape. They describe the same DNA at two resolutions, and — per (c)'s caution — every one of the 2,485 is a separate row in the matrix.

##### (i) Why that site is at once the worst case for privacy and an opportunity

**Worst case for privacy, by construction.** Every one of the 88 haplotypes carries a **distinct** allele here. That is not a variant, it is a **fingerprint**: learning a person's allele at this one site distinguishes them from all 87 others in the cohort outright. Put it in the self-information terms of §1.1.3. *Measured in this project:* the target `HG00438` carries alleles **89** and **88** at this record, each with a frequency of **`f_v` = 1.540 × 10⁻⁴**. That figure is the value `docs/memory.md` gives for an allele carried by a single haplotype at a **novel** site — one matched in no external panel — where the denominator is padded to `NOVEL_DENOMINATOR = 2 × (N_1000G + N_pangenome)` = 2 × (3,202 + 44) = **6,492** (derived, and it carries the panel: our pinned **1000 Genomes 30x GRCh38 panel** of 3,202 samples plus the 44 pangenome individuals; under the shipped Phase 3 panel of 2,504 samples the denominator would be a different number). 1 / 6,492 = 1.540 × 10⁻⁴ (derived), and −log(1.540 × 10⁻⁴) = **8.78 nats** of self-information per allele (derived). For scale, *measured in this project* across the sites where `HG00438` carries a non-reference allele, the median self-information is **0.63 nats** where the variant matched the 1000 Genomes 30x panel and **0.43 nats** where it matched the PanGenie callset instead. So this one record carries **roughly fourteen times** the self-information of a typical 1000-Genomes-matched non-reference call (derived: 8.78 / 0.63 = 13.9; against the PanGenie-matched median of 0.43 nats the ratio is 20.4, also derived). On the arithmetic of §1.1.3 and §1.3.1 it is, in `docs/memory.md`'s own words, the most identifying record on the chromosome.

**And yet** — this is the part to carry into Section 3 — the score PanMixer computes for the block containing it **never reads this record at all**. *Measured in this project:* the record is **not an anchor** (§1.1.4: it is not in the PanGenie callset), and the block it lands in — **block 1639, 92 records of which 82 are anchors** — therefore takes the haplotype-model branch, where the score is computed over the block's **anchors only**. The 8.78 nats are computed and then never used.

> **Caution** — State that finding in its correct form, which is narrower than it first sounds. It is **not** a claim that the published results are wrong; we have reproduced none of the paper's figures. It is a claim about **what the metric measures**: the per-block privacy score is the self-information of the block's **anchor alleles**, not of the block. A reader who takes the name at face value will over-read it. §A.4.3 gives this its own treatment with the measured comparison. **And a retraction belongs here.** An earlier version of this primer explained this site's zero score by a `−log(0)` code path: a frequency of exactly 0 makes `−log f_v` infinite, and the released code detects the overflow and clamps the value to **0** — the minimum rather than the maximum — at the line that carries a literal `#FIX ME` comment. **That explanation is withdrawn.** *Measured in this project:* across all **319,092** chr21 sites where `HG00438` carries a called allele, **zero** have `f_v = 0`, so that branch never fires on our data. The original reasoning failed because it assumed a very rare allele would come back from the frequency table with a frequency of exactly 0. It cannot here, for two reasons: the table counts **the target's own allele** in every one of its three branches, and at a **novel** site it floors the frequency by padding the denominator to `NOVEL_DENOMINATOR` = 6,492, so an allele carried by a single haplotype comes out at **1.540 × 10⁻⁴** rather than 0 (measured in this project; the division is derived). The clamp is therefore a **latent** defect with a panel-dependent trigger, not the mechanism operating at this site; the mechanism here is non-anchor exclusion (§A.4.3).

**Simultaneously an opportunity — to evaluate, not yet a claim.** Three specific reasons, each with its own caveat, and **nothing here is measured about our mechanism, because our mechanism does not exist yet:**

- A **path-level** guarantee of the kind Section 4 builds is stated in terms of how distinguishable two whole releases are, and does not weaken as an allele gets rarer. This record, by contrast, contributes nothing to the score PanMixer computes for its block — *measured in this project*, by the mechanism established just above: the record is **not an anchor**, and its block takes the anchor-only HMM branch. (Do **not** state this as "the frequency-based score collapses where rarity is extreme": `docs/memory.md` Correction 7 retracts exactly that explanation for this site, having measured that the record receives a perfectly ordinary `f_v`. The zero contribution is an anchor-set effect, not a frequency effect.) If that contrast holds up, this site class is where the two approaches should differ most visibly.
- The **segmentation** could be better. This one region spans **370** overlapping block-dictionary entries (measured in this project), because block occupancy is decided by a record's `POS` alone (see Step 12: every variant is placed by a binary search on its position) while its span is 257,439 bp. A **snarl-tree** segmentation would make it **one** unit instead of 370, with units disjoint by construction and every boundary a real graph cut point — so switching donors at a boundary would be **graph-valid by construction**. The standing objection is that snarls are not LD blocks; `docs/memory.md` records that this objection is weak for a copying model, because the linkage is carried by the donor panel and the transition probabilities rather than by where the boundaries fall. This is **not decided** — what the model's positions should be is an explicitly reserved question.
- **Cost should probably be sensitive to scale.** Rewriting this record changes up to half a megabase; §A.5's per-record cost charges the same for it as for a one-base SNP. *Measured in this project* at capacity 0.1 on chr21 for `HG00438`, this record **was** rewritten on both strands — alleles **[89, 88] → [35, 64]** — and the largest single allele rewrite in that run was **257,485 bp**. Whether a length-weighted cost is the right correction is an open design question, not a decision.

Three caveats stay attached to all of that: whether a sampled release can even contain a **valid** traversal here (an independent per-variant sampler would happily produce a combination of child alleles that no real walk realises); whether cost should be length-weighted given that one unit can span 257 kb; and how to avoid **double-counting** the top-level record and its 2,464 children, given (c). None is settled.

> **Note** — This site recurs throughout the document, so it is worth a bookmark: `grch38#chr21:14569980`, row 51791, 90 declared alleles, 88 distinct, 257,439 bp of `REF`, 28 MB of text, 2,464 children, 370 blocks, one tandem-duplicated allele. Whenever a later section says "the hypervariable site," this is it.

#### In plain words (biology)

You have two copies of every autosome, one from each parent; a haploid human genome is about 3.1 billion base pairs, so a diploid cell holds about twice that, spread across 22 numbered autosomes plus X or Y — 24 distinct chromosome sequences in the species, 23 pairs in a person. At millions of positions people carry different DNA letters. Those positions are variants, the alternatives are alleles, and we label them 0 for "same as the reference" and 1, 2, … for the alternatives; a variant with just one alternative is biallelic, with more it is multi-allelic. A genotype tells you *which two* alleles you carry at a position but not which physical chromosome copy each sits on; a haplotype tells you the whole ordered run of alleles along one copy, and working out which is which is called phasing. Haplotypes are far more identifying than genotypes. Most variants are rare, and rare variants are exactly the ones that pick a person out of a crowd, because the information in seeing an allele of frequency f is about −log f — bits if you use log base 2, nats if you use natural log. Because chromosomes are reshuffled only once or twice per generation, nearby variants travel together; that correlation is linkage disequilibrium, it lets us chop the genome into blocks, and it means deleting one variant does not really hide it because its neighbours predict it. Genetic distance for this purpose is measured in centiMorgans, not base pairs. Sequencing machines emit short reads, which are either aligned onto a reference or assembled into long contiguous sequence; alignment can only find what the reference has a slot for, which is why a single linear reference systematically under-detects variation in people unlike its donors. A pangenome graph fixes that by storing many genomes at once as paths through one shared graph, where each branch point (bubble) is a variant. A VCF file is a big table with one row per variant and one column per sample holding small integers, and because each bubble maps to one row, the graph and the table are two views of the same data. On our rebuilt chr21 that table is 340,824 rows by 88 haplotypes, and it is carved into 100,757 block-dictionary entries — of which 88.4% *of entries* are singletons, which is only about 26% *of variants*. That correspondence is exact for a SNP and awkward for everything else: an insertion or a deletion has to be written with a shared padding base so neither `REF` nor `ALT` is empty, so a record really describes a *span* rather than a point, and in the graph both are just a bubble where one branch carries sequence the other lacks — which is why the graph is the primary object and the VCF a projection of it. Bubbles also nest, one inside a branch of another, and on our chr21 **31,489 records are nested** (measured in this project), **9.2%** of records (derived), which is why a haplotype can legitimately have **no value at all** at a variant — at least five different things produce that absence, from "the DNA does not exist on that chromosome" through an assembly gap to a `CONFLICT` and a merely haploid genotype, and the conversion into the matrix reduces all five to the same `-1`. The graph can also contain **cycles**, so a walk may leave a node, come back and leave again, which is how a tandem duplication is stored compactly — and because the two copies are not identical, the repeat cannot be collapsed into a clean count. At the extreme sits one chr21 site where all **88** haplotypes carry **88 distinct alleles** in a single **28 MB** VCF line with **2,464** nested children (measured in this project): the most identifying record on the chromosome, and the case every later section is tested against. §1.1.8 works all of this through.

---

### 1.2 The computer science you need

#### 1.2.1 A cohort as a matrix

Once you accept the VCF↔graph correspondence, the data model is embarrassingly simple: **a matrix of small non-negative integers**, rows indexed by variant site, columns indexed by haplotype, each entry the allele index (0 = reference, k ≥ 1 = the k-th alternate). For our chr21 pangenome that is a 340,824 × 88 matrix. For the 1000 Genomes 30x GRCh38 panel on chr21 it is 1,002,753 × 6,404. Whole-genome, at cohort scales that are now routine, it is tens of millions of rows by tens of thousands of columns.

That is not large by modern standards *if you only ever stream it* — but the operations we care about are not streaming. We want to ask "what happens inside this block," "which panel row corresponds to this graph row," "what is the frequency of this allele." So the data is kept sorted by coordinate, compressed in seekable blocks (`bgzip`), and indexed (`tabix`) for coordinate-range random access; and pipelines build in-memory **dictionaries keyed by the tuple (POS, REF, ALT)** to match a row in one file to the corresponding row in another.

Why that key, and why it is fragile: `POS` alone is not enough because two different variants can start at the same position, so the key must also pin down the exact `REF` and `ALT` allele **strings**, which must match character for character. And `POS` is only meaningful relative to a build. If one file is GRCh37 and the other GRCh38, the `POS` values name different sites and the dictionary silently comes up nearly empty.

*Measured in this project, on chr21, using PanMixer's own matching code (`get_mappings`) against the identical set of 340,824 pangenome records, varying nothing but the panel's build:*

| Metric | 1000 Genomes Phase 3 (wired in the shipped pipeline, GRCh37) | 1000 Genomes 30x (substituted here, GRCh38) |
|---|---|---|
| exact (POS, REF, ALT) matches | **1,025** | **258,610** |
| relaxed (POS, REF) matches | 2,176 | 268,851 |
| match rate out of 340,824 pangenome records | 0.30% | 75.88% |

A **252-fold** gap, produced by the reference build alone. Nothing crashes; no error is raised; the join just quietly returns almost nothing, and every downstream artifact built on it is built on 1,025 spurious rows instead of 258,610 real ones. *(Recorded in `docs/memory.md` under 2026-09-18, in `docs/bugs.md`, and in `docs/data.md`; an independent `bcftools` comparison reproduced the same counts.)* The residual ~24% non-match on the correct build is expected — the pangenome carries graph-only and structural variants absent from a SNV/indel/SV panel — but it has not been characterised and should not be assumed harmless.

#### 1.2.2 Hidden Markov models, from scratch

This subsection is long and deliberately slow. The hidden Markov model is the single
piece of machinery both PanMixer (Section 3) and our own mechanism (Section 4) are
built on, so it is worth doing properly rather than quickly. Everything here is
self-contained: if you can multiply probabilities together, you can follow it.

##### (a) First, a Markov chain — nothing hidden yet

A **Markov chain** is a sequence of random states where the next state depends only
on the *current* state, and not on how you got there. That "only the present
matters" rule is the **Markov property**.

Write the states as 1, 2, …, K. The chain needs two things:

- An **initial distribution** — call it pi(i) — the probability the chain starts in
  state i. These sum to 1 over all i.
- A **transition matrix** A, where A(i, j) = P(next state is j | current state is i).
  Each *row* sums to 1: from state i you must go somewhere.

That is the whole definition. If you know pi and A, you can generate sequences: draw
the first state from pi, then repeatedly draw the next state from the row of A
corresponding to where you are.

##### (b) Now hide it

In a **hidden Markov model (HMM)** you never see the states. You see something
*generated* from each state, called an **emission** or an **observation**. So there
are two parallel sequences:

```
hidden:    Z_1  ->  Z_2  ->  Z_3  ->  ...  ->  Z_T      (states; you never see these)
             |        |        |                 |
observed:   X_1      X_2      X_3               X_T      (emissions; these you see)
```

The arrows across the top are transitions (governed by A). The arrows going down are
emissions. Two independence assumptions define the model, and both matter later:

1. **Markov property (horizontal):** Z_t depends only on Z_{t-1}.
2. **Conditional independence of emissions (vertical):** X_t depends only on Z_t.
   Given the state at position t, the observation at t tells you nothing extra about
   any other position.

##### (c) The four ingredients, in this project's terms

Here is the translation that makes the abstraction concrete. Throughout this project:

| HMM ingredient | Symbol | What it is here |
|---|---|---|
| **Number of positions** | T | Sites (or blocks) along a chromosome, left to right |
| **States** | 1…K | The K haplotypes in the **copying panel** — "which panel haplotype am I currently copying from?" See the Definition below: in this project the copying panel is the *cohort's own* haplotypes, not an external reference panel |
| **Initial distribution** | pi(i) | Which panel haplotype the copying starts on. Usually uniform, 1/K each |
| **Transitions** | A(i, j) | The chance of *switching* which haplotype you copy from between consecutive positions. This is how **recombination** enters the model |
| **Emissions** | e_t(i) | Given that you are copying haplotype i at position t, how likely is the allele you actually observed. This is how **mutation and sequencing error** enter |

The mental image — the **Li-Stephens copying model** from §1.2.3 — is that a new
haplotype is assembled by copying along one panel haplotype for a while, occasionally
jumping to another, with occasional single-letter mismatches. The hidden state is
"who am I copying right now"; the observation is "what letter came out".

##### (d) The probability of one specific path

Before any algorithm, be clear on what the model actually assigns probabilities to.
A **path** is one complete assignment of hidden states, z_1, z_2, …, z_T. The model
gives the joint probability of a path *together with* the observations:

> **P(path and observations)**
> = pi(z_1) · e_1(z_1) · A(z_1, z_2) · e_2(z_2) · A(z_2, z_3) · e_3(z_3) · … · A(z_{T-1}, z_T) · e_T(z_T)

*In words: the chance of starting where you started, times the chance that state
emitted what you saw, times — for each step afterwards — the chance of moving where
you moved and the chance that new state emitted what you saw.* Read left to right,
it is exactly the generative story, one factor per arrow in the diagram above.

##### (e) The question we want answered, and why brute force fails

Two things are wanted. **First**, the total probability of the observations,
P(X_1…X_T) — this is the "how well does this model explain what I saw" number, and in
Section 3 it becomes a privacy score. **Second**, a *random path drawn in correct
proportion to its probability* — this is what actually produces a released haplotype
in Section 4.

The definition of the first is a sum over every possible path:

> P(observations) = sum over all paths z of P(path z and observations)

The obvious approach is to enumerate. In the toy example below, K = 3 and T = 3, so
there are 3³ = **27** paths — easy. But the count is K^T, and on our real chr21 data
K = 88 haplotypes and T ≈ 100,757 blocks, giving 88^100757 paths. That is not a large
number; it is an absurd one, vastly beyond the count of atoms in the observable
universe. Enumeration is not slow, it is impossible.

The **forward algorithm** computes the identical sum in time proportional to T·K²,
and under the structure we use, T·K. This is the classic dynamic-programming trade:
instead of walking every path separately, notice that all paths passing through the
same state at the same position share their entire future, so their pasts can be
added up *once* and reused.

##### (f) The forward algorithm, derived

Define the **forward variable**:

> **alpha_t(i) = P(the first t observations, AND the state at position t is i)**

Note carefully what this is: a *joint* probability of "what I have seen so far" and
"where I am now". It is not a conditional probability and it does not sum to 1 across
states — it sums to the probability of the observations so far, which shrinks as t
grows.

**Base case (t = 1).** To be in state i at position 1 and have seen X_1:

> alpha_1(i) = pi(i) · e_1(i)

**Recursive step.** To be in state j at position t having seen the first t
observations, you must have been in *some* state i at position t−1, moved to j, and
then emitted what was observed:

> alpha_t(j) = [ sum over i of alpha_{t-1}(i) · A(i, j) ] · e_t(j)

Why this is valid, spelled out: the bracket is the total probability of arriving in
state j at position t having accounted for everything up to t−1 — summed over every
way of arriving. It is legitimate to multiply that single bracket by e_t(j) only
because of the vertical independence assumption: given that we are in state j, the
emission at t does not care which route we took to get there. And it is legitimate to
reuse alpha_{t-1}(i) without re-deriving history only because of the Markov property:
given we are in state i at t−1, the past is irrelevant to the future.

**Termination.** Once you reach the end, sum over states to release the constraint of
being anywhere in particular:

> P(observations) = sum over i of alpha_T(i)

##### (g) The worked micro-example, with every multiplication shown

A panel of **K = 3** haplotypes at **T = 3** biallelic sites (alleles written 0 and 1):

| | site 1 | site 2 | site 3 |
|---|---|---|---|
| **h1** | 0 | 1 | 0 |
| **h2** | 0 | 1 | 1 |
| **h3** | 1 | 0 | 1 |

We observe a **new** haplotype: **`0, 0, 1`**. Read down the columns and check: no
single panel haplotype matches it at all three sites. h1 matches at sites 1 and 2 but
not 3; h2 matches at 1 and 3 but not 2; h3 matches only at 3. The model will have to
explain this new haplotype as a *mixture*.

Parameters:
- **Initial:** uniform, pi(i) = 1/3 for each of the three haplotypes.
- **Transitions:** stay in the same state with probability **a = 0.90**; switch to
  each of the two others with probability **b = 0.05**. Check the row sums to 1:
  0.90 + 0.05 + 0.05 = 1.00. ✓
- **Emissions:** if the state's allele equals the observed allele, probability
  **0.99** (a "match"); otherwise **0.01** (a "mismatch"). The 0.01 is the model's
  allowance for mutation and sequencing error — without it, any mismatch would make a
  path strictly impossible rather than merely unlikely.

First, write down the emission values at each site, since they are used repeatedly.
Compare each panel row against the observed allele for that column:

| | observed | e(h1) | e(h2) | e(h3) |
|---|---|---|---|---|
| **site 1** | 0 | h1 has 0 → match → **0.99** | h2 has 0 → match → **0.99** | h3 has 1 → mismatch → **0.01** |
| **site 2** | 0 | h1 has 1 → mismatch → **0.01** | h2 has 1 → mismatch → **0.01** | h3 has 0 → match → **0.99** |
| **site 3** | 1 | h1 has 0 → mismatch → **0.01** | h2 has 1 → match → **0.99** | h3 has 1 → match → **0.99** |

**Step 1 — alpha_1 (base case).** alpha_1(i) = pi(i) · e_1(i):

- alpha_1(h1) = (1/3) × 0.99 = **0.330000**
- alpha_1(h2) = (1/3) × 0.99 = **0.330000**
- alpha_1(h3) = (1/3) × 0.01 = **0.003333**
- Running total S_1 = 0.330000 + 0.330000 + 0.003333 = **0.663333**

*Interpretation: after one site, h1 and h2 are equally good and h3 is about 99 times
worse, exactly as the single observed allele implies.*

**Step 2 — alpha_2.** Do one entry the long way first, to see the sum over i
explicitly. For j = h3:

> bracket = alpha_1(h1)·A(h1,h3) + alpha_1(h2)·A(h2,h3) + alpha_1(h3)·A(h3,h3)
> = 0.330000 × 0.05 + 0.330000 × 0.05 + 0.003333 × 0.90
> = 0.016500 + 0.016500 + 0.003000
> = 0.036000
>
> alpha_2(h3) = 0.036000 × e_2(h3) = 0.036000 × 0.99 = **0.035640**

Note what just happened: h3 was a *terrible* candidate after site 1 (0.003333), but
site 2 is the one site where h3 matches and the others do not. Most of h3's new mass
came not from h3's own past but from h1 and h2 **switching into** it — 0.0165 + 0.0165
of the 0.036 total. That is recombination doing its job.

The other two entries, same method:

- bracket for h1 = 0.330000 × 0.90 + 0.330000 × 0.05 + 0.003333 × 0.05 = 0.297000 + 0.016500 + 0.000167 = 0.313667 → × 0.01 = **0.003137**
- bracket for h2 = the same by symmetry = 0.313667 → × 0.01 = **0.003137**
- S_2 = 0.003137 + 0.003137 + 0.035640 = **0.041913**

**Step 3 — alpha_3.**

- bracket for h1 = 0.003137 × 0.90 + 0.003137 × 0.05 + 0.035640 × 0.05 = 0.002823 + 0.000157 + 0.001782 = 0.004762 → × 0.01 = **0.000048**
- bracket for h2 = 0.003137 × 0.05 + 0.003137 × 0.90 + 0.035640 × 0.05 = 0.004762 → × 0.99 = **0.004714**
- bracket for h3 = 0.003137 × 0.05 + 0.003137 × 0.05 + 0.035640 × 0.90 = 0.000157 + 0.000157 + 0.032076 = 0.032390 → × 0.99 = **0.032066**

**Termination.**

> P(observations) = 0.000048 + 0.004714 + 0.032066 = **0.03682760**

##### (h) Checking it against brute force

Because this example is tiny, the forward algorithm's shortcut can be checked against
the definition directly. Enumerating all 27 paths, computing each one's joint
probability with the master formula from (d), and adding them gives
**0.0368276033** — identical to the forward algorithm's answer to every digit shown.
*(Verified in this project; the enumeration and the recursion were run side by side.)*

The enumeration also shows *where* that probability lives, which the single summary
number hides:

| path | joint probability | share of the total |
|---|---|---|
| h1 → h3 → h3 | 0.014554 | **39.52%** |
| h2 → h3 → h3 | 0.014554 | **39.52%** |
| h3 → h3 → h3 | 0.002646 | 7.19% |
| h2 → h2 → h2 | 0.002646 | 7.19% |
| h1 → h3 → h2 | 0.000809 | 2.20% |
| h2 → h3 → h2 | 0.000809 | 2.20% |
| *(the other 21 paths combined)* | — | 2.20% |

Two paths carry 79% of the probability, and both have the same shape: start on a
haplotype that matches site 1, switch to h3, stay there. The model has discovered the
mosaic structure by itself.

##### (i) The stay/switch shortcut: why the cost is T·K and not T·K²

Computing each alpha_t(j) as written requires summing over all K predecessors, and
there are K values of j, so each step costs K² multiplications — T·K² overall. With
K = 88 and T = 100,757 that is about 780 million operations. Survivable, but there is
a much better way when the transition matrix has the **stay-or-switch** form used
throughout this project: A(i, j) = a if j = i, and b if j ≠ i.

Look at the bracket again, and add and subtract the diagonal term:

> sum over i of alpha_{t-1}(i)·A(i, j)
> = [ b × (sum of ALL alpha_{t-1}) ] + (a − b) × alpha_{t-1}(j)

*Derivation in one line: pretend every state contributes b — that is b times the grand
total. Exactly one state, j itself, should have contributed a rather than b, so patch
that single entry by adding (a − b)·alpha_{t-1}(j).*

Confirm it on the number we computed the long way, alpha_2(h3), with S_1 = 0.663333:

> b·S_1 + (a − b)·alpha_1(h3) = 0.05 × 0.663333 + 0.85 × 0.003333
> = 0.033167 + 0.002833 = **0.036000** ✓ — the same bracket as before.

The grand total is computed **once per position**, then each state needs one multiply
and one add. The cost drops from T·K² to **T·K**: about 8.9 million operations instead
of 780 million on chr21. This is why these methods run genome-wide, and it is the one
piece of PanMixer's HMM code that is worth inheriting (Section 4).

##### (j) A practical problem: the numbers vanish

Notice the alphas shrinking: 0.66 → 0.042 → 0.037, and that is after three sites. Each
position multiplies in another emission probability, so the values decay roughly
geometrically. After a few hundred sites they fall below the smallest number a
double-precision float can represent (about 10^-308) and silently become zero —
**numerical underflow**, which destroys the computation without raising an error.

Two standard fixes, both used in this project's code:

1. **Normalize at each position.** After computing alpha_t, divide it by its own sum
   c_t and store that sum separately. The normalized vector is then a genuine
   probability distribution over states, and the total is recovered at the end as
   log P(observations) = sum over t of log c_t.
2. **Work in logarithms** throughout, replacing multiplication by addition and using a
   numerically careful "log-sum-exp" routine for the additions.

##### (k) What the forward pass does *not* give you

alpha_T normalized tells you the distribution of the state at the **final** position.
It is tempting to think you could get a whole path by doing that at every position
independently — computing each position's posterior marginal and drawing from each
separately. **That is wrong**, and it is worth seeing why.

Those per-position marginals here (computed by enumeration) are:

| | h1 | h2 | h3 |
|---|---|---|---|
| position 1 | 43.03% | 49.38% | 7.59% |
| position 2 | 0.92% | 8.01% | 91.07% |
| position 3 | 0.13% | 12.80% | 87.07% |

Drawing independently from these three rows could easily produce h2 → h1 → h3: two
switches, one of them into a haplotype that matches nothing, a path with negligible
joint probability. The marginals are each individually correct but they carry no
information about which *combinations* are coherent. A path must be drawn **jointly**.

*(An aside that shows the marginals are genuinely global: at position 1, h1 and h2
are not equally likely — 43.03% versus 49.38% — even though they are indistinguishable
at site 1 itself. The asymmetry comes from site **3**, where h2 matches and h1 does
not. Evidence from the far end of the sequence has propagated back to position 1. That
is a feature of the whole-sequence computation, not of position 1.)*

##### (l) Backward sampling, derived

The goal: draw a complete path in exact proportion to its posterior probability
P(z_1…z_T | observations). The trick is to factor that joint distribution backwards:

> P(z_1…z_T | X) = P(z_T | X) · P(z_{T-1} | z_T, X) · P(z_{T-2} | z_{T-1}, X) · …

This is just the chain rule of probability, applied right to left — always valid. Now
two simplifications make each factor computable.

**First:** the last factor is already in hand. P(z_T | X) is just alpha_T normalized.

**Second, the key step:** consider P(z_t | z_{t+1}, all observations). Once you
condition on the state at t+1, the observations *after* t+1 tell you nothing more
about z_t — any influence they have on z_t has to pass *through* z_{t+1}, and z_{t+1}
is already known. This is the Markov property again, running in reverse. So the future
observations drop out:

> P(z_t | z_{t+1}, X_1…X_T) = P(z_t | z_{t+1}, X_1…X_t)

And that we can evaluate with Bayes' rule:

> P(z_t = i | z_{t+1} = j, X_1…X_t)
> ∝ P(z_t = i, X_1…X_t) × P(z_{t+1} = j | z_t = i)
> = **alpha_t(i) × A(i, j)**

*In words: the weight of state i at position t is how plausible i was given everything
seen up to t (that is exactly what alpha stores), multiplied by how willing i is to
hand over to the state we have already committed to at t+1.* Normalize those K weights
and draw.

That is the whole algorithm — **forward filtering, backward sampling (FFBS)**: sweep
forward storing every alpha column, then sweep backward drawing one state at a time.
Note the storage requirement this implies: you must **keep all T alpha columns**, not
just the running one. Section 4 returns to this, because PanMixer's forward pass keeps
only the running column and therefore cannot sample this way at all.

##### (m) Backward sampling, worked

**Draw the last state.** Normalize alpha_3 = (0.000048, 0.004714, 0.032066) by its sum
0.036828:

> P(z_3 = h1) = 0.000048 / 0.036828 = **0.13%**
> P(z_3 = h2) = 0.004714 / 0.036828 = **12.80%**
> P(z_3 = h3) = 0.032066 / 0.036828 = **87.07%**

Suppose the draw gives **h3**.

**Draw position 2, given h3 at position 3.** Weights are alpha_2(i) × A(i, h3):

- h1: 0.003137 × 0.05 = 0.000157
- h2: 0.003137 × 0.05 = 0.000157
- h3: 0.035640 × 0.90 = 0.032076
- sum = 0.032390 → probabilities **0.48%**, **0.48%**, **99.03%**

Suppose the draw gives **h3** again.

**Draw position 1, given h3 at position 2.** Weights are alpha_1(i) × A(i, h3):

- h1: 0.330000 × 0.05 = 0.016500
- h2: 0.330000 × 0.05 = 0.016500
- h3: 0.003333 × 0.90 = 0.003000
- sum = 0.036000 → probabilities **45.83%**, **45.83%**, **8.33%**

Suppose the draw gives **h1**. The sampled path is **h1 → h3 → h3**.

**The consistency check.** Multiply the three conditional probabilities just used:

> 0.4583 × 0.9903 × 0.8707 = **0.3952**

Compare with the brute-force table in (h): the path h1 → h3 → h3 has posterior
probability **39.52%**. They agree exactly. This is not a coincidence — it is the
chain-rule factorization from (l) being reassembled, and it is the concrete
demonstration that FFBS samples from the true joint posterior rather than from some
convenient approximation.

*(Independently confirmed in this project by drawing 400,000 paths with FFBS and
comparing the empirical frequency of all 27 paths against the exact posterior: the
largest discrepancy was 0.07 percentage points, consistent with Monte-Carlo noise.)*

##### (n) Why a *sample*, and not the single best path

There is a well-known alternative, the **Viterbi algorithm**: the same dynamic program
with the sum replaced by a maximum, plus back-pointers, returning the single
highest-probability path. Here it would return h1 → h3 → h3 (or its tied twin), and
return it **every single time** — it is an **argmax**, deterministic by construction.

For this project that distinction is not a detail, it is the whole point. A
deterministic function of a private input leaks that input: run it twice, get the same
answer; the randomness that a privacy guarantee is purchased with simply is not there.
§1.3 makes this precise and Section 4 builds its guarantee on sampling. Section 3
measures a real case where a sampler *looked* random but was effectively deterministic,
which is the same failure wearing a disguise.

##### (o) Reading the result as biology

The sampled path h1 → h3 → h3 means: *copy haplotype h1 at site 1, then recombine, then
copy h3 for sites 2 and 3.* Emitting each state's own allele gives the synthetic
haplotype `0, 0, 1` — which is precisely the observed sequence, assembled as a
**mosaic of two panel haplotypes with one switch**, even though no single panel member
matched it.

That is the Li-Stephens idea in one sentence, and it is the engine underneath both
tools in this document. It also fixes what "broken" looks like: if the switch
probability were effectively zero, every sampled path would sit on one haplotype for
its entire length and the output would be a verbatim copy of one panel member. The
mosaic would be inert. **Section 3 measures exactly this failure in released code, on
real data.**

##### (p) The algorithm at a glance

| Step | Input | What happens | Output | Why it matters |
|---|---|---|---|---|
| 0. Set up | Panel (K haplotypes × T sites); the observed allele sequence | Fix pi (start), A (stay/switch), e (match/mismatch) | A fully specified HMM | A encodes recombination; e encodes mutation and sequencing error |
| 1. Forward filter | The HMM + observations | Sweep left to right: alpha_t(j) = [sum_i alpha_{t-1}(i)·A(i,j)] · e_t(j), normalizing each column to avoid underflow | **All T columns** of alpha, kept; the last column sums to P(observations) | Gives the likelihood, and stores everything sampling needs. O(T·K) under stay/switch |
| 2. Draw the final state | The last alpha column | Normalize and draw once | The path's state at position T | Weighted by the *entire* observation sequence, not just the last site |
| 3. Backward sample | alpha column t, plus the state already drawn at t+1 | Draw state t with probability proportional to alpha_t(i) · A(i, z_{t+1}) | The state at t; repeat down to t = 1 | Produces a *jointly* coherent path; independent per-site draws would not |
| 4. Read off the haplotype | The sampled state sequence | Emit each sampled state's panel allele at its site | A synthetic haplotype: a mosaic of panel haplotypes | This is the object released by the mechanism in Section 4 |

> **Key idea** — An HMM here answers "which panel haplotype is this new genome
> copying from, position by position?" You never see the answer directly, only the
> alleles it produced. The **forward algorithm** sums over the astronomically many
> possible answers efficiently, by noticing that paths through the same state can be
> pooled — giving both a likelihood and a stored table. **Backward sampling** then uses
> that table to draw one complete answer at random, correctly weighted, by starting at
> the end and walking back, each step asking "given where I go next, and everything I
> saw up to here, where was I?" The result is a mosaic: copy one haplotype for a
> stretch, switch, copy another. Sampling rather than taking the best path is what
> later makes a privacy guarantee possible.

#### 1.2.3 The Li–Stephens copying model

> **Definition** — **Three different things get called "a panel", and confusing them is the single most common way to misread this project.**
>
> - The **pangenome graph** is the data structure. It is not a panel at all; it is the object built *from* a cohort of assembled genomes, and a genome is a path through it (§1.1.6).
> - The **cohort** is that set of assembled genomes: here **44 individuals = 88 haplotypes**, measured in this project. The graph was built from them and they are public.
> - The **copying panel** is the set of haplotypes an HMM's states range over — the K in this section. In PanMixer the copying panel is **the cohort's own 88 haplotypes** when a block is scored, and **86** when a replacement is sampled, because the target's own two are removed first. It is *not* an external collection.
> - The **external reference panel** is 1000 Genomes: 3,202 samples on the build-correct data, measured in this project. PanMixer never copies from it. It supplies three other things — the LD block boundaries (Step 8), the population allele frequencies for one branch of the frequency table (Step 13), and the database the linkage attack searches (Step 21).
>
> So when this section says "the K haplotypes in the panel", read **copying panel**, and on this project's data read **the cohort's own haplotypes, minus the target**. §A.3 states the further restriction that matters for PanMixer: its chain is per block, over that block's anchor positions only.


The **Li–Stephens model** is the HMM above given a specific biological interpretation: **the states are the copying panel's haplotypes, and the hidden state at site t says "which of them is my chromosome copying from right here."** Because real chromosomes are mosaics of ancestral chromosomes (§1.1.4), any new haplotype ought to look like a *patchwork of existing ones with occasional switches*, plus a sprinkling of mismatches for new mutations and errors. That is exactly this model: **transitions = recombination, emissions = mutation.**

Two parameters carry all the biology:

- The **switch rate** must grow with the **genetic distance in centiMorgans** between consecutive sites, scaled by the **effective population size Ne** and by the panel size K. Far apart in cM ⇒ likely to switch; adjacent ⇒ almost certainly stay.
- The **mismatch rate** allows the copied allele to differ from the observed one.

**Effective population size Ne**, operationally: *Ne is the size of an idealised, randomly mating population that would exhibit the same amount of genetic drift — and the same rate of shared ancestry — as the real population.* It is far smaller than the census population because of past bottlenecks; humans are conventionally modelled with **Ne = 10,000**. Inside the Li–Stephens switch rate it acts as the constant that converts genetic distance into an expected number of ancestral recombination events. (PanMixer's published Fig. 6 uses Ne = 10,000 together with a standard recombination-rate constant r = 1.26; Section 3 works through that formula.)

Li–Stephens is the workhorse behind statistical phasing, imputation, and haplotype-based ancestry inference, and it is the model both tools in this project use. This also gives you the single most useful unifying observation available: **imputation is this same model run in the other direction.** Phasing asks "which panel haplotypes am I copying, given my genotypes"; imputation asks "given which panel haplotypes I appear to be copying, what must my allele be at a site I did not observe." So the reconstruction attack of §1.3.2 is not a separate technology — it is the machinery of this subsection pointed at a release.

Everything hinges on the switch rate being computed from **real genetic distance**. If the "distance" fed in is, say, *the difference of row indices in a VCF file*, the model is no longer a model of recombination at all, and the mosaic can collapse to verbatim copying. *Measured in this project:* executing the shipped `HaplotypeHMM` class on real chr21 data for subject **HG00438**, **300 of 300** randomly chosen HMM-path blocks returned a block **identical to a single donor haplotype** — no recombination at all. *(Recorded in `docs/memory.md` under 2026-09-18 and in `docs/bugs.md`.)* Section 3 dissects the four compounding causes.

#### 1.2.4 Optimisation: the knapsack problem and LP relaxation

The **0/1 knapsack problem**: you have n items; item i has value v_i and weight w_i; you have a capacity C; choose a subset maximising total value subject to total weight ≤ C. Each item is in or out — no fractions. This is NP-hard in general, though solvable by **dynamic programming** (building a table of best-value-for-each-partial-capacity, each cell computed once from earlier cells) in O(n·C) time when the weights are integers. That is called *pseudo*-polynomial: polynomial in the numeric *value* C, not in the number of digits needed to write C.

The **LP relaxation** — "LP" = **linear program**, an optimisation problem with a linear objective and linear constraints, solvable in polynomial time — replaces the binary decision variable x_i ∈ {0,1} with a continuous one, x_i ∈ [0,1]: you may take a *fraction* of an item. This makes the problem easy: sort items by **density** v_i / w_i ("value per unit of weight") and greedily take the best ones until capacity is exhausted, splitting at most one item. Two facts make the relaxation useful: the relaxed optimum is an **upper bound** on the true optimum (the relaxed problem permits everything the original does, and more), and the relaxed solution is nearly integral — **at most one item is fractional** — so rounding it yields a good feasible answer.

Why this appears here: one natural way to frame "how much privacy protection can I buy for a fixed amount of damage to the data" is a knapsack. Items = candidate edits, one per genomic block; **value** = privacy risk removed; **weight** = utility destroyed; **capacity** = how much utility loss you are willing to tolerate. That is precisely the framing used by PanMixer (the published tool dissected in Section 3).

*Measured in this project, running the shipped optimiser on our rebuilt chr21 for subject **HG00438** (seeds 123 and 456, results byte-identical for a fixed seed):*

| Capacity setting | Moves made | Utility loss spent | As a percentage of the budget denominator |
|---|---|---|---|
| 0.1 | **47,070** | 928.92 | **10.0%** |
| 0.5 | **177,865** | 3,077.60 | **33.1%** |

*(Both rows are measured in this project; the source is `logs/pm_obf.log`, logged in `docs/memory.md` under 2026-09-26.)*

Two things to read off this. At capacity 0.1 the budget binds exactly — the mechanism spends what it is allowed. At capacity 0.5 the knapsack **saturates below its stated budget**: it runs out of moves worth taking before it runs out of capacity, so capacity stops being a meaningful knob at the top end. Second, and separately, the *denominator* against which those percentages are computed is itself a constructed quantity — twice the sum of all per-move utility losses (9,288.85 for this subject-chromosome) — and the total achievable loss is strictly less than it. So "33.1% of the budget" is a fraction of a number the mechanism can never reach. Both facts have to be stated together, or the saturation looks mysterious.

#### In plain words (computer science)

Strip away the biology and a cohort of genomes is a big matrix of small integers: one row per variant, one column per haplotype, entries naming which allele that haplotype carries. Because the matrix is huge and we constantly need to line up rows between different files, everything is sorted by genomic coordinate, indexed for random access, and joined through dictionaries keyed on (position, reference allele, alternate allele) — a join that silently returns almost nothing if the two files use different coordinate systems, which we measured here as 1,025 matches against a GRCh37 panel versus 258,610 against a GRCh38 one, on identical pangenome records. The key statistical tool is the hidden Markov model: a chain of unobserved states that emit the things you do observe, with a forward pass that accumulates probabilities left to right and a backward pass that samples a coherent whole path. In genomics the states are "which existing haplotype am I copying from right now," which is the Li–Stephens model: a new chromosome is a patchwork of old ones with occasional switches caused by recombination and occasional mismatches caused by mutation — and the very same model, run in the other direction, is what imputation and reconstruction attacks use. The switch probability has to be driven by genetic distance in centiMorgans; if it is driven by something else, the model quietly degenerates into copying a single haplotype whole, which is what we measured on 300 of 300 blocks. When the transition matrix is just "stay or switch," the forward pass costs O(T·K) rather than O(T·K²) — one shared scalar per step plus a per-state correction — which is what makes these methods cheap. Finally, when you must pick a subset of edits under a budget, that is a knapsack problem: NP-hard in general, easy to relax into a linear program that gives an upper bound and a near-integral solution.

---

### 1.3 The security and the mathematics

#### 1.3.1 Why genomic data is unusually sensitive

Most personal data is sensitive because of what it says about you today. Genomic data is worse along four axes:

1. **It is uniquely identifying.** Consider independent, roughly balanced biallelic sites: each one splits the population approximately in half, so each contributes up to one **bit** of information, and n bits distinguish up to 2ⁿ distinct patterns. Since 2³³ ≈ 8.6 billion exceeds the world population, on the order of a few dozen independent common variants are in principle enough to single out one human being. Real variants are correlated rather than independent (§1.1.4), so the practical number is larger — but it is still *dozens*, not thousands. **A "small" genomic release is not small.**
2. **It is immutable.** You can change a password, a credit card number, or an address. You cannot change your genome. A breach today is a breach forever.
3. **It implicates people who never consented.** You share about half your variants with a parent, a child, or a full sibling, and about a quarter with a grandparent or half-sibling — including relatives not yet born. Publishing your genome partially publishes theirs.
4. **It predicts.** Genotype carries information about ancestry, about physical traits, and about disease risk — including conditions that have historically attracted discrimination in employment and insurance.

The field also carries specific historical weight. The published PanMixer paper (`paper/s41467-026-77591-0_reference.pdf`) opens its introduction by grounding the privacy concern in historical cases: Henrietta Lacks, whose cells (HeLa) were taken and distributed without consent and contributed to major biomedical advances while creating lasting mistrust; and the unauthorised use of Havasupai Tribe samples for research beyond the purpose consented to. Its stated conclusion is that for the pangenome to represent all populations, especially historically under-represented ones, privacy has to be built in from the start. The practical point for this project is the same: privacy protection is not only about preventing harm to today's donors, it is about making it *possible* for people — especially from populations with the most reason for mistrust — to contribute at all.

#### 1.3.2 The four attack shapes

- **Re-identification.** Given a record stripped of names, determine whose it is. In genomics this usually means matching the record against some external source keyed on genetic content.
- **Linkage attack.** The specific mechanism behind most re-identification: take the released dataset, take an **auxiliary** database the attacker already holds (another genomic panel, a direct-to-consumer ancestry database, a forensic database, a published study), and **join them on shared attributes**. Because genomes are so identifying, a join on a few hundred variants is effectively a join on a primary key. In this project that auxiliary database is called the *attack database*, and the attack works by scoring how well each released haplotype matches each database individual and checking whether the true individual comes out on top.
- **Membership inference.** A weaker but often sufficient question: *was person X in this dataset at all?* If the dataset is "people with schizophrenia," membership alone is the sensitive fact. **Standard field background:** the landmark result here is usually credited to **Homer and colleagues (2008)**, who showed that even *aggregate allele frequencies* — no individual records released at all — can reveal whether a known individual was part of the pool. This upended the prior assumption that summary statistics are inherently safe. It is a real and widely reproduced result, but it is background you should verify in the original literature rather than something measured in this project. (Note: the published PanMixer paper adds a membership-inference analysis that its earlier preprint did not contain.)
- **Reconstruction / imputation attack.** Fill in what was withheld. Because of LD (§1.1.4), variants that were masked, removed, or perturbed can often be predicted from their neighbours plus a public reference panel — using, as §1.2.3 noted, the very same Li–Stephens machinery. This is why "we deleted the sensitive variants" is not a defence, and why any perturbation scheme must be evaluated against an attacker who runs imputation on the release.

#### 1.3.3 Randomized mechanisms and indistinguishability

A **randomized mechanism** M is an algorithm that takes a private input and produces a *random* output. It is fully described by the conditional distribution of its output given its input: for input p, write **Q_p** for "the probability distribution over outputs when the input is p." The same input run twice gives different outputs; that is the point, not a defect.

The entire privacy question then becomes: **if two different private inputs p and q would produce output distributions that look nearly the same, then seeing the output tells you almost nothing about which input it came from.** So we need a way to say "nearly the same" for probability distributions — and we need it to mean something operational, not merely aesthetic.

#### 1.3.4 A hypothesis-testing interlude (four words you will need)

The next two subsections speak the language of statistical testing, which is not part of the assumed background, so here are the four words, each defined once:

- A **simple hypothesis** fully specifies a probability distribution, with no unknown parameters left over. "The data came from distribution P" is simple; "the data came from *some* normal distribution" is not.
- The **false-positive rate** (also called the significance level, or Type I error rate) is the probability that a test raises a flag when the flagged hypothesis is actually false.
- The **power** of a test is the probability that it correctly raises a flag when the hypothesis *is* true. Good tests have high power at a fixed false-positive rate.
- A **null distribution** is what a statistic would look like if the thing you are looking for were *not* there — e.g. what a matching score would look like if the target individual were not in the dataset. A **z-score** is (observed value − mean) / standard deviation: how many standard deviations from typical an observation is, relative to that null.

Throughout this document, "an attacker" is a hypothesis test: it observes a release and must decide between "this came from person A" and "this came from person B," or between "person A is in the dataset" and "person A is not."

#### 1.3.5 Total variation distance, and why the optimal attack succeeds with probability (1 + TV)/2

For two probability distributions P and Q over the same finite set Y, the **total variation distance** is

> **TV(P, Q) = maximum over subsets A of Y of | P(A) − Q(A) | = (1/2) × sum over y in Y of | P(y) − Q(y) |.**
>
> *In words: TV is the largest possible disagreement between the two distributions about the probability of any single event.* TV = 0 means the distributions are identical; TV = 1 means they have disjoint supports and can be told apart perfectly. (Our outcome sets are finite, so "maximum over events" is exact; in a more general setting one would say *supremum*, the least upper bound, which is the same thing here.)

Now the interpretation that makes TV the right currency. Suppose an adversary knows both P and Q, sees one sample y, and must guess which distribution produced it. Suppose the two hypotheses are equally likely a priori (**equal priors**, probability 1/2 each). The best strategy is obvious: guess whichever distribution assigned the *higher* probability to the observed y. (Ties may be broken arbitrarily — when P(y) = Q(y), either guess is right with the same probability, so the answer is unaffected.) The success probability is then

> P(success) = sum over y of (1/2) · max{ P(y), Q(y) }.

Use the identity max{a, b} = (a + b)/2 + |a − b|/2:

> P(success) = (1/2) · [ (1/2)·sum_y (P(y) + Q(y)) + (1/2)·sum_y |P(y) − Q(y)| ]
> = (1/2) · [ (1/2)(1 + 1) + (1/2)(2·TV) ]
> = **(1 + TV(P, Q)) / 2.**
>
> *In words: a blind coin flip gets you 1/2; the total variation distance is exactly the extra edge the data hands the very best possible attacker, and half of it shows up in the success rate.*

That is the whole argument, and it is worth internalising because it converts an abstract distance into an operational promise: **if you can prove TV ≤ tau, you have proved that no attacker — however clever, however well-resourced, using any algorithm at all — can do better than (1 + tau)/2 at telling the two inputs apart from one release.** With tau = 0.10 the ceiling is 0.55: barely better than the 0.50 you get by flipping a coin without looking at the data at all.

Two cautions, both of which become conditions later. First, this is the bound for **one** release; seeing several independent samples lets an attacker do better, which is why "one release per genome" will appear as an explicit requirement. Second, "equal priors" is an assumption; an attacker who starts with a strong prior does better than (1 + tau)/2 overall — but that improvement comes from the prior, not from the release, and the release's marginal contribution is still bounded.

#### 1.3.6 Likelihood ratios, differential privacy, and the privacy budget

The **likelihood ratio** of an observed output y under two hypotheses is LR(y) = P(y) / Q(y) — *how many times more probable this particular observation is under one hypothesis than the other.* The **Neyman–Pearson lemma** says that for two simple hypotheses, the optimal test always thresholds this ratio: no other test achieves higher power at the same false-positive rate. So controlling the likelihood ratio **pointwise** — for every possible output — is the strongest natural form of protection.

A pointwise bound, "LR(y) ≤ e^epsilon for **all** y," is strictly stronger than a bound on TV, which is an *average* over outputs: a mechanism can have small TV while still possessing some rare output that is a thousand times more likely under one input than the other. On that rare day the attacker is certain. Bounding the ratio pointwise forbids that.

**Differential privacy (DP)** is the standard formalisation of exactly this idea. A mechanism M is **epsilon-differentially private** if for every pair of "neighbouring" datasets D and D′ (differing in one individual's record) and every set S of possible outputs,

> **P(M(D) ∈ S) ≤ e^epsilon × P(M(D′) ∈ S).**
>
> *In words: whether or not any one person's record is in the data can change the probability of any outcome by at most a factor of e^epsilon.*

The parameter **epsilon** is the **privacy budget** (or privacy loss): the logarithm of the worst-case likelihood ratio. epsilon = 0 means the output is literally independent of any individual's data (perfectly private, and usually useless); small epsilon (say ≤ 1) is strong; large epsilon (≥ 10) is mostly ceremonial. Budgets **compose**: answering an epsilon₁-private query and an epsilon₂-private query about the same data costs epsilon₁ + epsilon₂ in total, which is why "how many releases?" is a first-class question and not a detail. *(Recall from §1.0 that this epsilon is not PanMixer's ε_j, which is a per-block privacy risk in nats, nor ε_private = 0.002, which is an operating threshold that paper selects.)*

**Local differential privacy (LDP)** removes the trusted curator: each individual randomizes their *own* record before it ever leaves their hands, so no one — not even the data collector — ever sees the truth. The textbook example is **randomized response**: asked a yes/no question, you flip a coin; heads, answer truthfully; tails, flip again and answer "yes" on heads, "no" on tails. Then P(answer = yes | truth = yes) = 3/4 and P(answer = yes | truth = no) = 1/4, so the likelihood ratio is exactly 3 and the mechanism is ln(3) ≈ 1.10-locally-private. Crucially the aggregate is still recoverable: if the observed yes-rate is r, the estimated true rate is 2r − 1/2, because the noise distribution is known and unbiased. The lesson generalises — *calibrated noise with a known distribution destroys individual-level certainty while preserving population-level signal.* The mechanism proposed in this project is closest in spirit to the local model: it protects one individual's own data at the moment of release, rather than protecting contributors to an aggregate held by a trusted curator.

The two currencies relate cleanly. If a mechanism's likelihood ratio is bounded by e^epsilon everywhere, then its TV distance satisfies

> **TV ≤ (e^epsilon − 1)/(e^epsilon + 1) = tanh(epsilon/2).**

Run that backwards: to guarantee TV ≤ tau it suffices to take **epsilon = 2·arctanh(tau)**, since tanh(arctanh(tau)) = tau exactly. Keep that identity in your pocket — it is precisely why the function **arctanh** appears in this project's mechanism, and it means the pairwise TV bound and the pointwise likelihood-ratio bound are two faces of one calibration rather than two independent claims.

#### 1.3.7 The distinction that matters most: a proved bound versus a failed attack

This is the conceptual point separating privacy *research* from privacy *theatre*, and you should be able to state it in your sleep.

| | Worst-case distributional guarantee | Empirical attack evaluation |
|---|---|---|
| **Form of the claim** | "For *all* attackers, *all* side information, *all* pairs of inputs, success ≤ (1+tau)/2" | "*This* attack, with *this* auxiliary database, on *this* data, achieved x% success" |
| **Quantifier** | Universal (∀) | Existential (∃, and only over what was actually tried) |
| **Can it be falsified later?** | No — it is a theorem, given its assumptions | Yes — a better attack or a better auxiliary database can appear tomorrow |
| **Failure mode** | Assumptions violated in implementation (e.g. more than one release; an unbounded utility function; an output support that secretly depends on the private input) | Overfitting to the tested attack; a false sense of safety |
| **What it cannot tell you** | How much utility survives, or how bad realistic attacks actually are (bounds are often very loose) | Anything at all about attackers you did not run |

Both are necessary. A theorem with no empirical evaluation may be vacuously safe and useless; an empirical evaluation with no theorem is only a promise that the attacks the authors thought of did not work. The honest position is to state both, and never to let one masquerade as the other.

A direct consequence for this project, worth stating now so that it does not surprise you in Section 4: **the two approaches compared in this document measure privacy in formally incomparable units.** PanMixer's privacy quantity is a per-block score computed from **pointwise mutual information (PMI)** — PMI(x, y) = log[ P(x, y) / (P(x)·P(y)) ], *the log of how often two things co-occur divided by how often they would co-occur by chance*; for an obfuscated block it reduces to the **self-information** −log p(h_bj) of §1.1.3, an average-case information-theoretic quantity in nats. Our tau is a **worst-case bound on the distance between two output distributions**. Neither can be converted into the other. The only place the two can be compared head-to-head is on **measured empirical attack success** — run the same attack against both releases and count.

#### In plain words (security and mathematics)

Genomes are unusually dangerous data: a few dozen independent common variants carry enough bits to pick one person out of the world's population, you cannot change your genome after a leak, your genome partly discloses your relatives' genomes including relatives not yet born, and it predicts ancestry, traits and disease risk. The attacks worth naming are re-identification (whose record is this?), linkage (join the release against some other database the attacker already holds), membership inference (was this person in the dataset at all? — the classic result, credited to Homer and colleagues in 2008, is that even published allele-frequency summaries can answer this), and reconstruction (predict what was withheld using the correlations between nearby variants, which is just imputation pointed the wrong way). The defensive idea is to make the release *random*, so that two different secret inputs produce output distributions that look almost the same. "Almost the same" is measured by total variation distance — the biggest disagreement the two distributions can have about any event — and it has an exact operational meaning: with equal priors, the best possible attacker guesses right with probability (1 + TV)/2, so bounding TV by 0.10 caps *any* attacker at 55% against the 50% they get by guessing blind. A stronger relative is bounding the likelihood ratio at every single possible output, which is what differential privacy does; its parameter epsilon is the log of that worst-case ratio, smaller is more private, and budgets add up across releases, so the number of releases is part of the guarantee. Local differential privacy is the version where you randomize your own data before handing it over, as in the randomized-response coin-flip trick, which destroys individual certainty while leaving population averages recoverable. The two currencies line up exactly at epsilon = 2·arctanh(tau). Finally, never confuse a theorem ("no attacker can beat this bound") with an experiment ("the attacks we tried failed"): the first is universal but can be loose and depends on its assumptions surviving the implementation, the second is concrete but can be overturned by tomorrow's better attack. You need both, and you must always label which one you are claiming.

---

### 1.4 What this project is trying to achieve

There are **two distinct privacy settings** hiding inside the phrase "pangenome privacy," and the whole project turns on keeping them apart.

**Setting A — protecting a contributor to the graph.** Here the sensitive individual is *inside* the public resource: they donated their genome, it was assembled, and their haplotypes became paths in the released pangenome graph. The graph is the product, and their data is part of the product. The privacy question is: can we perturb the resource so that the contributor cannot be re-identified, and their sensitive variants cannot be read off, while the resource stays scientifically useful for the things people build pangenomes to do — estimating allele frequencies, studying linkage disequilibrium, mapping reads with less reference bias? This is the setting PanMixer addresses, on the HPRC draft pangenome (44 individuals, 88 haplotypes in the configuration this project rebuilt). Its release object is the cohort's variant file with the protected contributor's sample column rewritten — note that the graph *topology* is not modified — and its notion of success is that measured linkage attacks and genome-reconstruction attempts fail while the downstream statistics stay close to the originals.

**Setting B — protecting an external individual whose path is released.** Here the sensitive individual is *outside* the public resource entirely. The pangenome graph G and the public cohort D that built it are fixed and public, and they stay fixed and public — nothing is ever written back into them. A new person, not in D, has their genome mapped onto G, and the object we want to release is **their path through the graph**. That path can reveal rare alleles, unusual combinations of alleles, and other identifying features, even though the graph itself is entirely public. The goal is a **randomized path release mechanism** with a *proved* guarantee — specifically, that for any two admissible input paths the resulting output distributions differ by at most a pre-specified total variation distance tau, so that by the (1 + TV)/2 argument of §1.3.5 no equal-prior binary attacker can identify the source of a release with probability above (1 + tau)/2. This is the setting of this project's own proposal, *Private Genome Path Release against a Public Pangenome Graph* (`paper/Private_Genome_Path_Release.pdf`).

Its machinery is exactly the machinery assembled above. A public, target-independent cohort haplotype HMM supplies a baseline distribution R_D over cohort-supported paths. A bounded utility function u(p, y) ∈ [0, 1] nudges that distribution toward paths resembling the target. The size of the nudge is calibrated by eta_tau = arctanh(tau), so that the nudge can never move the distribution further than tau. And the release is drawn by exact forward-filtering/backward-sampling at cost O(T·K) — on chr21, roughly 100,757 blocks × 88 haplotype states ≈ 8.9 million operations, which is nothing.

**How much of PanMixer can be reused for Setting B?** A rough estimate, not a line-by-line count, and stated as such: its **data**, its **block structure**, the **algebra of its forward recursion**, and its **evaluation suite** are reusable — that is the majority of its code by volume. The parts that actually constitute the new mechanism — an across-block sampler with a backward pass, the privacy score, and the release rule — are not present in it and are new work. Treat "roughly 55–65% reusable by volume, well under 10% of the mechanism" as an order-of-magnitude estimate with no counting method behind it; the qualitative statement is the reliable one. *(The reuse assessment and the corrections behind it are recorded in `docs/memory.md` under 2026-09-18.)*

For orientation, the symbols that recur in Sections 2 through 4, stated in words. Re-read §1.0 before carrying any of these across into PanMixer's notation.

| Symbol | Read it as |
|---|---|
| **G** | the fixed, public pangenome graph |
| **D** | the public cohort of genomes the graph was built from |
| **N** | the number of cohort individuals (44 here), so 2N = 88 haplotypes |
| **p** | the private input: the target individual's mapped path through G |
| **q** | any other admissible input path (used only in the pairwise guarantee) |
| **y** | a candidate output path |
| **R_D** | the public, target-independent baseline distribution over cohort-supported paths |
| **Q_p** | the output distribution of the mechanism when the private input is p |
| **tau (τ)** | the privacy parameter, between 0 and 1: the maximum allowed total variation distance between output distributions for any two inputs |
| **eta_tau (η_τ)** | arctanh(tau): the **tilt strength** — how hard the mechanism is permitted to lean toward the target. *Not* PanMixer's η_j, which is a utility loss |
| **epsilon_tau (ε_τ)** | 2·arctanh(tau): the corresponding pointwise likelihood-ratio bound, in log units. *Not* PanMixer's ε_j, which is a per-block privacy risk |
| **u(p, y)** | the utility of releasing y when the truth is p, bounded between 0 and 1 |
| **T** | the number of blocks (or sites) along the chromosome |
| **K** | the number of cohort haplotype states (88 here) |
| **beta_t (β_t)** | the non-negative weight of block t in the utility, summing to 1 across blocks |
| **phi_t (φ_t)** | the per-block agreement score between input and output paths, between 0 and 1 |

Sections 2 through 4 take all of this apart: **Section 2** examines the problem itself in detail, with no solution attached; **Section 3** dissects PanMixer's theory and its shipped implementation step by step, including ten places where the released code departs from the published Methods, measured here; and **Section 4** does the same for this project's mechanism and compares the two.

---

### In plain words — Section 1

A human carries two copies of each autosome — 24 distinct chromosome sequences exist in the species, 23 pairs in a person, about 3.1 billion base pairs per haploid copy and roughly twice that per diploid cell — and at millions of positions those copies differ from each other and from other people's. A variant is such a position, an allele is one of the alternatives there, a genotype is the unordered pair you carry, and a haplotype is the ordered run of alleles along one physical chromosome copy. Haplotypes are much more identifying than genotypes, and rare alleles much more identifying than common ones, because the information in an observation is minus the logarithm of its frequency — bits with log base 2, nats with natural log. Because chromosomes are reshuffled only once or twice per generation, nearby variants travel together (linkage disequilibrium, measured as r²), which both lets us cut the genome into blocks and means that deleting a variant does not hide it, since its neighbours predict it; genetic distance for this purpose is counted in centiMorgans, not base pairs. Genomics has traditionally aligned everyone to one linear reference sequence, which systematically under-detects variation in people unlike that reference's donors, because a read carrying sequence the reference lacks has nowhere to align. A pangenome graph replaces that with many genomes woven into one graph where every branch point is a variant and every individual is a path — and a graph can be flattened into a VCF table, one row per branch point, one column per haplotype, holding small integers. That table is the data model: measured here on chromosome 21 it is 340,824 rows by 88 haplotypes (44 individuals, after the pipeline drops the CHM13 reference cell line — which is human — and chrX), carved into 100,757 block-dictionary entries of which 88.4% **of entries** are singletons holding one variant each, accounting for 26.1% **of variants**, while the 10.4% of entries carrying two or more anchors hold 71.4% of variants. Coordinates matter absolutely: the shipped pipeline joins a GRCh38 pangenome against a GRCh37 panel (1000 Genomes Phase 3, 2,504 samples), which we measured as 1,025 matching records where the build-correct 30x GRCh38 panel (3,202 samples, 1,002,753 chr21 records) gives 258,610. The statistical engine throughout is the hidden Markov model in its Li–Stephens form, where the hidden state is "which existing haplotype am I copying from right now," the transitions encode recombination measured in centiMorgans, the emissions encode mutation, and a forward pass plus a backward sampling pass draws a coherent mosaic haplotype in O(T·K) time — unless the switch rate is computed from the wrong quantity, in which case the mosaic collapses into verbatim copying, which is what we measured on 300 of 300 chr21 blocks for subject HG00438. On the privacy side, genomes are immutable, implicate relatives who never consented, and are so identifying that a few dozen independent variants pin down a person; the attacks to know are re-identification, linkage against an external database, membership inference, and reconstruction by imputation. The defence is randomization: make the release a random draw whose distribution barely depends on the secret. Total variation distance measures "barely," and it has an exact operational meaning — with equal priors the best possible attacker succeeds with probability (1 + TV)/2, since the optimal test picks the larger likelihood and max(a,b) = (a+b)/2 + |a−b|/2. Bounding the likelihood ratio at every output, which is what differential privacy's epsilon does, is the stronger pointwise version, and the two calibrations coincide exactly at epsilon = 2·arctanh(tau). Keep two kinds of claim rigidly apart: a theorem that bounds *every* attacker under stated assumptions, and an experiment showing that *the attacks we ran* failed. And keep the notation of the two works apart too: PanMixer writes η for utility loss and ε for privacy risk, while this project writes η_τ for tilt strength and ε_τ for a budget. This project lives in two settings — protecting someone whose genome is inside the public graph (PanMixer's problem) and protecting an external person whose mapped path through a fixed public graph is what gets released (our problem) — and the rest of the document takes each apart in turn. Two refinements sit on top of that picture, and both are load-bearing later. First, the graph-to-table correspondence is exact only for a SNP: an insertion or a deletion has to be written with a shared padding base, so a record describes a **span** rather than a point; bubbles **nest**, so on our chr21 **31,489 records** sit inside another record's allele (measured in this project, **9.2%** of records, derived) and the same DNA is described twice at two resolutions; a haplotype can legitimately have **no value at all** at a variant, for at least five distinguishable reasons that the conversion into the matrix collapses into a single `-1`; the graph can contain **cycles**, which is how a tandem duplication is stored without writing the second copy out at all; and at the extreme one chr21 site carries **88 distinct alleles across 88 haplotypes** in a single **28 MB** VCF line with **2,464** nested children (all measured in this project) — the most identifying record on the chromosome, and the case every later section is tested against. §1.1.8 works all of that through. Second, and following from it, an LD block is a statement about **correlation**, not about **geometry**: **15,540 of 100,757** block-dictionary entries (**15.4%**, measured in this project) contain a record whose sequence span crosses a block boundary, so blocks are **not** disjoint in sequence, and every later piece of accounting that adds up block by block is leaning on an assumption the data does not quite grant.