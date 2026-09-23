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

*Measured in this project:* on our rebuilt chr21 there are **340,824 variant records**, of which **274,458 are biallelic SNPs (80.5%)**; the remainder are indels, multi-allelic sites, and structural variation. The 274,458 figure is the count of set entries in the biallelic-SNP mask that the shipped preprocessing produces (`build_biallelic_snp_mask`); the chr21 artifact set is registered in `docs/data.md` and the preprocessing run is logged in `logs/pm_chr21.log`.

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

#### In plain words (biology)

You have two copies of every autosome, one from each parent; a haploid human genome is about 3.1 billion base pairs, so a diploid cell holds about twice that, spread across 22 numbered autosomes plus X or Y — 24 distinct chromosome sequences in the species, 23 pairs in a person. At millions of positions people carry different DNA letters. Those positions are variants, the alternatives are alleles, and we label them 0 for "same as the reference" and 1, 2, … for the alternatives; a variant with just one alternative is biallelic, with more it is multi-allelic. A genotype tells you *which two* alleles you carry at a position but not which physical chromosome copy each sits on; a haplotype tells you the whole ordered run of alleles along one copy, and working out which is which is called phasing. Haplotypes are far more identifying than genotypes. Most variants are rare, and rare variants are exactly the ones that pick a person out of a crowd, because the information in seeing an allele of frequency f is about −log f — bits if you use log base 2, nats if you use natural log. Because chromosomes are reshuffled only once or twice per generation, nearby variants travel together; that correlation is linkage disequilibrium, it lets us chop the genome into blocks, and it means deleting one variant does not really hide it because its neighbours predict it. Genetic distance for this purpose is measured in centiMorgans, not base pairs. Sequencing machines emit short reads, which are either aligned onto a reference or assembled into long contiguous sequence; alignment can only find what the reference has a slot for, which is why a single linear reference systematically under-detects variation in people unlike its donors. A pangenome graph fixes that by storing many genomes at once as paths through one shared graph, where each branch point (bubble) is a variant. A VCF file is a big table with one row per variant and one column per sample holding small integers, and because each bubble maps to one row, the graph and the table are two views of the same data. On our rebuilt chr21 that table is 340,824 rows by 88 haplotypes, and it is carved into 100,757 block-dictionary entries — of which 88.4% *of entries* are singletons, which is only about 26% *of variants*.

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
| **States** | 1…K | The K haplotypes in the reference panel — "which panel haplotype am I currently copying from?" |
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

> **The short version.** An HMM here answers "which panel haplotype is this new genome
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

The **Li–Stephens model** is the HMM above given a specific biological interpretation: **the states are panel haplotypes, and the hidden state at site t says "which panel haplotype is my chromosome copying from right here."** Because real chromosomes are mosaics of ancestral chromosomes (§1.1.4), any new haplotype ought to look like a *patchwork of existing ones with occasional switches*, plus a sprinkling of mismatches for new mutations and errors. That is exactly this model: **transitions = recombination, emissions = mutation.**

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

*(Raw output in `logs/pm_obf.log`; summarised in `docs/memory.md` under 2026-09-18.)*

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

Sections 2 through 4 take all of this apart: **Section 2** examines the problem itself in detail, with no solution attached; **Section 3** dissects PanMixer's theory and its shipped implementation step by step, including four defects measured here; and **Section 4** does the same for this project's mechanism and compares the two.

---

### In plain words

A human carries two copies of each autosome — 24 distinct chromosome sequences exist in the species, 23 pairs in a person, about 3.1 billion base pairs per haploid copy and roughly twice that per diploid cell — and at millions of positions those copies differ from each other and from other people's. A variant is such a position, an allele is one of the alternatives there, a genotype is the unordered pair you carry, and a haplotype is the ordered run of alleles along one physical chromosome copy. Haplotypes are much more identifying than genotypes, and rare alleles much more identifying than common ones, because the information in an observation is minus the logarithm of its frequency — bits with log base 2, nats with natural log. Because chromosomes are reshuffled only once or twice per generation, nearby variants travel together (linkage disequilibrium, measured as r²), which both lets us cut the genome into blocks and means that deleting a variant does not hide it, since its neighbours predict it; genetic distance for this purpose is counted in centiMorgans, not base pairs. Genomics has traditionally aligned everyone to one linear reference sequence, which systematically under-detects variation in people unlike that reference's donors, because a read carrying sequence the reference lacks has nowhere to align. A pangenome graph replaces that with many genomes woven into one graph where every branch point is a variant and every individual is a path — and a graph can be flattened into a VCF table, one row per branch point, one column per haplotype, holding small integers. That table is the data model: measured here on chromosome 21 it is 340,824 rows by 88 haplotypes (44 individuals, after the pipeline drops the CHM13 reference cell line — which is human — and chrX), carved into 100,757 block-dictionary entries of which 88.4% **of entries** are singletons holding one variant each, accounting for 26.1% **of variants**, while the 10.4% of entries carrying two or more anchors hold 71.4% of variants. Coordinates matter absolutely: the shipped pipeline joins a GRCh38 pangenome against a GRCh37 panel (1000 Genomes Phase 3, 2,504 samples), which we measured as 1,025 matching records where the build-correct 30x GRCh38 panel (3,202 samples, 1,002,753 chr21 records) gives 258,610. The statistical engine throughout is the hidden Markov model in its Li–Stephens form, where the hidden state is "which existing haplotype am I copying from right now," the transitions encode recombination measured in centiMorgans, the emissions encode mutation, and a forward pass plus a backward sampling pass draws a coherent mosaic haplotype in O(T·K) time — unless the switch rate is computed from the wrong quantity, in which case the mosaic collapses into verbatim copying, which is what we measured on 300 of 300 chr21 blocks for subject HG00438. On the privacy side, genomes are immutable, implicate relatives who never consented, and are so identifying that a few dozen independent variants pin down a person; the attacks to know are re-identification, linkage against an external database, membership inference, and reconstruction by imputation. The defence is randomization: make the release a random draw whose distribution barely depends on the secret. Total variation distance measures "barely," and it has an exact operational meaning — with equal priors the best possible attacker succeeds with probability (1 + TV)/2, since the optimal test picks the larger likelihood and max(a,b) = (a+b)/2 + |a−b|/2. Bounding the likelihood ratio at every output, which is what differential privacy's epsilon does, is the stronger pointwise version, and the two calibrations coincide exactly at epsilon = 2·arctanh(tau). Keep two kinds of claim rigidly apart: a theorem that bounds *every* attacker under stated assumptions, and an experiment showing that *the attacks we ran* failed. And keep the notation of the two works apart too: PanMixer writes η for utility loss and ε for privacy risk, while this project writes η_τ for tilt strength and ε_τ for a budget. This project lives in two settings — protecting someone whose genome is inside the public graph (PanMixer's problem) and protecting an external person whose mapped path through a fixed public graph is what gets released (our problem) — and the rest of the document takes each apart in turn.