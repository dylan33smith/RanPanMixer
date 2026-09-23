#!/bin/bash
# Reproduce the PUBLISHED read-mapping BASELINE (unobfuscated chr21 graph):
# paper reports 77.83% perfect / 95.62% gapless / 77.01% MAPQ60.
# Follows tools/downstream/utility_out/vg_prep.py + quick_align.py, with one
# deviation: the contig rename is done FIRST, because the shipped order runs
# `bcftools norm -f` while the VCF still uses grch38#chr21 naming, which no
# reference FASTA in the repo provides.
set -euo pipefail
E=/home/ds85/miniconda3/envs/panmixer/bin
export PATH=$E:$PATH
P=/data/ds85/RanPanMixer/external/PanMixer
VG=$P/downloaded_tools/vg
REFDIR=$P/starting_data/references
OUT=/data/ds85/RanPanMixer/runs/vg_baseline_ucsc
mkdir -p $OUT
step(){ echo; echo "######## $* :: $(date +%H:%M:%S)"; }

step "0. rename map + inputs"
echo "grch38#chr21 chr21" > $OUT/rename_chrs.txt
ls -la $REFDIR/chr21_ucsc.fa $P/starting_data/chr21/pangenome.vcf.gz

step "1-5. VCF prep (view -c1 | norm | trim-alt | SNPs only)"
bcftools annotate --rename-chrs $OUT/rename_chrs.txt $P/starting_data/chr21/pangenome.vcf.gz -Ou \
 | bcftools view -c 1 -Ou \
 | bcftools norm -N -f $REFDIR/chr21_ucsc.fa -Ou \
 | bcftools view --trim-alt-alleles -Ou \
 | bcftools view -t chr21 -v snps -Ou \
 | bcftools sort -T $OUT/sorttmp -Oz -o $OUT/baseline.final.vcf.gz
tabix -p vcf -f $OUT/baseline.final.vcf.gz
echo "records in graph VCF: $(bcftools index -n $OUT/baseline.final.vcf.gz)"

step "6. vg autoindex (giraffe)"
rm -f $OUT/vg_idx*
$VG autoindex --workflow giraffe --prefix $OUT/vg_idx --ref $REFDIR/chr21_ucsc.fa --vcf $OUT/baseline.final.vcf.gz -t 32
ls -la $OUT/vg_idx*

step "7. giraffe align + vg stats, 5 external donors (paper's list)"
GBZ=$OUT/vg_idx.giraffe.gbz
DIST=$(ls $OUT/vg_idx*.dist | head -1)
MIN=$(ls $OUT/vg_idx*.min | head -1)
ZIP=$(ls $OUT/vg_idx*.zipcodes 2>/dev/null | head -1 || true)
echo "index: GBZ=$GBZ DIST=$DIST MIN=$MIN ZIP=${ZIP:-none}"
for S in HG00138 HG00635 HG01112 HG02698 NA18853; do
  echo "--- $S $(date +%H:%M:%S)"
  FQ=$P/read_fastqs/${S}.fastq.gz
  [ -s "$FQ" ] || cat $P/read_fastqs/${S}_1.fastq.gz $P/read_fastqs/${S}_2.fastq.gz > "$FQ"
  ARGS=(-Z "$GBZ" -d "$DIST" -m "$MIN" -f "$FQ" -t 32)
  [ -n "${ZIP:-}" ] && ARGS+=(-z "$ZIP")
  $VG giraffe "${ARGS[@]}" > $OUT/${S}.gam
  $VG stats -a $OUT/${S}.gam > $OUT/${S}.stats.txt
  rm -f $OUT/${S}.gam
  grep -E "Total aligned|Total perfect|Total gapless|mapping quality" $OUT/${S}.stats.txt | head -6
done
step "DONE"
