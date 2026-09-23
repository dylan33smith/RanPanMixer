#!/bin/bash
# chr21 preprocessing for PanMixer, re-pointed at the GRCh38 30x panel.
# The pinned checkout is NOT modified: we swap which file sits in the
# `1000g_phased.vcf.gz` slot and drive the shipped .sbatch scripts through the
# local sbatch shim. See docs/bugs.md for why the shipped GRCh37 panel is unusable.
set -euo pipefail
R=/home/ds85/projects/Pangenome/RanPanMixer
P=/data/ds85/RanPanMixer/external/PanMixer
export PATH=$R/scripts:/home/ds85/miniconda3/envs/panmixer/bin:$PATH
export PYTHON=/home/ds85/miniconda3/envs/panmixer/bin/python
export PYTHONNOUSERSITE=1
cd $P/starting_data

step () { echo; echo "######## $* :: $(date +%H:%M:%S)"; }

step "0. put the GRCh38 30x panel in the 1000g_phased slot"
if [ ! -L chr21/1000g_phased.vcf.gz ]; then
  mv -n chr21/1000g_phased.vcf.gz     chr21/DEPRECATED_1000g_phase3_grch37.vcf.gz
  mv -n chr21/1000g_phased.vcf.gz.tbi chr21/DEPRECATED_1000g_phase3_grch37.vcf.gz.tbi
  ln -sf 1000g_30x_phased.vcf.gz      chr21/1000g_phased.vcf.gz
  ln -sf 1000g_30x_phased.vcf.gz.tbi  chr21/1000g_phased.vcf.gz.tbi
fi
ls -la chr21/1000g_phased.vcf.gz

step "1. get_num_alleles (num_alleles.npy)"
sbatch --array=21 scripts/get_num_alleles.sbatch

step "2. get_blocks (plink LD blocks on the GRCh38 panel)"
sbatch --array=21 scripts/get_blocks.sbatch

step "3. convert_2_npy (VCFtoNP x3 -- the slow one)"
sbatch --array=21 scripts/convert_2_npy.sbatch

step "4. get_mappings  <-- THE BUILD CANARY"
sbatch --array=21 scripts/get_mappings.sbatch

step "5. build_biallelic_snp_mask"
sbatch --array=21 scripts/build_biallelic_snp_mask.sbatch

step "6. segment_blocks (blocks_dict.json)"
sbatch --array=21 scripts/segment_blocks.sbatch

step "7. get_af (allele_frequencies.npy)"
sbatch --array=21 scripts/get_af.sbatch

step "8. apply_mask (ORPHAN upstream -- never called by the driver, but the"
step "   gap-score attack needs its output; see docs/bugs.md)"
sbatch --array=21 scripts/apply_mask.sbatch

step "DONE"
ls -la chr21/
