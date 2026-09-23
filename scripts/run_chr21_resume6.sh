#!/bin/bash
set -euo pipefail
R=/home/ds85/projects/Pangenome/RanPanMixer
P=/data/ds85/RanPanMixer/external/PanMixer
export PATH=$R/scripts:/home/ds85/miniconda3/envs/panmixer/bin:$PATH
export PYTHON=/home/ds85/miniconda3/envs/panmixer/bin/python
export PYTHONNOUSERSITE=1
cd $P/starting_data
step () { echo; echo "######## $* :: $(date +%H:%M:%S)"; }
step "6. segment_blocks (blocks_dict.json)"; sbatch --array=21 scripts/segment_blocks.sbatch
step "7. get_af (allele_frequencies.npy)";    sbatch --array=21 scripts/get_af.sbatch
step "8. apply_mask (orphan upstream)";       sbatch --array=21 scripts/apply_mask.sbatch
step "DONE"
