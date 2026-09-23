set -euo pipefail
export PATH=/home/ds85/miniconda3/envs/panmixer/bin:$PATH
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
cd /data/ds85/RanPanMixer/external/PanMixer
OUT=/data/ds85/RanPanMixer/runs/obf_test
mkdir -p $OUT
# determinism check: same subject, same capacity, same seed, twice
for tag in seedA seedArepeat; do
  echo "### $tag  $(date +%H:%M:%S)"
  python tools/panmixer/obfuscate.py starting_data HG00438 0.1 21 $OUT/$tag --seed 123
done
echo "### seedB (different seed)  $(date +%H:%M:%S)"
python tools/panmixer/obfuscate.py starting_data HG00438 0.1 21 $OUT/seedB --seed 456
echo "### cap0.5  $(date +%H:%M:%S)"
python tools/panmixer/obfuscate.py starting_data HG00438 0.5 21 $OUT/cap0p5 --seed 123
echo "### DONE $(date +%H:%M:%S)"
