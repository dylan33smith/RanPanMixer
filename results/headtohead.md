# ARM 1 vs PanMixer — ONE BUILD, SAME 44 INDIVIDUALS, chr21

Provenance: scripts/arm1_release.py --ceilings (308 cells, exit 0, 9:38.88) +
PanMixer obfuscate.py on /home/ds85/pm_ourbuild (352 runs, exit 0, 16:13.71, seed 123).
Ours is an EXACT expectation; theirs is ONE DRAW per subject. Denominator: the
two-sided mask on both sides (the headline convention, decision 2026-10-27).

  our endpoint: arm='arm1' scope='chr21' acc_ceilings=['0.50', '0.55', '0.60', '0.65', '0.70', '0.75', '0.80'] n=44

**OURS — ARM 1, chr21, K=86 LOO, EXACT expectation.** Metrics are ROWS, `acc_ceiling` is the COLUMN. n = 44 targets.

| metric | 0.50 | 0.55 | 0.60 | 0.65 | 0.70 | 0.75 | 0.80 |
|---|---|---|---|---|---|---|---|
| `target_fidelity` median | 0.176250 | 0.176095 | 0.175940 | 0.175785 | 0.175631 | 0.175476 | 0.175321 |
|   min | 0.158669 | 0.158496 | 0.158322 | 0.158149 | 0.157975 | 0.157802 | 0.157628 |
|   max | 0.186487 | 0.186337 | 0.186187 | 0.186037 | 0.185887 | 0.185737 | 0.185586 |

**THEIRS — PanMixer on OUR BUILD, chr21, same 44 subjects, seed 123, two-sided denominator.** `capacity` is the COLUMN.

| metric | 0.01 | 0.05 | 0.1 | 0.2 | 0.3 | 0.5 | 0.7 | 0.9 |
|---|---|---|---|---|---|---|---|---|
| `changed_frac` median | 0.011632 | 0.050219 | 0.100303 | 0.173454 | 0.176759 | 0.176790 | 0.176790 | 0.176790 |
|   min | 0.010921 | 0.045560 | 0.085687 | 0.154711 | 0.156076 | 0.156108 | 0.156108 | 0.156108 |
|   max | 0.012536 | 0.057606 | 0.111590 | 0.183902 | 0.188384 | 0.188469 | 0.188469 | 0.188469 |
| `total_alleles` median | 654,211 | 654,196 | 654,072 | 653,020 | 652,880 | 652,880 | 652,880 | 652,880 |
| n | 44 | 44 | 44 | 44 | 44 | 44 | 44 | 44 |

  THEIR SATURATION on our build: 0.176790 (17.6790%) median over 44
  their published / own-build figure at capacity 0.3: 0.164090 (16.4090%)
  our median at acc_ceiling 0.80: 0.175321
  our median at acc_ceiling 0.50: 0.176250  (eta = 0, UNTILTED)

  subjects on BOTH sides: 44 (ours 44, theirs 44)

  BUILD EFFECT, HG00438 at capacity 0.3, like for like:
    their build : 0.164090   (logs, 2026-10-07)
    our build   : 0.163284
    difference  : -0.000806  (-0.0806 points)
  ⚠ this is the ONLY like-for-like build comparison available: their own-build sweep is HG00438 alone, so it must NOT be set against our 44-subject median.

  PAIRED: our `target_fidelity` at acc_ceiling 0.80 vs their `changed_frac` at saturation (capacity 0.9), same subject, same build, n = 44
    we are MORE FAITHFUL (lower) in 27/44 subjects
    median difference : -0.001023  (-0.1023 points)
    mean difference   : -0.001197
    sd of difference  : 0.002437
    range             : -0.005946 .. +0.004346
    sign test (two-sided, exact): p = 1.742e-01
    paired t-like z   : -3.26

  OUR KNOB, acc_ceiling 0.50 -> 0.80 (untilted -> strongest), per subject:
    median movement: +0.000935  (+0.0935 points, toward MORE faithful)
    subjects moving in the faithful direction: 44/44
    range: +0.000877 .. +0.001054

  ⚠ Our 0.50 column is eta = 0: NO protection applied. Its 0.176250 is reached with ZERO target information and is NOT ours to quote as a privacy result.
