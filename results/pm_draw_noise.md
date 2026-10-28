  cells: 352 (expect 352)   draws per cell: [5]

**THEIR DRAW NOISE, 5 seeds per cell, n = 44 subjects per capacity.**

| metric | 0.01 | 0.05 | 0.1 | 0.2 | 0.3 | 0.5 | 0.7 | 0.9 |
|---|---|---|---|---|---|---|---|---|
| within-cell sd, median | 6.80e-06 | 5.25e-05 | 1.20e-04 | 7.96e-04 | 9.86e-04 | 9.89e-04 | 9.89e-04 | 9.89e-04 |
| within-cell sd, max | 3.54e-05 | 1.26e-04 | 2.58e-04 | 1.91e-03 | 1.95e-03 | 1.95e-03 | 1.95e-03 | 1.95e-03 |
| between-subject sd | 3.24e-04 | 2.95e-03 | 5.53e-03 | 8.88e-03 | 1.01e-02 | 1.01e-02 | 1.01e-02 | 1.01e-02 |

  AT SATURATION: within-cell (draw) sd 9.888e-04 vs between-subject sd 1.007e-02 -- draw noise is 9.8% of the biological spread
  SEM of a 5-draw mean: 4.422e-04 (one draw: 9.888e-04)

  PAIRED vs ONE DRAW (seed 123):
    more faithful in 27/44   median -0.001023  mean -0.001197  sd 0.002437
    sign test p = 1.742e-01      paired z = -3.26

  PAIRED vs MEAN of 5 DRAWS:
    more faithful in 33/44   median -0.000972  mean -0.001528  sd 0.001929
    sign test p = 1.260e-03      paired z = -5.25
