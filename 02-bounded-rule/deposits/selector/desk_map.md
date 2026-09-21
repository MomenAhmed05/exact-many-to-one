# Desk map: which rule to use

Selector family: **plug-in iff benchmark df >= t and ratio spread <= s**, else the guaranteed envelope (K <= 3) or an unmeasured cell (K = 5/8).

- random row-level split (seed 20260916), 165 train / 165 test rows
- train-fitted thresholds (`df >= 7`, spread <= 25): train accuracy **0.776**, held-out accuracy **0.752**
- test-fitted thresholds (`df >= 11`, spread <= 25): held-out accuracy **0.764**
- legacy group-disjoint split, same family: train 0.783 / test 0.646 (the number quoted in the first report; NOT a property of the family)
- desk families inside +/-2% under the plug-in: **False** (median over 4 seeds per cell)

## Desk families

| family | bench df | median |plug-in error| | seeds | inside +/-2% |
|---|---|---|---|---|
| A | 3 | 12.32% | 4 | False |
| A | 5 | 0.51% | 4 | True |
| A | 7 | 1.42% | 4 | True |
| A | 11 | 2.03% | 4 | False |
| A | 19 | 1.52% | 4 | True |
| A | 35 | 1.81% | 4 | True |
| B | 3 | 14.53% | 4 | False |
| B | 5 | 1.77% | 4 | True |
| B | 7 | 0.32% | 4 | True |
| B | 11 | 0.61% | 4 | True |
| B | 19 | 0.70% | 4 | True |
| B | 35 | 1.53% | 4 | True |
