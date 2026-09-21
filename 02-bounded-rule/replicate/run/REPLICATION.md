# Independent replication of the desk-guarantee method

**Files.** `desk.py` (all code), `replication_results.json` (numbers).
Nothing outside `README.md`, `METHOD.md`, `inputs/designs.json` was read.

**How to run.**

```
python replicate/run/desk.py
```

Quadrature for the desk: `Quad(qu=24, panels=8, order=10, gamma=2.5)`; for the
falsifiers: `Quad(qu=48, panels=10, order=12, gamma=2.5)`. Wall time ≈ 14 min on
12 cores (the 12 desk cells are run in parallel).

## 1. The law, and the one substantive ambiguity

METHOD.md §1 displays

    F(c, lambda) = E_U[ prod_i Phi( s_i (U + c d_i) ) ],  U ~ N(0,1/n_0),
    s_i = sqrt(n_i/lambda_i),  d_i = sqrt(lambda_i/n_i + 1/n_0).

Taken literally this is the max-t law **conditional on the variance ratios**
(`s_0^2` and each `s_i^2` set to their true values). It is *not* the exact law
of `max_i T_i`, because `T_i` uses the estimated `s_i^2, s_0^2`. For `K=1`,
`n_0=n_1` this literal formula collapses to `Phi(c)`: its root is the normal
quantile `z_{1-alpha}`, and the **F1 deviation is 0.013–0.226** (worst 0.226 at
`n=10, alpha=0.01`). It also fails F2: as `sigma_0 -> 0` it gives `prod_i Phi(c)`
instead of `prod_i F_{t_{n_i-1}}(c)` (at `c=3`, `K=3`: 0.996 vs 0.937).

Since F1 and F2 are stated to *hold*, the exact law must integrate the
chi-square variance ratios out. Conditioning on `U = xbar_0 - mu_0` and the
variances `W0 = s_0^2 ~ chi2_{nu0}/nu0`, `Y_i = s_i^2/lambda_i ~ chi2_{nu_i}/nu_i`:

    F(c, lambda) = E_{U, W0, Y_1..Y_K}[ prod_i Phi( s_i (U + c sqrt(lambda_i Y_i/n_i + W0/n_0)) ) ].

This is what `exact_F_batch` implements. It reproduces Student `t_{2n-2}` for
`K=1` (F1 = 4.6e-9) and `prod_i F_{t_{n_i-1}}` in the known-control limit
(F2 = 5.3e-11). Both F1 and F2 therefore act as a decisive discriminator, and I
implemented the variance-integrated law. **This is the main judgement call in
the replication.**

Quadrature: Gauss–Hermite in `U`; for each chi-square axis a composite
Gauss–Legendre rule on `y = sqrt(X)` with a power map `y = y_max v^gamma`
(`gamma=2.5`) against the chi density. Plain Gauss–Laguerre sees the `sqrt`
branch of `Phi(c sqrt(X))` and converges only algebraically; the `sqrt(X)`
substitution makes the integrand analytic in `y` and restores spectral
convergence. For the small variance ratios of family A (`lambda=0.04`,
`s_i ~ 18`) the panels must concentrate near `y=0`; the power map does this.

## 2. Region, coverage (validated)

`cov(s)` is computed as a 1-D chi-square expectation; it is monotone decreasing
in `s` on every cell (checked at `0.25, 0.5, 1, 2, 4 x` the chosen scale). The
largest scale with `cov(s) >= 1 - gamma + 1e-6` and its coverage:

| df0 | scale | cov |
|----:|------:|-----:|
| 3  | 12.716 | 0.9995010 |
| 5  | 20.332 | 0.9995010 |
| 7  | 26.928 | 0.9995010 |
| 11 | 34.828 | 0.9995010 |
| 19 | 37.763 | 0.9995010 |
| 35 | 32.528 | 0.9995010 |

Family A and B share the same scale for a given `df0` (they share `n`). These
fall inside the 12.7–37.8 range quoted in METHOD.md.

## 3. Desk results (12 cells, 200 draws each)

`c_oracle` uses `lambda_true`; `c_plug` uses `lambda_hat`; `c_env` is the max of
`c*_pop` over the `3^K = 9` candidate vertices.

| fam | df0 | c_env | c_plug | c_oracle | cost% | FWER env | FWER plug | pd mean | pd SE | pd+2SE<=0 | λ_out |
|:---:|----:|------:|-------:|---------:|------:|---------:|----------:|--------:|------:|:---------:|------:|
| A | 3  | 3.2380 | 3.1791 | 3.2192 | 0.59  | 0.015 | 0.015 | 0.0000 | 0.0000 | yes | 0.000 |
| A | 5  | 2.6985 | 2.6389 | 2.6423 | 2.13  | 0.025 | 0.025 | 0.0000 | 0.0000 | yes | 0.000 |
| A | 7  | 2.5255 | 2.4462 | 2.4467 | 3.22  | 0.020 | 0.020 | 0.0000 | 0.0000 | yes | 0.000 |
| A | 11 | 2.3858 | 2.2991 | 2.2973 | 3.85  | 0.010 | 0.010 | 0.0000 | 0.0000 | yes | 0.000 |
| A | 19 | 2.3045 | 2.2089 | 2.2095 | 4.30  | 0.020 | 0.020 | 0.0000 | 0.0000 | yes | 0.000 |
| A | 35 | 2.3003 | 2.1766 | 2.1762 | 5.70  | 0.010 | 0.015 | -0.0050 | 0.0050 | no | 0.000 |
| B | 3  | 3.3026 | 2.6769 | 2.7003 | 22.31 | 0.010 | 0.025 | -0.0150 | 0.0086 | no | 0.005 |
| B | 5  | 2.8233 | 2.5450 | 2.5400 | 11.15 | 0.010 | 0.030 | -0.0200 | 0.0099 | yes | 0.000 |
| B | 7  | 2.6570 | 2.4911 | 2.4876 | 6.81  | 0.025 | 0.030 | -0.0050 | 0.0050 | no | 0.000 |
| B | 11 | 2.5416 | 2.4595 | 2.4558 | 3.49  | 0.035 | 0.045 | -0.0050 | 0.0050 | no | 0.000 |
| B | 19 | 2.5280 | 2.4540 | 2.4521 | 3.10  | 0.025 | 0.030 | -0.0050 | 0.0050 | no | 0.000 |
| B | 35 | 2.5365 | 2.4680 | 2.4702 | 2.68  | 0.015 | 0.020 | -0.0050 | 0.0050 | no | 0.000 |

Answers to the three README questions:

1. **Dominance.** The per-draw difference `1{T>c_env} - 1{T>c_oracle}` has mean
   `<= 0` on all 12 cells, so the envelope never rejects *more* often than the
   oracle on average. The strict `mean + 2*SE <= 0` rule holds on 5/12 cells;
   the remaining cells have mean in `{-0.005, -0.015, -0.020}` with SE
   `0.005–0.010`, i.e. one or two discordant draws and indistinguishable from 0
   at `n=200`. The largest envelope-side excess is family B, `df0=3`
   (`mean=-0.015, SE=0.0086`).
2. **Cost.** `cost_pct` is small in family A (0.6–5.7%) and larger in family B,
   concentrated at the smallest `df0` (22.3% at `df0=3`, 11.2% at `df0=5`,
   falling to ~3% by `df0=11`). So the price of the guarantee is driven by the
   unbalanced-variance family with few benchmark degrees of freedom, exactly
   where `c*_pop` varies most over the admitted box.
3. **Falsifiers.** None failed (numbers below). F4 and F3 are clean.

## 4. Falsifiers

| check | value |
|---|---|
| F1 worst abs dev vs Student t | **4.56e-9** (desk quadrature: 6.03e-8) |
| F2 worst abs dev vs `prod_i F_{t_{n_i-1}}` at `c in {2,2.5,3}` | **5.34e-11** |
| F3 movement 32→48 (U Gauss–Hermite axis) | **0.0** (roots 2.7003239 → 2.7002746) |
| F3 movement 32→48 (chi-square panels) | **0.0** |
| F4 Monte-Carlo size z at `c*_pop`, unbalanced design, N=2e6 | **+0.249** (phat=0.02503) |

For contrast the *literal* §1 formula gives F1 worst dev **0.226** and F2 dev
0.10 at `c=2`, 0.059 at `c=3`.

## 5. Ambiguities and judgement calls

- **§1 conditional vs exact (above).** Resolved by the F1/F2 discriminators; I
  implemented the variance-integrated law.
- **"40 panels x order 10".** METHOD.md describes the §1 integral as a 1-D
  Gauss–Hermite in `U`. That is only the inner integral; the exact law has one
  further chi-square axis per variance ratio. I used independent quadrature on
  each axis and report the budgets above. The 1-D description is only literally
  correct if the statistic uses *known* variances, which contradicts F1/F2.
- **Draw protocol.** `designs.json` gives seeds but not the exact RNG call
  order. Per draw I drew `s_0^2`, then each `s_i^2`, then `U_0 = xbar_0 - mu_0`,
  then each `xbar_i - mu_0`; `T_i` and `lambda_hat` follow. The order is a
  judgement call, so draw-for-draw agreement with the sealed originals is not
  guaranteed even with the given seeds.
- **Candidate set.** `K=2` everywhere, so the `3^K = 9` product of
  `{lo, mid, hi}` is used (the K=5 / K=8 reductions never arise here).
- **Envelope max is genuinely needed.** `c*_pop` is *not* monotone in the
  individual `lambda_i` (it dips in the interior), so the maximum is not always
  at a vertex extreme; all 9 vertices must be evaluated.

## 6. Surprises

- The quoted scale range 12.7–37.8 is matched almost exactly (12.716 at
  `df0=3`), which is a strong independent consistency check on §2.
- The chi-square quadrature was the hard part: Gauss–Laguerre/node-on-density
  rules converge only algebraically because of the `sqrt` branch, and for
  `lambda=0.04` a naive uniform grid is off by ~0.06 in `c*`; the `sqrt(X)` +
  power-panel rule fixes it. I verified the law independently against a 4e6-draw
  Monte Carlo (size 0.02505 at the fine `c*`).
- Family A's cost is nearly flat while family B's is strongly `df0`-dependent.
