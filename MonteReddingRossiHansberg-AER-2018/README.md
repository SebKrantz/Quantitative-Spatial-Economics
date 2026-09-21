# Commuting, Migration, and Local Employment Elasticities

**Ferdinando Monte, Stephen J. Redding & Esteban Rossi-Hansberg** — *American Economic Review* 108(12), December 2018, 3855–3890. [doi:10.1257/aer.20151507](https://doi.org/10.1257/aer.20151507)

*JEL:* J23, J61, R23, R32, R41 · *Keywords:* commuting; migration; local labour markets; gravity; quantitative spatial equilibrium.

Julia port of the MATLAB **Toolkit for Quantitative Spatial Models** by Gabriel M. Ahlfeldt and Tobias Seidel (MIT licence, `LICENSE-MRRH2018-toolkit`), which implements MRRH2018 in the **Seidel & Wickerath (2020, *RSUE* 85)** German-county variant. Upstream checkout: `../Ahlfeldt/MRRH2018-toolkit/` (read-only).

Run it:

```bash
julia MonteReddingRossiHansberg-AER-2018/monte_redding_rossihansberg_commuting_model.jl
```

About 70 seconds. Prints a 21-check validation report and writes six PDFs to `graphs/`.

---

## 1. What the Paper Does

The paper asks **how much local employment moves when a place gets a positive labour-demand shock** — and answers that there is no single number. The elasticity is an *endogenous variable* that differs across locations according to how open their labour market is to commuting.

The mechanism is one modelling choice. A worker draws a Fréchet shock over every **(residence, workplace) pair**, so she chooses where to live and where to work *jointly*. Commuting and migration become the same decision, and a productivity shock in county *i* can raise employment in *i* either by pulling in residents (migration) or by pulling in commuters from neighbours — without anyone moving house, and without bidding up *i*'s land prices.

Three results follow:

1. **Heterogeneity.** Shocking each of 3,111 US counties with a 5% productivity shock one at a time gives an employment elasticity ranging from about **0.5 to 2.5**, mean 1.52. The resident elasticity is tighter, about 0.2 to 1.2. Since employment and residents can only differ through commuting, that gap *is* the commuting channel.
2. **A summary statistic.** Standard controls (size, wages, land area, neighbours) explain only about half the variation. The **residence own-commuting share** $\lambda^R_{ii|i}$ alone gives $R^2 = 0.89$ (Table 2, column 5). Reduced-form work on local labour markets can control for this heterogeneity with one variable.
3. **Aggregates.** Commuting matters for welfare, not just for incidence. A 12% fall in commuting costs — the median change implied by observed 1990–2010 flows — raises welfare by **3.3%**, larger than standard estimates of the gains from trade for the US.

The paper backs this with quasi-experimental evidence independent of the model: winner-versus-runner-up counties in million-dollar-plant competitions show larger employment responses where commuting markets are more open.

---

## 2. The Model (Section I)

Locations $n, i \in N$; $n$ indexes residence/consumption, $i$ indexes workplace/production.

### 2.1 Preferences and the two-margin choice

A worker $\omega$ living in $n$ and working in $i$ has Cobb-Douglas utility over goods and residential land, scaled by an idiosyncratic amenity draw and an iceberg commuting cost:

$$U_{ni\omega} = \frac{b_{ni\omega}}{\kappa_{ni}}\Big(\frac{C_{n\omega}}{\alpha}\Big)^{\alpha}\Big(\frac{H_{n\omega}}{1-\alpha}\Big)^{1-\alpha}, \qquad G_{ni}(b) = e^{-B_{ni}b^{-\epsilon}} \tag{1, 2}$$

Indirect utility is $b_{ni\omega} w_i / (\kappa_{ni}P_n^{\alpha}Q_n^{1-\alpha})$ (eq. 9). With Fréchet draws over all $N^2$ pairs, the unconditional choice probability is a **double** gravity equation:

$$\lambda_{ni} = \frac{B_{ni}\big(\kappa_{ni}P_n^{\alpha}Q_n^{1-\alpha}\big)^{-\epsilon} w_i^{\epsilon}}{\sum_{r}\sum_{s} B_{rs}\big(\kappa_{rs}P_r^{\alpha}Q_r^{1-\alpha}\big)^{-\epsilon} w_s^{\epsilon}} \equiv \frac{\Phi_{ni}}{\Phi} \tag{10}$$

Employment and residents are its two marginals, $L_i = \bar L\sum_n \lambda_{ni}$ and $R_n = \bar L\sum_i \lambda_{ni}$ (eq. 11), and the conditional workplace choice and average residential income are

$$\lambda_{ni|n} = \frac{B_{ni}(w_i/\kappa_{ni})^{\epsilon}}{\sum_s B_{ns}(w_s/\kappa_{ns})^{\epsilon}}, \qquad \bar v_n = \sum_i \lambda_{ni|n} w_i \tag{12, 14}$$

**This second margin is what distinguishes MRRH2018 from everything else in this repository.** Every other model here gives a worker one location.

### 2.2 Goods, trade, land

Monopolistic competition with CES varieties, fixed cost $F$ in labour, iceberg trade cost $d_{ni}$. The trade share and price index are the usual gravity pair:

$$\pi_{ni} = \frac{L_i(d_{ni}w_i/A_i)^{1-\sigma}}{\sum_k L_k(d_{nk}w_k/A_k)^{1-\sigma}}, \qquad P_n = \frac{\sigma}{\sigma-1}\Big(\frac{L_n}{\sigma F \pi_{nn}}\Big)^{\frac{1}{1-\sigma}}\frac{d_{nn}w_n}{A_n} \tag{6, 8}$$

Income equals expenditure and land clears:

$$w_i L_i = \sum_n \pi_{ni}\,\bar v_n R_n, \qquad Q_n = (1-\alpha)\frac{\bar v_n R_n}{H_n} \tag{7, 5}$$

Note where the two margins meet: **expenditure is indexed by residence** ($\bar v_n R_n$) while **the wage bill is indexed by workplace** ($w_i L_i$). Commuting drives a wedge between where income is earned and where it is spent, and that wedge is the whole model.

Equilibrium is six vectors plus a scalar $\{w_n, \bar v_n, Q_n, L_n, R_n, P_n; \bar U\}$ solving eqs. (5)–(8), (11) and labour-market clearing $\bar L = \sum_n R_n = \sum_n L_n$ — seven conditions, exactly identified.

### 2.3 The SW2020 extension (what the toolkit, and this port, actually run)

Seidel & Wickerath add two endogenous objects and two equations:

$$A_i = \bar A_i L_i^{\nu} \quad\text{(agglomeration)}, \qquad H_n = \bar H_n P_{H,n}^{\delta} \quad\text{(housing supply)}$$

so the trade-share numerator picks up $L_i^{1-(1-\sigma)\nu}$. Setting $\nu = \delta = 0$ recovers MRRH2018 exactly, and the port ships both parameter sets.

### 2.4 Quantification: the toolkit's actual contribution

MRRH2018 and SW2020 both quantify from **observed bilateral commuting flows**, rationalising the many zero flows with $\kappa_{ni} = \infty$. That is fine for a level exercise and useless for anything else: an infinite cost is inconsistent with any smooth travel-time matrix, it rules out extensive-margin counterfactuals (a new road turning a zero flow positive), and in most countries the bilateral matrix simply does not exist.

The toolkit's `getBiTK.m` (Codebook **Algorithm 1**) replaces it. Given wages, residence population, workplace employment and *any* smooth cost matrix, iterate

$$\lambda_{ni|n} \propto B_i\,w_i^{\epsilon}\,\kappa_{ni}^{-\epsilon\mu}, \qquad \hat L_i = \sum_n \lambda_{ni|n}R_n, \qquad B_i \leftarrow B_i\cdot\frac{L_i^{\text{obs}}}{\hat L_i}$$

to a fixed point. Workplace amenities $B_i$ absorb whatever the cost matrix cannot explain, and the model is quantified from **employment by workplace, employment by residence, a wage or rent index, land area, and a bilateral cost matrix** — no flow matrix. *That* is why this model is worth porting, and it is the port's primary entry point.

If wages are also unobserved, feed $w = \mathbf{1}$: the recovered $B_i$ is then an ARSW2015-style transformed wage and $w_i = B_i^{1/\epsilon}$. The port exposes this as `wages = :none`.

### 2.5 Counterfactuals in changes

Exact hat algebra ($\hat x = x'/x$) maps changes in fundamentals $\{\hat{\bar A}_n, \hat B_{ni}, \hat\kappa_{ni}, \hat d_{ni}\}$ into changes in everything endogenous, **without the levels of the unobserved fundamentals**. Eight equations, Codebook Algorithms 3–10, solved by one outer fixed point on $\{\hat w, \hat\lambda\}$ (Algorithm 11). Welfare collapses to

$$\hat{\bar U} = \hat\Phi^{1/\epsilon} = \hat B_{ni}^{1/\epsilon}\big(\hat\kappa_{ni}\hat P_{Q,n}^{\alpha}\hat P_{H,n}^{1-\alpha}\big)^{-1}\hat w_i\,\hat\lambda_{ni}^{-1/\epsilon}$$

which is **location-invariant by construction** — free mobility equalises expected utility, so all $J^2$ cells must agree. The port uses their spread as a free convergence diagnostic (it comes out at $3\times10^{-8}$).

---

## 3. Calibration

| Symbol | Code | County / SW2020 (default) | MRRH2018's own US values | Role |
|---|---|---|---|---|
| $\alpha$ | `alpha` | 0.70 | 0.60 | expenditure share on tradables ($1-\alpha$ on housing) |
| $\epsilon$ | `epsilon` | 4.60 | 3.30 | Fréchet shape over (residence, workplace) pairs |
| $\mu$ | `mu` | 0.47 | 1.3424 (so $\mu\epsilon = \phi = 4.43$) | distance/travel-time elasticity of $\kappa$ |
| $\delta$ | `delta` | 0.38 | 0 (inelastic land) | housing supply elasticity |
| $\sigma$ | `sigma` | 4 | 4 | CES across varieties (Broda–Weinstein) |
| $F$ | `fixC` | 1 | 1 | fixed production cost, labour units |
| $\nu$ | `nu` | 0.05 | 0 | agglomeration, $A_i = \bar A_i L_i^{\nu}$ |
| $\psi$ | `psi` | 0.42 | 0.43 | distance elasticity of $d_{ni}$ |
| $J$ | — | 401 Kreise | 3,111 counties | locations |

Two composites are the ones actually disciplined by data:

* **Trade gravity** $\psi(1-\sigma) = -1.26$, against MRRH2018's CFS estimate of $-1.29$.
* **Commuting gravity** $-\epsilon\mu = -2.162$, against MRRH2018's US estimate of $-\phi = -4.43$. The observed German flows give $-2.087$ with residence and workplace fixed effects, implying $\mu = 0.454$ — a close match to the toolkit's calibrated 0.47, and the only check on $\mu$ that this data supports.

Commuting is far more local than trade in both countries, which is the paper's point about moving people versus moving goods.

---

## 4. Results

### 4.1 Table 5 — the headline welfare experiment

The paper's shock is on the **relative** ease of commuting (its Head–Ries measure, eq. 23, normalises by own commuting), so $\hat\kappa$ applies to off-diagonal pairs only and the diagonal stays at 1.

| $\hat\kappa$ (off-diagonal) | MRRH2018 (US) | OwnData (DE) | ReadData (DE) | ReadData (DE), US params |
|---|---|---|---|---|
| 0.79 (p75 reduction) | +6.89% | +5.22% | +11.34% | +10.18% |
| **0.88 (median, −12%)** | **+3.26%** | **+2.32%** | **+5.18%** | **+4.85%** |
| 0.96 (p25 reduction) | +0.89% | +0.63% | +1.43% | +1.39% |
| 1.13 (increase) | −2.33% | −1.40% | −3.27% | −3.42% |

**The level is not reproduced and was not expected to be** — this is German data, not US. The two German tracks *bracket* the paper at every percentile, and the signs, ordering and rough proportions all hold. Switching the ReadData track to MRRH2018's own parameters moves +5.18% only to +4.85%, so the gap is **data, not parameters**. The direction of each track is understood:

* **ReadData lands above** because German counties are slightly *more* open to commuting than US ones — median $\lambda_{nn|n}$ of 0.643 against the US 0.69 in 2000 — so a proportional cut in off-diagonal costs buys more.
* **OwnData lands below** because its *predicted* flows are far too concentrated on the diagonal (median 0.837 against the observed 0.643): straight-line distance with a very small within-county diagonal makes counties look more self-contained than they are. A user with real travel times should pass them to `build_costs`.

The residual is general-equilibrium dampening, which depends on the number and size distribution of locations (401 Kreise against 3,111 counties) and cannot be decomposed further without MRRH2018's own data.

### 4.2 Section III — local employment elasticities

Shocking counties one at a time with the paper's 5% productivity shock, on a 50-county subsample:

| | employment elasticity | resident elasticity | $R^2$ on $\lambda^R_{ii\mid i}$ |
|---|---|---|---|
| MRRH2018 (US, 3,111 counties) | [0.50, 2.50], mean 1.52 | [0.20, 1.20] | 0.890 |
| **German counties, MRRH2018 parameters** | **[1.12, 2.18], mean 1.79** | **[0.47, 1.05], mean 0.73** | **0.851** |
| German counties, toolkit parameters | [2.82, 4.20], mean 3.72 | [1.26, 3.60], mean 2.11 | 0.277 |

With the paper's own parameters this is a close quantitative reproduction on different data: the elasticity is heterogeneous, employment is more dispersed than residents, the slope on the own-commuting share is negative (more open → larger response), and one variable explains 85% of the variation against the paper's 89%.

The toolkit row shows what the SW2020 parameters do: $\delta = 0.38$ (elastic housing), $\epsilon = 4.6$, $\nu = 0.05$ and $\alpha = 0.70$ all weaken the congestion force or strengthen mobility, and all push the same way. The $R^2$ collapses to 0.28 because with elastic housing the housing-supply margin competes with commuting as a source of heterogeneity — which is, in miniature, the paper's own Section III.B point about Saiz elasticities.

### 4.3 Didactic counterfactual — a new inner-German border

The toolkit's teaching exercise (not a result of either paper): multiply the cost of every East–West route by 1,000.

| Shock | Welfare |
|---|---|
| Trade border only | −2.12% |
| Commuting border only | −0.18% |
| Both | −2.24% |

Residents fall sharply in the East (mean log change −0.55) and rise slightly in the West. A commuting border alone costs an order of magnitude less than a trade border, because only a thin strip of counties actually commutes across the line, whereas every county trades across it. Halving the trade-cost elasticity to $\psi = 0.21$ more than doubles the welfare loss (−5.17%), since cheaper trade means more of it to lose.

---

## 5. Validation

The upstream toolkit ships **no saved outputs** (`data/output/` holds only `Dummy.txt`), so the port is validated against the paper and against internal identities. All 21 checks pass on every run; the script prints them.

| Check | Target | Achieved |
|---|---|---|
| V1 Conditional commuting probabilities sum to 1 per residence | exact | 2.7e-15 |
| V2 Unconditional $\lambda_{ni}$ sums to 1 | exact | 0.0 |
| V3 Labour market clearing $\sum L = \sum R = \bar L$ | exact | 1.4e-16 |
| V4 Trade shares sum to 1 per destination | exact | 2.4e-15 |
| V5 Algorithm 1 matches observed workplace employment | upstream rule < 1e-3 | 9.9e-4, max rel dev 7.9e-8 |
| V6 Algorithm 1 fixed point is damping-invariant | — | 1.6e-8 between $\zeta$ = 1.0 and 0.5 |
| V7 Productivity inversion: income == expenditure | 6-dp rounding rule | 4.7e-7, 71 iterations |
| **V8 Levels round-trip from a cold start** | ~1e-6 | **1.4e-7** |
| V8b Local uniqueness under a 50% log-normal perturbation | — | 3.9e-8, returns to the same point |
| **V9 Zero-shock exact hat: all hats == 1** | exact | **2.6e-8, welfare +0.0000%** |
| V10 Welfare is location-invariant across all $J^2$ cells | — | spread 3.2e-8 |
| **V11 Model commuting gravity == $-\epsilon\mu$** | −2.162 | **−2.16200000, diff 0.0** |
| **V12 Model trade gravity == $\psi(1-\sigma)$** | −1.26 | **−1.26000000, diff 2e-16** |
| V13 Observed German commuting gravity | — | −2.087, implies $\mu$ = 0.454 vs calibrated 0.47 |
| V14 Structural trade gravity vs MRRH2018's estimate | −1.29 | −1.26 |
| **V15 Headline welfare, $\hat\kappa$ = 0.88** | +3.26% | **+2.32% / +5.18%, brackets the paper** |
| V16/V17 Full Table 5: signs, ordering, monotonicity | — | all match |
| V18 Uniform commuting shock vs closed form $\hat U = 1/\hat\kappa$ | 1.13636364 | 1.13636364, diff 2e-10 |
| **V19 Local employment elasticity heterogeneity** | [0.5, 2.5] | **[1.12, 2.18] with the paper's parameters** |
| **V20 $R^2$ on the own-commuting share** | 0.89 | **0.851** |

Two of these deserve a word. **V11/V12 are identities, not estimates**: model-generated $\log\lambda_{ni}$ is exactly a residence effect plus a workplace effect plus $-\epsilon\mu\log\text{dist}$, so the two-way fixed-effect slope *must* return $-\epsilon\mu$ to machine precision. Anything else would be a transposed matrix. They are the cheapest possible test that the $[n,i]$ orientation is right everywhere. And **V8 is a genuine test, not a tautology**: the data is an exact fixed point by construction, but finding its way back from a uniform cold start through 521 iterations of a strongly non-linear system is not.

---

## 6. Known Upstream Defects, and What This Port Does About Each

| # | Defect | Location | Handling |
|---|---|---|---|
| 1 | $\psi$ switched 0.42 → 0.21 mid-script | `Counterfactuals.m:139`, restored at `:221` | **Not a bug.** Reading the whole script, it is a bracketed "half the trade cost" sensitivity run. $\psi = 0.42$ is the baseline; $\psi = 0.21$ ships as `LOW_TRADECOST_PARAMS`. The interruption hazard (an aborted run leaves `data/output/*.mat` in the low state) is commented. |
| 2 | GRID track uses $\nu = 0$, county track $\nu = 0.05$, silently | `GRID_MRRH2018_toolkit.m:55` | Both named: `COUNTY_PARAMS` and `GRID_PARAMS`. Neither is inherited. The headline counterfactual is reported under both (+2.32% vs +2.29%). |
| 3 | `getBiTK.m` takes the full step despite Codebook Algorithm 1 step 5 specifying $\zeta < 1$ | `getBiTK.m:70` | `damp` keyword added, default 1.0 (exact upstream reproduction). V6 verifies $\zeta = 0.5$ reaches the same fixed point to 1.6e-8. |
| 4 | `counterFactsTK.m` outer loop is `while true` — an extreme shock hangs MATLAB | `counterFactsTK.m:87` | `maxiter` (default 5,000) plus `@warn` and return of the last iterate. |
| 5 | `solveProductTradeTK.m` returns the **transposed** trade-share matrix | `:93, :99` vs `:66, :69` | Invisible only because $d_{ni}$ is symmetric — the convergence block applies `dni.^(1-σ)` without transposing it, so an asymmetric cost matrix would silently transpose the model. The port uses one orientation $\pi[n,i]$ throughout and asserts the row sums. |
| 6 | Normalisations commented out in three updaters | `updateEmplTK.m:28-31` and others | Left out, matching shipped behaviour. Noted in comments. |
| **7** | **`solveProductTradeTK.m:107` drops both the $L_n^{-\nu}$ agglomeration correction and $d_{nn}$ from the price index** | `:107` | **Found in this port, not previously documented.** Understates $P_n$ by up to a factor of 2.9 on this data. Harmless upstream — $P_n$ is only ever mapped in `Descriptives.m`, never used in a counterfactual, and `updatePricesTK.m:25` has both terms and is correct. **Not** harmless here, because $P_n$ enters the forward solver through $\lambda_{ni}$. The correct expression is used; `price_index_upstream` is kept for the comparison the script prints. |

---

## 7. Relation to the Other Models in This Repository

| Model | Shares | Differs |
|---|---|---|
| `Redding-JIE-2016` | CES gravity $\pi_{ni}\propto\text{cost}^{1-\sigma}$, imperfect labour mobility | Migration only — one location per worker. Porting MRRH2018 into it means adding a second, nested Fréchet margin. |
| `QSE-ARE-2017` | Helpman (1998): tradables + non-tradable housing, migration | No commuting at all. A simpler ancestor that MRRH2018 nests in spirit. |
| `QRE-HoRaUE-2025` | the gravity/Fréchet DNA | The Handbook workhorse nests eight models through **one** location-choice margin and four elasticities. **MRRH2018 is not among them and cannot be** — a second Fréchet draw over workplaces is a genuine extra dimension, not a corner case. Worth flagging for anyone extending that file. |
| `AllenArkolakis-RES-2022`, `FuchsFoongWong-MMN-2026` | exact-hat counterfactuals, damped fixed points | Those **endogenise** transport costs (route choice, Leontief-inverse congestion, multimodal CES). Here $\kappa_{ni}$ and $d_{ni}$ are exogenous icebergs and a counterfactual supplies $\hat\kappa/\hat d$ from outside. Bolting AA2022's route-choice block onto this model is the obvious merge. |
| `BrockhausHinzSerfaty-2026` | multi-region trade, exact hat algebra | Multi-sector Caliendo–Parro with endogenous route/mode/toll choice; no labour mobility or commuting. |
| `Ahlfeldt/ARSW2015-toolkit` | the residence × workplace Fréchet engine itself | Single metro area, no inter-regional goods trade, plus endogenous building height and a land-use choice. MRRH2018 drops the vertical margin and adds multi-region trade. Complementary, not overlapping. |

In one line: **this is the only model in the repository with a commuting margin, and the only one where a worker chooses two locations.**

### Notation collisions

The same Greek letters mean different things across these models. The file header carries the full warning; the two dangerous ones:

* **$\lambda$** is a commuting *probability* here and the traffic *congestion elasticity* (0.07–0.09) in `AllenArkolakis-RES-2022`.
* **$\alpha$** is the expenditure share on tradables (0.70) here and the *productivity externality* (0.1) in AA2022 — different kinds of parameter entirely.

Also: $\epsilon = 4.6$ here is a Fréchet shape over commuting pairs, not the migration elasticity 3–5 used elsewhere; and $\pi_{ni}$ is a *trade* share here but a *commuting* probability in ARSW2015.

---

## 8. Contents of This Folder

| Path | What |
|---|---|
| `monte_redding_rossihansberg_commuting_model.jl` | **The port.** Single self-contained file, 1,325 lines, runs with `julia <file>.jl`. |
| `TRANSLATION_PLAN.md` | Equation-by-equation mapping, parameter table, solver architecture, validation targets, every decision taken, and the post-implementation record. |
| `graphs/` | Six PDFs written by the run (see below). |
| `matlab_source/` | All of upstream's `progs/*.m` and `scripts/*.m`. `progs/old/` is deliberately excluded. |
| `reference_data/` | The 401-German-county input CSVs. |
| `MRRH2018-paper.pdf`, `MinerU_markdown_MRRH2018-paper.md` | The AER paper. |
| `MRRH2018-appendix.pdf`, `MinerU_markdown_MRRH2018-appendix.md` | Online appendix — the derivations live here. |
| `MRRH2018-codebook.pdf` | Ahlfeldt & Seidel's codebook: the full equilibrium and counterfactual systems plus Algorithms 1–11 in pseudo-code. Read this first. |
| `LICENSE-MRRH2018-toolkit` | Upstream MIT licence. |

Figures: `fit_getbi.pdf` (Algorithm 1 fit and own-commuting shares), `gravity.pdf` (commuting and trade gravity, model vs observed), `fundamentals.pdf` (inverted $\bar A$, $B_i$, CMA, trade openness), `table5_commuting_costs.pdf` (the headline experiment against the paper's points), `border_counterfactual.pdf` (border-distance discontinuities), `local_elasticities.pdf` (elasticity distribution and its relation to commuting openness).

Upstream's choropleths are not reproduced: `MAPIT.m` needs the MATLAB Mapping Toolbox and the `VG250_KRS_clean_final` shapefile, neither of which is vendored here.

### Using this with other data

The OwnData pathway is the point. Replace exactly five objects in `load_county_data`:

| Object | Shape | Note |
|---|---|---|
| `dist` | $J\times J$ | distance or travel time; **strictly positive diagonal** (an internal distance) |
| `L_obs` | $J$ | employment by **workplace** |
| `R_obs` | $J$ | employment by **residence** |
| `w_obs` | $J$ | wages by workplace, or pass `wages = :none` |
| `area` | $J$ | land area |

Nothing else changes. The bilateral commuting matrix is needed only for the ReadData comparison track and the descriptive gravity regressions.

---

## 9. Key References

* Monte, F., S. J. Redding & E. Rossi-Hansberg (2018). "Commuting, Migration, and Local Employment Elasticities." *American Economic Review* 108(12), 3855–3890.
* Seidel, T. & J. Wickerath (2020). "Rush hours and urbanization." *Regional Science and Urban Economics* 85, 103580. — the German county application this toolkit's data and code descend from.
* Ahlfeldt, G. M. & T. Seidel (2024). *Codebook for: Toolkit for Quantitative Spatial Models.* — `MRRH2018-codebook.pdf`.
* Ahlfeldt, G. M., S. J. Redding, D. M. Sturm & N. Wolf (2015). "The Economics of Density: Evidence from the Berlin Wall." *Econometrica* 83(6), 2127–2189. — the other residence × workplace Fréchet model, vendored at `../Ahlfeldt/ARSW2015-toolkit`.
* Head, K. & J. Ries (2001). "Increasing Returns versus National Product Differentiation." *American Economic Review* 91(4), 858–876. — the relative-flows inversion behind eq. (23).
* Broda, C. & D. E. Weinstein (2006). "Globalization and the Gains from Variety." *Quarterly Journal of Economics* 121(2), 541–585. — the source of $\sigma = 4$.
* Greenstone, M., R. Hornbeck & E. Moretti (2010). "Identifying Agglomeration Spillovers." *Journal of Political Economy* 118(3), 536–598. — the million-dollar-plant design used in Section IV.
