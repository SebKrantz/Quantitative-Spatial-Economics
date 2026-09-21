# Ahlfeldt & Barr (2022): The Economics of Skyscrapers

Julia port of the MIT-licensed **AB2022 Stata toolkit** by Gabriel M. Ahlfeldt.

> Ahlfeldt, Gabriel M. and Jason Barr (2022), "The economics of skyscrapers: A synthesis",
> *Journal of Urban Economics* 129, 103419. https://doi.org/10.1016/j.jue.2021.103419

```bash
julia ahlfeldt_barr_skyscraper_model.jl
```

Runs in a few seconds, prints the full validation suite against the toolkit's own Stata output,
and writes twelve PDFs to `graphs/`.

**Why this model is here:** it is the only model in this repository with an **endogenous
building-height margin**. Every other model — `AllenArkolakis-RES-2022`, `FuchsFoongWong-MMN-2026`,
`Santamaria-2022`, `QRE-HoRaUE-2025`, `Redding-JIE-2016`, `QSE-ARE-2017`,
`BrockhausHinzSerfaty-2026` — treats floor/land supply as fixed, as an exogenous elasticity, or not
at all, and every location is a point with only a horizontal trade/migration problem. See
[§9](#9-relation-to-the-other-models-in-this-repository).

---

## 1. What the Model Does

A stylized, **perfectly-open, one-dimensional monocentric city** — Alonso–Mills–Muth with an
endogenous *vertical* margin bolted on.

Space is a line `x ∈ [-50, 50]` km on a 0.01 km grid (10,001 points), `D = |x|`. Two competing land
uses, commercial (`C`) and residential (`R`), plus an agricultural outside option at exogenous rent
`r_a`. Developers choose building height to maximise profit per unit land given the floor-space
price they can charge; land goes to whichever use bids the highest land rent; a single city-wide
wage `y` and total employment `N` adjust until commercial labour demand equals residential labour
supply.

**There is no commuting-cost term.** Distance enters *only* through the exponential decay
`exp(-τ^U·D)` of the productivity/amenity fundamentals, and one wage applies everywhere. `τ^C`/`τ^R`
are a reduced-form decay of agglomeration/amenity value with distance, not a commuting primitive.
The Codebook lists a reservation utility `Ū` among the primitives, but it appears nowhere in the
toolkit code and is not one of its 15 solver arguments — openness is implemented by making `N` a
target object rather than fixing utility, so `Ū` is absorbed into `a^R(x)`. That is a simplification
of the toolkit, not a defect.

---

## 2. Model (Codebook §A.1, the "13 equations")

`U ∈ {C, R}`. Shifters (Codebook Table 1):

```
Ã^U(x) = ā^U · ã^U · N^{β^U} · exp(-τ^U·|x|)
```

`ã^U` is the location-specific amenity component — uniformly 1 for the stylized city, and the object
the Chicago inversion (§6.5) solves for. `β^C > 0` is an **agglomeration** elasticity: local
productivity rises with city size.

| # | Equation | |
|---|---|---|
| 1 | `a^C(x) = Ã^C(x)^{1/(1-α^C)} · y^{α^C/(α^C-1)}` | commercial floor-price shifter |
| 2 | `a^R(x) = Ã^R(x)^{1/(1-α^R)} · y^{1/(1-α^R)}` | residential floor-price shifter |
| 3–4 | `S̃^U(x) = min( S*^U(x), S̄^U )`, where `S*^U = (a^U/(c^U(1+θ^U)))^{1/(θ^U-ω^U)}` | constrained height |
| 5–6 | `r^U(x) = a^U/(1+ω^U)·(S̃^U)^{1+ω^U} − c^U·(S̃^U)^{1+θ^U}` | land bid rent |
| 7–8 | land goes to the highest bidder among `r^C`, `r^R`, `r_a` | highest-and-best use |
| 9–10 | `p̄^U(x) = a^U(x)/(1+ω^U) · S^U(x)^{ω^U}` | horizontal (per-floor) rent |
| 11 | `L(x) = α^C/(1-α^C) · p̄^C(x)/y · S^C(x)` | workplace employment (MRS) |
| 12 | `n(x) = S^R(x)/f̄^R(x)`, `f̄^R = (1-α^R)/p̄^R · y` | residence employment (Marshallian) |
| 13 | `∫_{-x0}^{x0} L(x)dx = ∫_{-x1}^{-x0} n(x)dx + ∫_{x0}^{x1} n(x)dx = N` | labour-market clearing |

The developer's problem ties eqs. 3–6 and 9–10 together. Maximising `p̄(S)·S − c·S^{1+θ}` under
`p̄ = a/(1+ω)·S^ω` gives revenue `a/(1+ω)·S^{1+ω}`, hence the FOC `a·S^ω = c(1+θ)·S^θ` and
`S* = (a/(c(1+θ)))^{1/(θ-ω)}` — exactly the `S*` above, and exactly the first term of the land rent.
This internal consistency is what settles the `ω` discrepancy in §5.

`x0` is the outer margin of the commercial zone (CBD radius), `x1` the outer margin of the urban
area.

---

## 3. Parameters

| Symbol | Code | Value | Meaning |
|---|---|---|---|
| `α^C` | `alpha_C` | 0.85 | Commercial share of the non-floor-space input in production |
| `α^R` | `alpha_R` | 0.66 | Residential expenditure share on the non-housing good |
| `β^C` | `beta_C` | 0.03 | Commercial **agglomeration** elasticity |
| `β^R` | `beta_R` | 0.00 | Residential agglomeration elasticity (off by default) |
| `τ^C` | `tau_C` | 0.01 | Commercial productivity decay, per km |
| `τ^R` | `tau_R` | 0.005 | Residential amenity decay, per km |
| `ω^C` | `omega_C` | 0.03 | Commercial rent elasticity of height, `ω^U = ω̃^U/(1-α^U)` |
| `ω^R` | `omega_R` | 0.07 | Residential rent elasticity of height |
| `θ^C` | `theta_C` | 0.5 | Commercial **construction-cost elasticity of height** |
| `θ^R` | `theta_R` | 0.55 | Residential construction-cost elasticity of height |
| `c^U` | `c_C`, `c_R` | 1.4 | Baseline construction cost |
| `ā^C` | `a_bar_C` | 2 | Fundamental commercial productivity |
| `ā^R` | `a_bar_R` | 1 | Fundamental residential amenity |
| `r_a` | `r_a` | 150 | Agricultural land rent (a per-observation *variable* in Stata) |
| `S̄^U` | `S_bar_C`, `S_bar_R` | 999 | Height limits (999 ≈ unconstrained) |
| — | `y` | 2.5 (start) | City-wide wage — **target object** |
| — | `L` | 1e6 (start) | Total employment `N` — **target object** |

The toolkit uses **three different parameterisations** and they are not interchangeable:

| Context | `c_C`, `c_R` | `r_a` | Source |
|---|---|---|---|
| `AB2022.ado` defaults (the `Params` defaults here) | 1.4 | 150 | `AB2022.ado:387-389` |
| `_2_ANALYSIS.do` baseline figures (`params_analysis()`) | 1.3 | 50 | `_2_ANALYSIS.do:23` |
| `_1_PROGS.do` CONV/EMP inversion (`params_conv()`) | 1.4 | 30 | `_1_PROGS.do:276,288` |

### 3.1 Notation-collision warning

**Read this before copying any calibrated value out of this folder.** These letters mean different
things here than elsewhere in this repository; a transplanted value will look plausible and be
silently wrong.

| Symbol | **AB2022 (here)** | Elsewhere in `Quantitative-Spatial-Economics/` |
|---|---|---|
| `θ` | construction-cost elasticity of height, **0.5–0.55** | Fréchet/trade elasticity, **6–10** |
| `β` | **agglomeration** elasticity, **+0.03** | **amenity externality**, **−0.3** — opposite sign *and* opposite economics |
| `α` | non-floor-space input/expenditure share, **0.66–0.85** | productivity externality, **0.1** |
| `ω` | rent elasticity of height | not used |
| `λ` | not used | traffic-congestion elasticity 0.07–0.09 — and, in `Ahlfeldt/ARSW2015`, an *agglomeration* elasticity at nearly the same magnitude |
| `σ`, `ε` | not used | CES substitution ≈4; migration elasticity 3–5 |

---

## 4. Solution Algorithm

Nested fixed-point iteration, mirroring the Codebook's Algorithms 1–6 one-for-one. No optimiser, no
analytical Jacobian.

| Codebook | Julia | What it does |
|---|---|---|
| Alg. 1 `SOLVER` | `solver!` | One deterministic recursive pass given `(y, N)`. No iteration inside. |
| Alg. 2 `WAGE` | `wage!` | Inner loop: while `\|L̂/N̂ − 1\| > 0.01`, set `y ← 0.5y + 0.5·y·(L̂/N̂)^0.01`, re-solve. |
| Alg. 3 `FINDEQ` | `findeq!` | Outer loop on `N`: solve, clear the wage, then `L ← 0.5L + 0.25(L̂ + N̂)`. |
| Alg. 4 `CONV` | `conv!` | Nudge `ã^U` toward `HEIGHT^U/S*^U`, re-solve. |
| Alg. 5 `EMP` | `emp!` | Scale `ã^R` by `(N̄/N)^0.01` toward a target population. |
| Alg. 6 `INVERT` | `invert!` | Drive `conv!` to `R² ≥ 0.999`, then `emp!` to `\|N − N̄\| ≤ 1000`. |

Both loops use the upstream **1% relative tolerance** with fixed damping. Starting values `y = 2.5`,
`N = 10^6`.

**Off-by-one in `findeq!`, faithfully preserved.** The objective is recomputed *before* `L` is
updated, and the `L` update is the last statement in the loop body. So on exit the stored `L` is one
update ahead of the endogenous columns it is paired with. This is not cosmetic — reproducing the
shipped `INVERTED.dta` requires recovering
`L_prev = 2·(L_stored − 0.25·(L̂ + N̂))`, without which `A_tilde_x_C` is off by
`β^C·ln(L_stored/L_prev) = 1.4e-04` instead of the `8.5e-08` noise floor.

---

## 5. The ω Discrepancy — and Which Convention the Shipped Data Actually Used

The two Stata sources implement **different formulas** for the horizontal bid rent (eqs. 9–10):

```
stata_source/AB2022.ado:97-98      p_bar_x_C = a_x_C * 1/(1-omega_C) * S_x_C^omega_C
stata_source/_1_PROGS.do:132-133   p_bar_x_C = a_x_C * 1/(1+omega_C) * S_x_C^omega_C
```

`AB2022-codebook.pdf` eqs. 9–10 read `a(x)/(1+ω^U)` — the Codebook agrees with `_1_PROGS.do`,
against the `.ado`.

**`(1+ω)` is theoretically correct**, and this is settled without reference to any data: it is the
only form whose FOC generates the `S*` that both files use, and the only one under which `p̄·S` is
the first term of the land rent (§2). The `.ado` is internally inconsistent.

Both are implemented, selected by `bid_rent_convention`. **The default is `:codebook`.**

### The decisive evidence

`reference_data/BASE.csv` **cannot settle the question** — it holds no solved values at all (§6.1).
`reference_data/INVERTED.csv` can, and does:

| Convention | Worst relative deviation across all 19 solved columns |
|---|---|
| `:ado` — `1/(1-ω)` | **2.49e-06** ✓ |
| `:codebook` — `1/(1+ω)` | 1.51e-01 ✗ |

**The shipped `INVERTED.dta` was generated with `(1−ω)`, the `.ado` convention** — despite the
Codebook, the walkthrough and the developer's own FOC all saying `(1+ω)`.

The reference data is in fact a **hybrid that neither current source file reproduces on its own**:

| Block | `INVERTED.csv` matches | Does *not* match |
|---|---|---|
| Land rent (eqs. 5–6) | `_1_PROGS.do:104-105` / Codebook | `AB2022.ado:69-70` (out by a factor ~4) |
| Bid rent (eqs. 9–10) | `AB2022.ado:97-98` | `_1_PROGS.do:132-133` / Codebook (out by 5.8%/13.1%) |

Most plausible history: the original solver had the `c(1+θ)` land-rent bracket *and* the `(1−ω)` bid
rent; the `.ado` (v0.93, 02/2024) later broke the land-rent bracket while `_1_PROGS.do` (01/2024)
fixed the bid rent to match the Codebook — and `INVERTED.dta` predates both edits.

### Size of the effect

Pointwise the ratio is exactly `(1+ω)/(1−ω)`: **1.0619** commercial, **1.1505** residential. In
general equilibrium (`_2_ANALYSIS.do` parameters), switching `:codebook → :ado`:

| | `:codebook` (1+ω) | `:ado` (1−ω) | change |
|---|---|---|---|
| max floor rent, commercial | 16.4721 | 17.8444 | +8.33% |
| max floor rent, residential | 13.5268 | 15.6048 | +15.36% |
| wage `y` | 2.3294 | 2.3342 | +0.20% |
| employment `N` | 1,726,309 | 2,009,515 | **+16.41%** |
| CBD radius `x0` | 4.47 km | 4.72 km | +5.59% |
| urban radius `x1` | 27.75 km | 28.16 km | +1.48% |

A 3-point bid-rent convention difference moves equilibrium city population by 16%.

---

## 6. Validation

Everything below is produced by running the file. Five tiers, all passing.

### 6.0 The precision floor

The reference CSVs carry ~8 significant decimal digits, printed from Stata's 4-byte floats — **and
the wage `y` is stored to only six** (`2.10701`). Since `a^C ∝ y^{-5.667}`, which is then raised to
`1/(θ^C-ω^C) = 2.128` for height and compounded again for labour demand, half an ulp of that printed
`y` is amplified by up to ~21×.

So the pass criterion is not an arbitrary constant. For each column the tolerance is **the band the
reference's own stored precision permits**: re-solve at `y ± half-ulp` and take the induced spread,
floored at `4·eps(Float32)/2 ≈ 2.4e-07` for columns with no `y`-sensitivity. A `dev/tol` ratio well
below 1 means the port agrees more closely than the data can resolve.

### 6.1 Tier A — `BASE.csv` (10,001 × 28): **PASS**

**Finding: `BASE.csv` is not a solved equilibrium.** `_1_PROGS.do:81` saves `BASE.dta` *between*
generating the placeholder variables (`:48-78`) and defining `SOLVER` (`:90`). 21 of its 28 columns
are entirely missing — 0 non-null out of 10,001 each.

| Columns | n_ref | max abs dev | max rel dev | Status |
|---|---|---|---|---|
| `x`, `D`, `y`, `L`, `r_a`, `a_rand_C`, `a_rand_R` | 10,001 | **0.000e+00** | **0.000e+00** | exact |
| `a_x_C`, `a_x_R`, `A_tilde_x_C`, `A_tilde_x_R`, `r_x_C`, `r_x_R`, `U`, `S_star_x_C`, `S_star_x_R`, `S_x_C`, `S_x_R`, `S_x`, `p_bar_x_C`, `p_bar_x_R`, `L_x_C`, `f_bar_x_R`, `n_x`, `URBAN`, `COM`, `SHADE`, `SHADEU` | 0 | — | — | all-missing both sides |

The 7 populated columns reproduce **exactly** (the grid `x = (n-5001)/100` is exactly representable).
The all-missing check is a real test — it catches a port that silently initialises to `0.0` instead
of `NaN` — but not an economic one. Note `r_a = 100` here, the `gen` placeholder from
`_1_PROGS.do:53`, which is none of the three solver parameterisations.

**Consequence:** `BASE.csv` cannot discriminate the ω convention and cannot supply land-use
boundaries. Both requirements move to `INVERTED.csv`.

### 6.2 Tier B — `INVERTED.csv` (4,001 × 49): **PASS**

The only solved gold standard in the folder. Taking its own exogenous state as given
(`a_rand_C`, `a_rand_R`, `y`, `r_a = 30`, and the `CONV`/`EMP` parameters), one `solver!` pass
reproduces every solved column. `L_stored = 1004690.56`, `L_prev = 999993.537290` (ratio 1.00469705).

| column | n_ref | max abs dev | max rel dev | tolerance | dev/tol | status |
|---|---|---|---|---|---|---|
| `A_tilde_x_C` | 4001 | 2.448e-07 | 8.513e-08 | 2.384e-07 | 0.36 | OK |
| `A_tilde_x_R` | 4001 | 1.426e-07 | 1.400e-07 | 2.384e-07 | 0.59 | OK |
| `a_x_C` | 4001 | 1.507e-05 | 7.458e-07 | 1.345e-05 | 0.06 | OK |
| `a_x_R` | 4001 | 6.112e-06 | 5.144e-07 | 6.980e-06 | 0.07 | OK |
| `r_x_C` | 4001 | 3.014e-03 | 2.402e-06 | 4.292e-05 | 0.06 | OK |
| `r_x_R` | 4001 | 5.407e-04 | 1.550e-06 | 2.254e-05 | 0.07 | OK |
| `U` | 4001 | 0.000e+00 | 0.000e+00 | 2.384e-07 | 0.00 | OK |
| `S_star_x_C` | 242 | 2.350e-04 | 1.638e-06 | 2.861e-05 | 0.06 | OK |
| `S_star_x_R` | 1866 | 4.933e-05 | 9.938e-07 | 1.454e-05 | 0.07 | OK |
| `S_x_C` | 242 | 2.350e-04 | 1.638e-06 | 2.861e-05 | 0.06 | OK |
| `S_x_R` | 1866 | 4.933e-05 | 9.938e-07 | 1.454e-05 | 0.07 | OK |
| `S_x` | 2108 | 2.350e-04 | 1.638e-06 | 2.861e-05 | 0.06 | OK |
| `p_bar_x_C` | 242 | 1.837e-05 | 8.244e-07 | 1.431e-05 | 0.06 | OK |
| `p_bar_x_R` | 1866 | 1.015e-05 | 5.848e-07 | 7.997e-06 | 0.07 | OK |
| `L_x_C` | 242 | 2.903e-02 | 2.491e-06 | 4.529e-05 | 0.06 | OK |
| `f_bar_x_R` | 1866 | 3.185e-08 | 5.263e-07 | 5.624e-06 | 0.09 | OK |
| `n_x` | 1866 | 2.708e-03 | 1.505e-06 | 2.016e-05 | 0.07 | OK |
| `URBAN` | 4001 | 0.000e+00 | 0.000e+00 | 2.384e-07 | 0.00 | OK |
| `COM` | 4001 | 0.000e+00 | 0.000e+00 | 2.384e-07 | 0.00 | OK |
| `x`, `D`, `r_a`, `a_rand_C`, `a_rand_R` | 4001 | 0.000e+00 | 0.000e+00 | — | — | exogenous, exact |

**Worst relative deviation 2.491e-06, worst dev/tolerance ratio 0.59.** The error chain is pure
input-quantization amplification: `a_rand` (float32) → `A_tilde` (8.5e-08) → `^6.667` → `a_x_C`
(7.5e-07) → `^2.128` → `S*` (1.6e-06) → `L_x_C` (2.5e-06). Every column agrees at 6–59% of what the
data's own precision permits.

Two replication details were mandatory, both non-obvious:

- **(a)** the off-by-one `L` (§4);
- **(b)** `S_star_x_C`/`S_star_x_R` are **never cleared** between solver passes — the Stata clear-list
  (`AB2022.ado:58`, `_1_PROGS.do:93`) lists 14 variables and omits both. They retain stale values at
  locations whose land use has since changed: in `INVERTED.csv`, `S_star_x_C` is non-missing at 1,387
  rows where `U ≠ 1`, and `S_star_x_R` at 2,135 rows where `U ≠ 2`. The port replicates this and
  compares `S_star` only within its own zone.

### 6.3 Tier C — general equilibrium: **PASS**

| | value | test |
|---|---|---|
| Labour demand `L̂` | 1,007,078.26 | |
| Labour supply `N̂` | 1,011,696.91 | |
| `\|L̂/N̂ − 1\|` | **0.00457** | < 0.01 ✓ |
| `\|L/(0.5(L̂+N̂)) − 1\|` | **0.00465** | < 0.01 ✓ |

`FINDEQ` stops at a 1% relative tolerance with fixed damping, so its terminal `(y, N)` is **not** a
precise fixed point — it is path-dependent, and a re-solve from a different start legitimately lands
elsewhere inside the 1% band. Exact reproduction of the *iteration path* that produced
`INVERTED.dta` is therefore not attainable. That is a property of the source, not of this port.

### 6.4 Tier D — land-use boundaries: **PASS**

`U` identical at **4001 / 4001** rows. `x0 = 0.08` km, `x1 = 4.46` km, matching at grid resolution.

`x1` is computed two ways upstream — `.ado:79-80` takes `max(x | U<3)`, `_1_PROGS:114-115` takes
`min(x | U==3)`, differing by exactly one grid step. The port uses the `_1_PROGS` form. Both are
reported diagnostics; neither enters any equation.

### 6.5 Tier E — Chicago skyline inversion: **PASS**

Merging `EMPIRICAL_CH_skyline.csv` (558 rows) and binning at 0.1 km gives `HEIGHT_C` at 242 rows and
`HEIGHT_R` at 1,883 — matching the reference's 242/1,866 in its `|x| ≤ 20` window.

| | result | target |
|---|---|---|
| height correlation `R²(HEIGHT_R, S_x_R)` | **0.99957** after 5 `CONV` steps | ≥ 0.999 ✓ |
| population gap | **641.2** after 205 `EMP` steps | ≤ 1000 ✓ |

Recovered equilibrium vs. the reference:

| | port | reference | gap |
|---|---|---|---|
| wage `y` | 2.1131 | 2.10701 | 0.29% |
| employment `N` | 1,000,641 | 999,994 | 0.06% |
| `x0` | 0.08 km | 0.08 km | exact |
| `x1` | 4.45 km | 4.46 km | one grid step |
| `ã^C` zero share | 94.0% | 94.0% | — |
| `ã^R` zero share | 52.9% | 53.4% | 0.5pp |
| `ã^C` mean (>0) | 0.9487 | 0.9428 | 0.6% |
| `ã^R` mean (>0) | 1.0913 | 1.0822 | 0.8% |

The recovered amenities are compared **in distribution, not pointwise**: the inversion path is as
tolerance-limited as the equilibrium (§6.3), so an exact pointwise match is neither expected nor
claimed. Landing within 0.3% on the wage and 0.06% on population after a 205-step inversion through
a 1%-tolerance solver is a strong result, not a weak one.

---

## 7. Counterfactuals

### Height limit (`_2_ANALYSIS.do:38`): 20 floors on both uses

| | baseline | height limit | change |
|---|---|---|---|
| max commercial height | 75.70 | 20.00 | −73.58% |
| max residential height | 36.05 | 20.00 | −44.52% |
| max floor rent (C) | 16.4721 | 16.0590 | −2.51% |
| max land rent (C) | 390.71 | 409.26 | +4.75% |
| wage `y` | 2.3294 | 2.2653 | −2.75% |
| employment `N` | 1,726,309 | 841,490 | **−51.25%** |
| CBD radius `x0` | 4.47 km | 6.36 km | +42.28% |
| urban radius `x1` | 27.75 km | 22.17 km | −20.11% |

A height limit does far more than cap the skyline: it halves city population and shrinks the urban
footprint by 20%, because constrained developers cannot bid as much for land. The **CBD spreads**
(`x0` +42%) as firms substitute from the vertical to the horizontal margin, while the city overall
contracts.

> **Caveat — read §8, defect 3.** The toolkit evaluates land rent at the *unconstrained* height, so
> the `+4.75%` land-rent figure and the boundary shifts in this table inherit that error. Under
> Codebook eq. 5 the same cap gives max land rent **206.79** rather than 409.26, `N` = 834,255,
> `x0` = 6.23 km, `x1` = 21.90 km. The file reports both.

### Commercial subcenter (`_4_INVERTEDCOUNTER.do`)

From the Chicago-inverted fundamentals, inject `ã^C ~ 1.15·U(0,1)` between 12 and 15 km, spread it
to 0.1 km bins, and re-solve:

| | inverted | with subcenter | change |
|---|---|---|---|
| max commercial height, `\|x\| ≤ 5` km | 189.65 | 150.75 | **−20.51%** |
| max commercial height, 11–16 km | 0.00 | 156.20 | new subcenter |
| wage `y` | 2.1131 | 2.2243 | **+5.26%** |
| employment `N` | 1,000,641 | 1,693,621 | +69.25% |

Reallocation *and* growth: the original CBD loses a fifth of its height to the new subcenter, while
the city as a whole expands because it has become more productive overall. The wage effect matches
`_4_INVERTEDCOUNTER.do`'s own annotation ("Wage increases by about 5%") almost exactly; the
population effect is larger than its "+30%" because **upstream calls `runiform()` with no seed**, so
the shock differs on every run and levels are not comparable across runs. This port seeds it
(`seed=2022`) so it at least reproduces itself.

Note `x1` is unchanged here: under the `_1_PROGS` definition (`min(x | U==3, x ≥ 0)`) it marks the
*first* agricultural location, so with a spatially disconnected subcenter it stops at the gap between
the CBD and the subcenter rather than at the outer edge of the urban area.

---

## 8. Upstream Defects

| # | Defect | Handling |
|---|---|---|
| 1 | **Bid-rent ω:** `.ado:97-98` uses `1/(1-ω)`, `_1_PROGS.do:132-133` and Codebook eqs. 9–10 use `1/(1+ω)`. | Both implemented, default `:codebook`. §5. |
| 2 | **Land-rent bracket:** `.ado:69-70` writes `c(1+ω)` where the Codebook and `_1_PROGS.do:104-105` write `c(1+θ)` — inconsistent with the `.ado`'s own `S*`, and out by a factor ~4 against the shipped data. | Treated as a plain bug; only the Codebook form implemented. |
| 3 | **Land rent ignores the height limit** — found by this port, not in the toolkit's notes. Both Stata files inline the closed form, so `S̄` never appears in eqs. 5–6, and land rent is evaluated at the *unconstrained* `S*` rather than Codebook eq. 5's `S̃ = min(S̄, S*)`. Irrelevant when the cap does not bind (so Tiers A–E are unaffected), but at the centre with `S̄^C = 20` the toolkit gives **409.26** where eq. 5 gives **204.90** — a **99.7% overstatement**. Since land rent drives the land-use allocation, the height-limit counterfactual's boundaries inherit it. | Both implemented via `constrained_land_rent`; **default `true` = Codebook eq. 5, i.e. corrected**. Pass `false` to reproduce the toolkit. The counterfactual reports both. |
| 4 | **`S_star` never cleared** between solver passes (`.ado:58`, `_1_PROGS.do:93`). | Replicated deliberately — the reference data depends on it. §6.2(b). |
| 5 | **No `maxiter` guard** on either solver loop. A parameterisation that never reaches the 1% tolerance without tripping the `L̂+N̂ == 0` check hangs Stata forever. | `maxiter` + `@warn` added to `wage!`, `findeq!`, `invert!`. |
| 6 | **Undeclared SSC dependency `grc1leg`**, despite the toolkit README claiming no user-written ado files. | Stata-figure-only; irrelevant to this port. Noted for completeness. |
| 7 | `.ado:238` defines a `FINDEQ` that is **never called** — `AB2022` re-implements the same loop inline at `:414-446`. Three copies exist upstream. | One here. |
| 8 | `_1_PROGS.do:110-111` inline comments for the `U==2`/`U==1` branches are **swapped**. | Comment-only; logic matches. |
| 9 | `_1_PROGS.do:286-289` — `EMP`'s `foreach name in R C` loop never uses `name`, so the body runs **twice** per call. | Replicated; it changes the scaling. |
| 10 | `Ahlfeldt/AB2022-toolkit/CLAUDE.md` states `CONV`'s λ = 0.5. It is **0.05** (`_3_INVERSION.do:91` → `_1_PROGS.do:304`). | Corrected here. `Ahlfeldt/` is read-only and was not modified. |

---

## 9. Relation to the Other Models in This Repository

**This is the only model here with an endogenous building-height margin.** In
`AllenArkolakis-RES-2022`, `FuchsFoongWong-MMN-2026`, `Santamaria-2022`, `QRE-HoRaUE-2025`,
`Redding-JIE-2016`, `QSE-ARE-2017` and `BrockhausHinzSerfaty-2026`, floor/land supply is fixed, an
exogenous elasticity, or absent; every location is a point with only a horizontal trade/migration
problem. The developer's height problem and the construction-cost elasticity `θ^U` are entirely
absent from those gravity/trade/traffic models.

The complementarity runs the other way too: **transport cost is exogenous here — indeed absent**
(§1) — which is precisely the margin those models endogenize.

| Mechanism | AB2022 (here) | Rest of this repository |
|---|---|---|
| Endogenous building height / FAR | **yes** | absent |
| Endogenous transport cost (route choice, congestion) | absent (no commuting term at all) | AA2022 (Leontief inverse, λ), FFW2026 (multimodal CES), Santamaria-2022 (Dijkstra + investment), Optimal Transport Networks |
| Multi-region CES goods trade | absent (one city) | all of them |
| Bilateral commuting | absent | absent — see `Ahlfeldt/ARSW2015`, `MRRH2018` |

Within Ahlfeldt's own ladder (`Ahlfeldt/CLAUDE.md`):

- **`ABSJ2027`** (Ahlfeldt, Baum-Snow & Jedwab, *ReStud* forthcoming) is AB2022's **direct
  successor**: same vertical margin and 1-D geometry, but *imperfectly* open (migration closure) and
  with an explicit urban growth boundary alongside height caps, calibrated across 12,873 cities.
- **`ARSW2015`** (Berlin blocks) and **`MRRH2018`** (German counties) **drop the vertical margin**:
  floor space there is *measured data*, not a developer's choice. They add bilateral commuting
  instead.

The obvious synthesis target: nest this bid-rent/height block inside each location of a
`Redding-JIE-2016`-style multi-region model, replacing AB2022's single-region labour clearing with a
gravity migration closure, so that local floor-space supply is a developer's problem rather than an
exogenous elasticity.

---

## 10. Contents of This Folder

```
AhlfeldtBarr-JUE-2022/
├── ahlfeldt_barr_skyscraper_model.jl   THE PORT -- single self-contained file
├── TRANSLATION_PLAN.md                 equation mapping, tolerances, decisions, outcome
├── README.md                           this file
├── graphs/                             12 PDFs: {height, floorrent, landrent}
│                                       x {baseline, heightlimit, chicago, subcenter}
├── AB2022-paper.pdf                    JUE 129, 103419 (21 pp)  + MinerU markdown
├── AB2022-appendix.pdf                 Online Appendix          + MinerU markdown
├── AB2022-codebook.pdf                 the authoritative economic documentation (Table 1, §A.1, §A.2)
├── LICENSE-AB2022-toolkit              MIT
├── stata_source/
│   ├── AB2022.ado                      the shipped self-contained solver (470 lines)
│   ├── AB2022.sthlp                    Stata help file
│   └── _0_META.do … _4_INVERTEDCOUNTER.do   the five walkthrough scripts
└── reference_data/
    ├── BASE.csv                        10,001 x 28 -- the UNSOLVED placeholder grid (§6.1)
    ├── INVERTED.csv                    4,001 x 49  -- the solved, Chicago-inverted model
    └── EMPIRICAL_CH_skyline.csv        558 rows -- observed Chicago skyline, the inversion target
```

Dependencies: `Statistics`, `Printf`, `Plots`, `DataFrames`, `CSV`. No `Project.toml` — this
repository deliberately has none; install packages globally. Julia 1.12.4.

The upstream checkout lives at `Ahlfeldt/AB2022-toolkit/` and is **read-only** (it is a separate git
repo with the author's remote, and is gitignored by the parent). Nothing under it was modified.

---

## 11. Key References

- **Ahlfeldt, G. M. & J. Barr (2022)**, "The economics of skyscrapers: A synthesis", *Journal of
  Urban Economics* 129, 103419.
- **Ahlfeldt, G. M. & J. Barr (2024)**, *Codebook for: The economics of skyscrapers*, v0.95 —
  `AB2022-codebook.pdf`, the authoritative statement of the 13 equations and Algorithms 1–6.
- **Ahlfeldt, G. M., N. Baum-Snow & R. Jedwab**, "The Skyscraper Revolution", *Review of Economic
  Studies*, forthcoming — the successor model (`Ahlfeldt/ABSJ2027-toolkit`).
- **Ahlfeldt, G. M., S. J. Redding, D. M. Sturm & N. Wolf (2015)**, "The Economics of Density:
  Evidence from the Berlin Wall", *Econometrica* 83(6) — drops the vertical margin, adds commuting.
- **Alonso, W. (1964)**; **Mills, E. S. (1967)**; **Muth, R. F. (1969)** — the monocentric-city
  foundation onto which the vertical margin is grafted.

Cite Ahlfeldt & Barr (2022) when using this code. The toolkit is MIT-licensed; see
`LICENSE-AB2022-toolkit`.
