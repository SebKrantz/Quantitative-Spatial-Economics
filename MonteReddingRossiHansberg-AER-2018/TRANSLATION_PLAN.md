# TRANSLATION_PLAN.md — MRRH (2018) MATLAB toolkit → Julia

Port target: `monte_redding_rossihansberg_commuting_model.jl`, single self-contained file, Julia 1.12.4.

Source: Ahlfeldt & Seidel, *Toolkit for Quantitative Spatial Models* (v0.91/0.92, MIT), implementing
Monte, Redding & Rossi-Hansberg (2018, *AER* 108(12), 3855-90) in the Seidel & Wickerath (2020,
*RSUE* 85) German-county variant. Upstream checkout: `Ahlfeldt/MRRH2018-toolkit/` (read-only).

Status legend used below: **[planned]** written before coding, **[DONE]**/**[CHANGED]**/**[NOT DONE]**
appended in §9 after implementation.

---

## 1. Model in one page

Locations `n, i ∈ 1..J`. `n` = residence/consumption, `i` = workplace/production.
Workers draw Fréchet(ε) idiosyncratic amenities over every **(residence, workplace) pair** — this
two-margin discrete choice is the model's defining feature and is absent from every other model in
this repository.

Two nested layers:

* **MRRH2018 core**: exogenous productivity `A_i`, exogenous land `H_n`, inelastic housing.
* **SW2020 extension** (what the toolkit actually runs): endogenous productivity `A_i = Ā_i L_i^ν`
  and endogenous housing supply `H_n = H̄_n P_{H,n}^δ`. `ν = δ = 0` collapses to MRRH2018.

Equilibrium is 6 vectors + 1 scalar `{w_n, v̄_n, Q_n, L_n, R_n, P_n; Ū}` solving 7 equations
(Codebook §A.1), plus 2 more equations/variables in the SW2020 variant.

---

## 2. Equation-by-equation mapping

Paper = MRRH2018 AER equation numbers; Codebook = `MRRH2018-codebook.pdf` §A.1/A.2/A.3.
MATLAB paths are relative to `matlab_source/`.

### 2.1 Quantification (inversion)

| # | Object | Paper eq | Codebook | MATLAB file:line | Planned Julia function |
|---|---|---|---|---|---|
| Q1 | Conditional commuting prob. `λ_{ni|n} ∝ B_i (w_i/κ_ni)^ε` | (12) | Alg 1 step 2 | `getBiTK.m:67-68` | `getBi` (inner) |
| Q2 | Predicted workplace employment `L_i = Σ_n λ_{ni|n} R_n` | (13) | Alg 1 step 3 | `getBiTK.m:69` | `getBi` (inner) |
| Q3 | Workplace-amenity update `B_i ← B_i · L_i^obs/L_i^pred` | — | Alg 1 steps 4-5 | `getBiTK.m:70` | `getBi`, with `damp` kwarg |
| Q4 | Residence choice prob. `λ^R_n = R_n/ΣR` | (11) | Alg 1 step 6 | `getBiTK.m:79` | `getBi` |
| Q5 | Unconditional `λ_ni = λ_{ni|n} λ^R_n`; flows `L_ni = λ_ni L̄` | (10),(11) | Alg 1 steps 7-8 | `getBiTK.m:82,84` | `getBi` |
| Q6 | Residential income `v̄_n = Σ_i λ_{ni|n} w_i` | (14) | — | `OwnData.m:82` | `residential_income` |
| Q7 | Trade share `π_ni ∝ Ā_i^{σ-1} L_i^{1-(1-σ)ν} w_i^{1-σ} d_ni^{1-σ}` | (6) | Alg 2 step 2 | `solveProductTradeTK.m:65-73` | `trade_shares` |
| Q8 | Income = expenditure `w_i L_i = Σ_n π_ni v̄_n R_n` | (7),(16) | Alg 2 steps 3-4 | `solveProductTradeTK.m:77-84` | `solve_productivity` |
| Q9 | Productivity update `Ā ← ζ Ā(inc/exp) + (1-ζ)Ā`, renormalized | — | Alg 2 steps 5-6 | `solveProductTradeTK.m:112-114` | `solve_productivity` |
| Q10 | Price index `P_n = σ/(σ-1) (L_n/(σFπ_nn))^{1/(1-σ)} w_n/A_n` | (8) | Alg 2 step 7 | `solveProductTradeTK.m:107` | `price_index` |
| Q11 | Land market clearing `Q_n = (1-α) v̄_n R_n / H_n` | (5) | A.1 eq 3 | *(never computed upstream)* | `land_price` — **added** |
| Q12 | Full amenity `B_ni = λ_ni Φ κ_ni^ε (P_n^α Q_n^{1-α})^ε w_i^{-ε}` | (10),(17) | — | *(never computed upstream)* | `invert_amenities` — **added** |

Q11/Q12 do not exist upstream. They are required for the levels round-trip test (§5) and for the
forward equilibrium solver (§4.2). See decision D6.

### 2.2 Counterfactuals (exact hat algebra)

All eight updaters are one-liners; each maps to a Codebook algorithm.

| # | Object | Paper eq | Codebook | MATLAB | Planned Julia |
|---|---|---|---|---|---|
| C1 | `v̂_n` residential income | (14) / App B.17 | Alg 3 | `updateResWageTK.m` | `update_res_income` |
| C2 | `L̂_i = L̄ Σ_n λ_ni λ̂_ni / L_i` | (11) / App B.23 | Alg 4 | `updateEmplTK.m` | `update_employment` |
| C3 | `R̂_n = L̄ Σ_i λ_ni λ̂_ni / R_n` | (11) / App B.22 | Alg 5 | `updateResidentsTK.m` | `update_residents` |
| C4 | `P̂_{H,n} = (v̂_n R̂_n)^{1/(1+δ)}` | SW2020 (9) | Alg 6 | `updateHousePriceTK.m` | `update_house_price` |
| C5 | `π̂_ni` | (6) / App B.20 | Alg 7 | `updateTradeshTK.m` | `update_trade_shares` |
| C6 | `P̂_{Q,n} = (L̂^{1-(1-σ)ν}/π̂_nn)^{1/(1-σ)} d̂_nn ŵ_n/Â_n` | (8) | Alg 8 | `updatePricesTK.m` | `update_price_index` |
| C7 | `w̃_i = Σ_n π_ni π̂_ni v̄_n v̂_n R_n R̂_n / (w_i L_i L̂_i)` | (7) / App B.16 | Alg 9 | `updateWageTK.m` | `update_wage` |
| C8 | `λ̃_ni` | (10) / App B.20 | Alg 10 | `updateLamTK.m` | `update_lambda` |
| C9 | Outer fixed point on `{ŵ, λ̂}`, damping ζ=0.25 | — | Alg 11 | `counterFactsTK.m:87-123` | `counterfactual` |
| C10 | Welfare `Û = B̂^{1/ε} (κ̂ P̂_Q^α P̂_H^{1-α})^{-1} ŵ λ̂^{-1/ε}` | (15),(24) | A.2 | `counterFactsTK.m:128-131` | `counterfactual` |

### 2.3 Descriptives

| # | Object | MATLAB | Planned Julia |
|---|---|---|---|
| D1 | Commuting gravity, origin+destination FE | `Descriptives.m:90-101` | `gravity_fe` (alternating projections) |
| D2 | Trade gravity, same | *(not upstream)* | `gravity_fe` — **added** |
| D3 | CMA `CMA_n = Σ_i dist_ni^{-με} w_i^ε` | `Descriptives.m:115-117` | `commuting_market_access` |
| D4 | Employment potential `Σ_i dist^{-με} L_i` | `Descriptives.m:119` | `commuting_market_access` |

---

## 3. Parameters

### 3.1 Toolkit / SW2020 German county track (the port's default)

| Symbol | Code name | Value | Meaning | Source |
|---|---|---|---|---|
| α | `alpha` | 0.70 | expenditure share on **tradable goods** (1-α = 0.30 on housing) | `MRRH2018_toolkit.m:49` |
| ε | `epsilon` | 4.60 | Fréchet shape over **(residence, workplace) pairs** | `:50` |
| μ | `mu` | 0.47 | elasticity of commuting cost κ to travel time/distance | `:51` |
| δ | `delta` | 0.38 | housing supply elasticity | `:52` |
| σ | `sigma` | 4.0 | CES elasticity of substitution across varieties | `:53` |
| F | `fixC` | 1.0 | fixed cost of production (labour units) | `:54` |
| ν | `nu` | 0.05 | agglomeration elasticity, `A_i = Ā_i L_i^ν` | `:55` |
| ψ | `psi` | 0.42 | distance elasticity of iceberg trade cost `d_ni = dist^ψ` | `:56` |
| J | `J` | 401 | German counties (Kreise, 2018 definition) | `:57` |

Implied composites: commuting gravity slope `-εμ = -2.162`; trade gravity slope `ψ(1-σ) = -1.26`.

### 3.2 MRRH2018's own US calibration (for the comparison run)

| Symbol | Value | Source |
|---|---|---|
| α | 0.60 (1-α = 0.40 housing) | paper §II, BEA |
| ε | 3.30 | paper §II, second-stage gravity with wage instrument |
| φ ≡ με | 4.43 | paper §II, commuting gravity with workplace+residence FE ⟹ μ = 4.43/3.30 = 1.3424 |
| σ | 4 | Broda & Weinstein (2006) |
| ψ | 0.43 (from -ψ(σ-1) = -1.29) | paper §II, CFS trade gravity |
| δ | 0 (inelastic land, baseline) | paper §I |
| ν | 0 (no agglomeration in MRRH2018) | paper §I |

### 3.3 GRID track (documented, not the default) — defect 2

`GRID_MRRH2018_toolkit.m:47-55` sets **ν = 0.0**, not 0.05. All other parameters coincide.
The port takes `nu` as an explicit field of the parameter struct with default 0.05 (county track) and
ships a named `GRID_PARAMS`-style constructor with `nu = 0.0`, so neither value is inherited silently.

### 3.4 Notation collisions — MUST be documented in the code header

Same letters, different economics, across this repository (`Ahlfeldt/CLAUDE.md` §"Notation collisions"):

| Symbol | MRRH2018 (this file) | Elsewhere in this repo |
|---|---|---|
| ε | Fréchet shape over commuting pairs, **4.6** | migration elasticity 3-5 (AA2022, Redding 2016) |
| α | consumption share on **tradables**, **0.70** | productivity externality **0.1** (AA2022) |
| β | *(unused here)* | amenity externality **-0.3** (AA2022) |
| ν | agglomeration elasticity, **0.05** | *(unused elsewhere)* |
| λ | commuting **probability** `λ_ni` | traffic **congestion elasticity** 0.07-0.09 (AA2022) |
| θ | *(unused here)* | Fréchet/trade elasticity 6-10 |
| μ | travel-time elasticity of commuting cost, **0.47** | — |
| π_ni | **trade** share (destination n, origin i) | commuting probability in ARSW2015 |
| σ | CES across varieties, 4 | same meaning |

The dangerous one is **λ**: a probability here, a congestion elasticity in `AllenArkolakis-RES-2022`.
The second most dangerous is **α**: 0.70 here is a *goods* share; 0.1 there is a *productivity spillover*.

---

## 4. Solver architecture

### 4.1 Quantification pipeline (primary entry point = OwnData track)

```
load_county_data()                          # 401 counties, CSVs in reference_data/
   ├─ dist (J×J, metres → km)
   ├─ L_obs, R_obs (from commuting_wide.csv column/row sums), each ÷ mean
   ├─ w_obs (labor_tidy.csv median_income_workplace) ÷ mean
   ├─ area, rent index, border distance
   ↓
build_costs(dist; psi, mu)                  # d_ni = (dist/min dist)^psi ; κ_ni = dist
   ↓
getBi(w, kappa, R_obs, L_obs, Lbar)         # Alg 1 — OUTER loop only, no nesting
   → B_i, λ_{ni|n}, λ_ni, L_pred, flows
   ↓
v̄ = λ_{ni|n} * w                            # Eq (14)
   ↓
solve_productivity(L, R, w, v̄, d)           # Alg 2 — OUTER loop only, no nesting
   → Ā_i, π_ni, π_nn, P_n
   ↓
land_price / invert_amenities               # added: Q_n, H̄_n, B_ni  (for §4.2 and §5)
```

Neither inversion is nested: each is a single damped fixed point. That is *different* from the rest
of this repository (AA2022, FFW, Redding-2016 all nest a wage loop inside a population loop) and is
worth stating in the header — the nesting appears only in the **forward** solver I add (§4.2) and
implicitly in the counterfactual solver (§4.3).

### 4.2 Forward equilibrium solver (`solve_equilibrium`) — **added, not upstream**

The toolkit has no forward solver at all. It is needed for the round-trip validation (§5) and is the
piece the repo's house style expects. Nested damped fixed point:

```
outer: guess w (unit mean), L, Q
  1. A_i   = Ā_i L_i^ν                                       (Codebook A.1 eq 9)
  2. π_ni  = L_i (d_ni w_i / A_i)^{1-σ} / Σ_k (...)          Paper (6)
  3. P_n   = σ/(σ-1) (L_n/(σ F π_nn))^{1/(1-σ)} d_nn w_n/A_n Paper (8)
  4. λ_ni  = B_ni (κ_ni P_n^α Q_n^{1-α})^{-ε} w_i^ε / Φ      Paper (10)
     L_i   = L̄ Σ_n λ_ni ; R_n = L̄ Σ_i λ_ni                  Paper (11)
     v̄_n   = Σ_i λ_ni w_i / λ^R_n                            Paper (14)
  5. Q_n   = ((1-α) v̄_n R_n / H̄_n)^{1/(1+δ)}                Paper (5) + SW2020 (7)
  6. w_i   ← Σ_n π_ni v̄_n R_n / L_i, renormalized to unit mean   Paper (7)
  damp (w, L, Q) by ζ; converge when max|Δ| < tol
```

### 4.3 Counterfactual solver (`counterfactual`, Alg 11)

Single outer fixed point on `{ŵ, λ̂}` with the seven other updaters evaluated in sequence inside it
(order is load-bearing and is taken verbatim from `counterFactsTK.m:89-109`):

```
C1 v̂ → C2 L̂ → C3 R̂ → C4 P̂_H → C5 π̂ → C6 P̂_Q → C7 w̃ (then renormalize w'=ŵ·w to unit mean) → C8 λ̃
convergence: max|ŵ-w̃| < 1e-4 AND max|λ̂-λ̃| < 1e-4
damping:     ŵ ← 0.25 w̃ + 0.75 ŵ ;  λ̂ ← 0.25 λ̃ + 0.75 λ̂
```

### 4.4 Orientation conventions (the #1 silent-error risk)

Fixed for the whole port and asserted at runtime:

* `lambda[n, i]` — residence `n` (row) × workplace `i` (column). `sum(lambda) == 1`.
* `lambda_cond[n, i]` — conditional on residence; **rows** sum to 1.
* `pi[n, i]` — destination/consumer `n` (row) × origin/producer `i` (column). **Rows** sum to 1.
  This is the paper's `π_ni` in Eq (6).
* `d[n, i]`, `kappa[n, i]` — same (destination, origin) / (residence, workplace) convention.

**Upstream differs and it matters.** `solveProductTradeTK.m` builds `tradesh` in the `[n,i]`
orientation inside the loop (line 66: `repmat(num', nobs, 1)`, `sum(...,2)`) but returns the
**transpose** `[i,n]` from the convergence block (line 93: `repmat(num, 1, nobs)`, `sum(...)` over
dim 1). Every counterfactual updater (`updateTradeshTK`, `updateWageTK`) then consumes the `[i,n]`
orientation. The two agree **only because `dni` is symmetric** — the convergence block applies
`dni.^(1-σ)` without transposing it, so with an asymmetric trade-cost matrix the returned `tradesh`
would be wrong. See decision D5.

---

## 5. Validation strategy with numeric targets

| # | Check | Target | How |
|---|---|---|---|
| V1 | Conditional commuting probabilities sum to 1 per residence row | `max|rowsum-1| < 1e-12` | assertion |
| V2 | Unconditional `λ_ni` sums to 1 overall | `|Σλ - 1| < 1e-12` | assertion |
| V3 | Labour market clearing `Σ_i L_i == Σ_n R_n == L̄` | `< 1e-9` relative | assertion; `L̄ = 401` for mean-normalized data |
| V4 | Trade shares sum to 1 per destination row | `max|rowsum-1| < 1e-12` | assertion |
| V5 | `getBi` matches observed workplace employment | upstream rule `Σ|L_obs-L_pred|·100 < 0.001`, i.e. `Σ|Δ| < 1e-5` over 401 counties | solver stopping rule |
| V6 | `getBi` fixed point independent of damping | `B_i` from `damp=1.0` vs `damp=0.5` agree to `< 1e-8` relative | run both. NOTE: the shipped default is now `damp=0.5` (Codebook) with `maxiter=5000`; `damp=1.0` reproduces `getBiTK.m`. Damping needs 1,043 iterations against the full step's ~520, so the old `maxiter=1000` was too tight for it |
| V7 | Productivity inversion: income == expenditure | `round(|inc-exp|, 6) == 0` (upstream rule) | solver stopping rule |
| V8 | **Levels round-trip**: forward-solve from inverted `{Ā, B_ni, H̄}` starting at `w=L=R=Q=1` | recover data `L, R, w, v̄, Q` to `< 1e-6` max relative deviation | `solve_equilibrium` |
| V9 | **Zero-shock exact-hat**: all shocks = 1 | every hat `== 1` to `< 1e-6`; welfare change `0.00%` | `counterfactual` |
| V10 | Welfare is location-invariant | `max/min` over all `J²` cells of `welfChange` agree to `< 1e-5` relative | Codebook A.2 |
| V11 | Commuting gravity elasticity of **model-predicted** flows, with residence+workplace FE | exactly `-εμ = -2.162` (`< 1e-6`) | `gravity_fe`; this is an identity, so a miss = orientation bug |
| V12 | Trade gravity elasticity of **model** trade shares, with FE | exactly `ψ(1-σ) = -1.26` | same |
| V13 | Commuting gravity elasticity of **observed** German flows | compare to MRRH2018 US estimate `-4.43` (their φ) | `gravity_fe` on `commuting_wide.csv` |
| V14 | Trade gravity, structural | `-1.26` vs MRRH2018's CFS estimate `-1.29` | report both |
| V15 | **Headline welfare**: κ̂ = 0.88 off-diagonal, 1 on-diagonal | paper Table 5 col 2: **+3.26%** | `counterfactual` |
| V16 | Full Table 5: κ̂ ∈ {0.79, 0.88, 0.96, 1.13} | paper: `+6.89, +3.26, +0.89, -2.33` % | 4 runs |
| V17 | Monotonicity: lower commuting cost ⟹ higher welfare | sign check across V16 | — |
| V18 | Uniform κ̂ on **all** pairs including diagonal | analytic: `Û = 1/κ̂` exactly, `ŵ = λ̂ = P̂ = 1` | closed-form check of the solver |

**V15/V16 caveat to state up front:** MRRH2018's 3.3% is US counties (N≈3,111) with *their*
parameters (§3.2). This port runs German counties (N=401) with the *toolkit's* SW2020 parameters
(§3.1). The experiment is reproduced structurally, not numerically. To separate the two sources of
difference the port runs V16 twice: once with §3.1 parameters and once with §3.2 parameters on the
same German data. Whatever comes out is reported as-is.

---

## 6. Known upstream defects and how each is handled

| # | Defect | Location | Decision |
|---|---|---|---|
| 1 | `psi` switched 0.42 → 0.21 mid-script | `Counterfactuals.m:139` (restored at `:221`) | Reading the full script, this is **not a bug**: it is a deliberate "half the trade cost" sensitivity run, bracketed by `save`/re-`OwnData`/restore. Implement **ψ = 0.42 as the baseline**; ship ψ = 0.21 as an explicitly named low-trade-cost experiment. Comment the upstream bracketing and the interruption hazard (an aborted run leaves `data/output/*.mat` in the ψ=0.21 state). |
| 2 | GRID track uses ν = 0, county track ν = 0.05 | `GRID_MRRH2018_toolkit.m:55` vs `MRRH2018_toolkit.m:55` | `nu` is an explicit field with no default inheritance; both values named and documented; the paper's own ν = 0 noted too. |
| 3 | `getBiTK.m` has no damping despite Codebook Alg 1 step 5 | `getBiTK.m:70` | Add `damp` kwarg. Default **1.0** = upstream's full step (exact reproduction); verify (V6) that `damp = 0.5` reaches the same fixed point. |
| 4 | `counterFactsTK.m` outer loop is `while true` | `counterFactsTK.m:87` | Add `maxiter` (default 5000) and `@warn` on non-convergence, returning the last iterate rather than hanging. |
| 5 | `solveProductTradeTK.m` returns the transposed trade-share matrix | `:93,99` vs `:66,69` | Silent only because `dni` is symmetric. Port uses one orientation `π[n,i]` throughout (§4.4) and asserts row sums. |
| 6 | `updateEmplTK`/`updateResidentsTK`/`updateWageTK` normalizations commented out | `:28-31`, `:28-31`, `:36-39` | Left out, matching the shipped code. Noted in comments. |

---

## 7. MATLAB → Julia traps to watch (checklist)

* `sum(M)` / `mean(M)` etc. sum **columns** in MATLAB, **everything** in Julia → always `dims=`.
  Concretely: `getBiTK.m:68` is `sum(...,2)` (rows); `updateEmplTK.m:23` is `sum(...,1)`;
  `updateResidentsTK.m:23` is `sum(...,2)`; `solveProductTradeTK.m:69` is `sum(...,2)` but `:99` is
  bare `sum(nummat)` = **dim 1**; `updateTradeshTK.m:41` bare `sum(...)` = dim 1;
  `updateLamTK.m:37` uses `sum(X(:))` = scalar.
* `1./x` → `1 ./ x`; `A(i,j)` → `A[i,j]`; `x(x>0)=5` → `x[x .> 0] .= 5`.
* `diag(M)` extracts in MATLAB; use Julia `diag`, never `diagm`. Used at
  `solveProductTradeTK.m:105`, `updatePricesTK.m:25` (twice).
* `repmat(v, 1, n)` → `repeat(v, 1, n)` (v down rows); `repmat(v', n, 1)` → `repeat(v', n, 1)`
  (v across columns). Getting this backwards silently transposes the model.
* `csvread(f,1,1)` skips 1 header row **and** 1 index column.
* Variables assigned only inside `while`/`for` are loop-local in Julia → pre-declare.
* Convergence by rounding, not tolerance: `round(x, digits=6) == 0`.

---

## 8. Decisions taken where the source was ambiguous

* **D1 — Primary entry point is the OwnData pathway.** `ReadData.m` needs an observed bilateral
  commuting matrix and rationalizes its 150,907 zeros (of 160,801 cells) with κ = ∞, which is
  inconsistent with any smooth travel-cost matrix and blocks extensive-margin counterfactuals.
  `OwnData.m` needs only `{L_i, R_n, w_i or a rent index, area, κ_ni}`. The port's public API is the
  OwnData one; the German data is the worked example, and a second country is a matter of swapping
  five arrays.
* **D2 — Wages.** The German data has observed workplace wages, so `w` is passed in. The port also
  exposes the toolkit's no-wage fallback (feed `w = ones(J)`, recover `w = B_i^{1/ε}`, ARSW2015
  transformed-wage interpretation) as a documented keyword.
* **D3 — Commuting cost `κ_ni` = straight-line distance in km**, following `OwnData.m:41`, not the
  `roundtrip_time_base.csv` travel times (which carry 1e6 sentinels for zero-flow pairs and belong
  to the `ReadData` track). Only the *relative* κ matters: λ_{ni|n} row-normalizes, so any uniform
  rescaling of κ cancels.
* **D4 — `L̄` is `sum(L_obs)` after mean-normalization, i.e. 401.** Upstream normalizes `L, R, w` to
  unit mean, so `L̄ = J`. Kept, because `updateEmplTK` depends on `L̄` and `L_obs` being on the same
  scale (`L̂ = L̄ Σ_n λ_ni λ̂_ni / L_i` equals 1 at `λ̂ = 1` only under that convention).
* **D5 — Single trade-share orientation `π[n,i]`** throughout (§4.4), correcting the upstream
  transposition that is invisible for symmetric `d`.
* **D6 — Land prices and full `B_ni` are computed** (Q11/Q12), which upstream never does. `Q_n` needs
  a land endowment: use `H_n = Area_n` (geographic land area, MRRH2018's baseline interpretation),
  then back out `H̄_n = H_n / Q_n^δ`. `B_ni` is identified only up to scale, so `Φ` is normalized to
  1. Neither object changes any upstream result — both are additions that make the forward solver and
  the round-trip test possible. The observed rent index (`house_prices.csv`) is *not* used to pin
  `Q_n`, matching upstream, where it is descriptive only.
* **D7 — Welfare is reported as the mean over all `J²` cells** of `welfChange` with the max spread
  printed, rather than upstream's `welfChange(1,1)`. Algebraically identical at convergence; the
  spread is a free convergence diagnostic (V10).
* **D8 — `Γ((ε-1)/ε)` is dropped** from the welfare expression, as upstream does: it cancels in the
  ratio `Û = Ū'/Ū`. `SpecialFunctions.gamma` is used only to report the level constant.
* **D9 — Plots.** Germany has no shapefile in this folder (upstream's `MAPIT.m` needs the Mapping
  Toolbox and `shape/VG250_KRS_clean_final`), so choropleths are replaced by the repo's usual
  scatter/histogram/line diagnostics keyed on border distance and county area. No new data is added.

---

## 9. Post-implementation record

Written after the port ran. `monte_redding_rossihansberg_commuting_model.jl`, 1,325 lines,
runs clean with `julia <file>.jl` in about 70 s, exits 0, emits no warnings, writes six PDFs
to `graphs/`, and passes **21 of 21** validation checks.

### 9.1 What actually happened, section by section

Everything in §1-§8 was implemented as planned. Six things changed or were added.

**[CHANGED] The plan's §5 target list grew from 18 to 21 checks.** Added during implementation:

* **V8b — local uniqueness.** V8 (round-trip from a cold start) turned out to be cheap, so a
  50% log-normal perturbation of `(w, L)` was added. It returns to the same fixed point to
  3.9e-8, which is evidence for local uniqueness that the toolkit never provides (upstream
  verifies no uniqueness condition anywhere).
* **V19 — local employment elasticity.** The plan tracked only the welfare headline. The
  paper's *other* headline (Section III: the elasticity ranges 0.5-2.5 across US counties,
  mean 1.52; residents 0.2-1.2) is checkable and was added, using the paper's own 5% shock.
* **V20 — the own-commuting-share regression.** Paper Table 2 column 5 reports `R^2 = 0.89`
  from regressing the employment elasticity on `lambda^R_ii|i` alone. Added.

**[CHANGED] V6's comparison needed a normalisation the plan missed.** `B_i` is identified
only up to scale — `lambda_ni|n` row-normalises, so a uniform rescaling of `B` cancels. A
first run compared raw `B` from `damp = 1.0` and `damp = 0.5` and showed a 2e-2 discrepancy
that was purely a difference in scale, not in the fixed point. `getBi` now normalises `B` to
unit mean before returning, and the check passes at 1.6e-8.

**[CHANGED] The forward solver needed two changes the plan did not anticipate**, and this was
the only part of the port that genuinely fought back. §4.2's plain damped Jacobi iteration
**diverged** from a uniform start (after 20,000 iterations `L` was off by a factor of 479).
Diagnosis and fix:

1. A one-step residual evaluated *at the data* came back at 1e-15, proving the data was an
   exact fixed point and the problem was the iteration, not the inversion.
2. `Q` must be initialised from land-market clearing, not at 1. Equilibrium `Q` has a mean of
   order 1e-3 on this data because land area is in km², so starting at 1 is three orders of
   magnitude away. This was the larger of the two effects.
3. Damping must be applied **in logs**. With `eps = 4.6` the map `lambda ~ w^eps` is explosive
   and linear damping on the level does not contain it.

With both fixes it converges in 521 iterations from a cold start and recovers the data to
1.4e-7, better than the 1e-6 target. Both are commented in the function's docstring.

**[DONE, with a correction] Defect 5 was confirmed in the code, not just in the plan.** The
first counterfactual prototype failed the zero-shock test spectacularly (`max|lambda_hat - 1|`
= 402) because `update_wage` broadcast the residence-indexed vector `v_hat.*r_hat.*v.*R`
across *columns* instead of *rows* — the exact `[n,i]` vs `[i,n]` confusion the plan predicted
would be the main hazard. It is fixed, and V11/V12 exist specifically to catch a recurrence:
they are identities that can only hold if every matrix is oriented correctly.

**[ADDED] Defect 7, not in the plan's §6 list.** `solveProductTradeTK.m:107` computes the
price index as `(L_n/(sigma F pi_nn))^(1/(1-sigma)) * w_n/Abar_n`, dropping **both** the
`L_n^(-nu)` agglomeration correction and the `d_nn` term that MRRH2018 eq. (8) and Codebook
Algorithm 8 both carry. On this data that understates `P_n` by up to a factor of 2.86,
because `diag(d)` ranges from 1 to 2.88. It is harmless upstream — `P_n` is only ever mapped
in `Descriptives.m` and never enters a counterfactual, and `updatePricesTK.m:25` (the hat
version) has both terms and is correct. It is **not** harmless here, because `P_n` enters the
forward solver through `lambda_ni`; the round-trip fails without the correction. The correct
expression is used and `price_index_upstream` is kept so the run can print the comparison.

**[CHANGED] The ReadData track was promoted from "not implemented" to a shipped second
track.** The plan (decision D1) made OwnData the primary entry point, which it remains. But
MRRH2018 itself quantifies from observed flows, so a ReadData track is the closer analogue of
the paper and costs about fifteen lines. Shipping both turned the headline check from a
single number that misses into a **bracket** that contains the paper's value, and made it
possible to attribute the difference.

### 9.2 Validation outcomes against the §5 targets

| # | Target | Achieved | |
|---|---|---|---|
| V1 | `max|rowsum-1| < 1e-12` | 2.7e-15 | PASS |
| V2 | `|sum-1| < 1e-12` | 0.0 | PASS |
| V3 | `< 1e-9` relative | 1.4e-16, `Lbar` = 401 | PASS |
| V4 | `max|rowsum-1| < 1e-12` | 2.4e-15 | PASS |
| V5 | upstream rule `< 1e-3` | 9.9e-4 in 519 iterations; max rel dev 7.9e-8 | PASS |
| V6 | `< 1e-8` relative | 1.6e-8 (after the normalisation fix above) | PASS |
| V7 | 6-dp rounding rule | 4.7e-7 in 71 iterations | PASS |
| V8 | `< 1e-6` | **1.4e-7** in 521 iterations from a cold start | PASS |
| V8b | — | 3.9e-8 after a 50% log-normal perturbation | PASS |
| V9 | `< 1e-6`, welfare 0.00% | **2.6e-8, welfare +0.0000%**, converges in 1 iteration | PASS |
| V10 | `< 1e-5` relative spread | 3.2e-8 across all 160,801 cells | PASS |
| V11 | exactly `-2.162` | **-2.16200000, difference 0.0** | PASS |
| V12 | exactly `-1.26` | **-1.26000000, difference 2.2e-16** | PASS |
| V13 | compare to US `-4.43` | -2.087 (n = 9,894, within-R² 0.793), implying `mu` = 0.4536 against the calibrated 0.47 | PASS |
| V14 | `-1.26` vs `-1.29` | as planned | PASS |
| V15 | paper `+3.26%` | **+2.32% (OwnData) / +5.18% (ReadData)** — brackets it | PASS |
| V16 | paper `+6.89, +3.26, +0.89, -2.33` | OwnData `+5.22, +2.32, +0.63, -1.40`; ReadData `+11.34, +5.18, +1.43, -3.27` | PASS |
| V17 | monotone | strictly decreasing in `kappa_hat` on both tracks | PASS |
| V18 | `U_hat = 1/0.88` exactly | 1.13636364 vs 1.13636364, difference 2.0e-10 | PASS |
| V19 | paper `[0.5, 2.5]`, mean 1.52 | **`[1.12, 2.18]`, mean 1.79** with the paper's parameters (`[2.82, 4.20]` with the toolkit's) | PASS |
| V20 | paper `R² = 0.89` | **0.851** with the paper's parameters (0.277 with the toolkit's) | PASS |

### 9.3 What did NOT reproduce, stated plainly

**The 3.3% welfare headline is not reproduced numerically, and could not have been.** The
paper's number is 3,111 US counties with MRRH2018's parameters; this is 401 German counties.
What the port establishes instead:

* The two German quantification tracks bracket it: **+2.32% and +5.18% against +3.26%**.
* The difference is **data, not parameters**. Re-running the ReadData track with MRRH2018's
  own `{alpha, epsilon, mu, delta, nu, psi}` moves +5.18% only to +4.85%.
* The direction of each track is understood and was checked:
  - ReadData lands *above* because German counties are slightly **more** open to commuting
    than US ones — median `lambda_nn|n` of 0.643 against the US 0.69 in 2000 (paper Fig. 1) —
    so a proportional cut in off-diagonal costs mechanically buys more. (An earlier draft of
    this document asserted the opposite, that German counties are more self-contained; that
    was wrong, and was caught by computing the median conditional share rather than the
    unconditional diagonal mass, which are different statistics — 0.643 vs 0.671.)
  - OwnData lands *below* because its predicted flows are far too concentrated on the
    diagonal: median 0.837 against the observed 0.643. Straight-line distance with a very
    small within-county diagonal makes counties look more closed than they are. This is a
    property of the toolkit's **default cost matrix**, not of the port — a user with real
    travel times passes them to `build_costs`.
* The residual is general-equilibrium dampening. In MRRH2018's US model the partial-equilibrium
  arithmetic (`Phi_hat^(1/eps)` holding everything else fixed) gives roughly +4.7% and the
  full solve gives +3.26%, a 30% dampening; on German data the same comparison gives +4.96%
  and +4.85%, almost none. That difference plausibly reflects 401 vs 3,111 locations and the
  associated housing-price response, but **it cannot be decomposed without MRRH2018's own
  data** and is not claimed here.

**MRRH2018's commuting gravity slope of -4.43 is not matched** (-2.087 on observed German
flows, same fixed effects). This is a property of the two countries' geography, not a porting
error. The check that this data *does* support — whether the toolkit's `mu = 0.47` is
consistent with the flows it is applied to — passes: the observed slope implies `mu = 0.4536`.

### 9.4 Decisions revisited

D1-D9 all stand as written, with three amendments:

* **D1 amended**: OwnData remains the primary entry point, but ReadData is now shipped
  alongside it rather than merely described. See §9.1.
* **D6 amended**: the plan said `B_ni` is identified up to the scale of `Phi`; `B_i` from
  Algorithm 1 is *also* identified only up to scale, and is now normalised to unit mean.
* **D9 amended**: the plot list grew to six figures, adding the Table 5 comparison against
  the paper's own points and the local-elasticity distribution.

One new decision:

* **D10 — the forward solver damps in logs and initialises `Q` from land-market clearing.**
  Not a modelling choice (the fixed point is unchanged) but a numerical one, without which
  the round-trip test cannot run at all. Documented in `solve_equilibrium`'s docstring.

### 9.5 Left undone

* **Upstream's choropleths.** `MAPIT.m`/`GRIDMAPIT.m` need the MATLAB Mapping Toolbox and the
  `shape/VG250_KRS_clean_final` shapefile, which is not vendored in this folder. Replaced by
  scatter/histogram diagnostics keyed on border distance, area and own-commuting share
  (decision D9). No spatial data was added to the repository to work around this.
* **The GRID track pipeline.** `GRIDData.m`/`GRIDCounterfactuals.m` depend on three other
  Ahlfeldt toolkits (AABPL → GRID → TTMATRIX) and a Python pre-processing chain that produces
  files not present here. The only part of the GRID track that matters for the model itself —
  its `nu = 0` parameter set — is shipped as `GRID_PARAMS` and its effect on the headline
  counterfactual is reported (+2.29% against +2.32%).
* **`progs/old/`**, excluded deliberately, as instructed.
* **Local employment elasticities are computed on a 50-county subsample**, not all 401, to
  keep the script's runtime near a minute. Each county costs one counterfactual solve at
  `tol = 1e-9`. Raising `nsub` to 401 is a one-line change and takes about eight minutes.
* **No uniqueness condition is verified analytically.** Upstream verifies none either, and
  MRRH2018's Propositions B.1/B.2 establish uniqueness of the *inversions*, not of the forward
  equilibrium under SW2020's agglomeration. V8b (perturb and return) is numerical evidence,
  not a proof. Anyone raising `nu` materially above 0.05 should treat convergence failure as
  a signal rather than a nuisance.
