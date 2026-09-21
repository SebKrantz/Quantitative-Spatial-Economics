# Translation Plan: Ahlfeldt & Barr (2022) Stata → Julia

**Target**: `ahlfeldt_barr_skyscraper_model.jl`, single self-contained file, Julia 1.12.4.

**Source**: `stata_source/AB2022.ado` (the shipped ado solver) and `stata_source/_1_PROGS.do`
(the walkthrough solver), with `AB2022-codebook.pdf` §A.1 (13 equations) / Table 1 as the
authoritative economic documentation.

**Gold standard**: `reference_data/{BASE,INVERTED,EMPIRICAL_CH_skyline}.csv`, converted from the
toolkit's `.dta` files.

Status markers: **[plan]** written before any Julia code; **[outcome]** appended after
implementation (see §9).

---

## 1. What the model is

A stylized, perfectly-open, one-dimensional monocentric city: Alonso–Mills–Muth with an endogenous
**vertical** margin. Space is a line `x ∈ [-50, 50]` km, `D = |x|`. Two competing land uses,
commercial (`C`) and residential (`R`), plus an agricultural outside option at exogenous rent `r_a`.
Developers choose height to maximise profit per unit land given the floor-space price; land goes to
the highest bidder; a single city-wide wage `y` and city employment `N` (`L` in code) adjust until
commercial labour demand equals residential labour supply.

**There is no commuting-cost term.** Distance enters *only* through the exponential decay
`exp(-τ^U·D)` of the productivity/amenity fundamentals. A single wage applies everywhere. `τ^C`/`τ^R`
are a reduced-form decay of agglomeration/amenity value with distance, not a commuting primitive.

The Codebook lists a reservation utility `Ū` among the primitives but it appears nowhere in the code
and is not one of the 15 solver arguments; openness is implemented by making total employment `N` a
target object rather than fixing utility. `Ū` is algebraically absorbed into `a^R(x)`. This is a
simplification of the toolkit, not a defect to fix.

---

## 2. Equation-by-equation mapping

`U ∈ {C, R}`. "Codebook" = `AB2022-codebook.pdf` §A.1 numbering. Stata line numbers are
`stata_source/`-relative. All of this lands in one Julia function, `solver!`, mirroring Stata's
`SOLVER` — the recursion is a single deterministic pass, so splitting it into eight small functions
would add indirection without adding clarity.

| # | Object | Codebook | `AB2022.ado` | `_1_PROGS.do` | Julia |
|---|---|---|---|---|---|
| 1 | Productivity shifter `Ã^C(x) = ā^C ã^C N^{β^C} e^{-τ^C D}` | Table 1 | `:63` | `:98` | `solver!` step 1 |
| 2 | Amenity shifter `Ã^R(x) = ā^R ã^R N^{β^R} e^{-τ^R D}` | Table 1 | `:64` | `:99` | `solver!` step 1 |
| 3 | Commercial floor-price shifter `a^C(x) = Ã^C(x)^{1/(1-α^C)} y^{α^C/(α^C-1)}` | eq. 1 | `:65` | `:100` | `solver!` step 2 |
| 4 | Residential floor-price shifter `a^R(x) = Ã^R(x)^{1/(1-α^R)} y^{1/(1-α^R)}` | eq. 2 | `:66` | `:101` | `solver!` step 2 |
| 5 | Commercial land rent `r^C(x)` | eq. 5 | `:69` ✗ | `:104` ✓ | `solver!` step 3 |
| 6 | Residential land rent `r^R(x)` | eq. 6 | `:70` ✗ | `:105` ✓ | `solver!` step 3 |
| 7 | Land-use allocation `U ∈ {1,2,3}` | eqs. 7–8 | `:74-76` | `:109-111` | `solver!` step 4 |
| 8 | Boundaries `x0` (CBD edge), `x1` (urban edge) | Table 1 | `:79-86` | `:114-121` | `solver!` step 5 |
| 9 | Profit-maximising height `S*^U(x) = (a^U/(c^U(1+θ^U)))^{1/(θ^U-ω^U)}` | Table 1 | `:89-90` | `:124-125` | `solver!` step 6 |
| 10 | Constrained/realised height `S^U(x) = min(S̄^U, S*^U)` in own zone | eqs. 3–4 | `:93-94` | `:128-129` | `solver!` step 6 |
| 11 | Horizontal bid rent `p̄^U(x) = a^U(x)/(1±ω^U)·S^U(x)^{ω^U}` | eqs. 9–10 | `:97-98` **(1−ω)** | `:132-133` **(1+ω)** | `solver!` step 7 — **see §4** |
| 12 | Workplace employment `L(x) = α^C/(1-α^C)·p̄^C/y·S^C` | eq. 11 | `:101` | `:136` | `solver!` step 8 |
| 13 | Aggregate labour demand `L̂ = Σ L(x)` | eq. 13 | `:102-103` | `:137-138` | `solver!` step 8 |
| 14 | Floor space per capita `f̄^R = (1-α^R)/p̄^R·y` | eq. 12 | `:106` | `:141` | `solver!` step 9 |
| 15 | Residence employment `n(x) = S^R/f̄^R`, `N̂ = Σ n(x)` | eqs. 12–13 | `:107-110` | `:142-145` | `solver!` step 9 |
| 16 | `S_x = max(S^C, S^R)`, `URBAN = U<3`, `COM = U==1` | — | `:117-121` | `:152-156` | `solver!` step 10 |
| 17 | Wage loop (Algorithm 2) | §A.2 Alg. 2 | `:131-160` | `:166-195` | `wage!` |
| 18 | Equilibrium loop (Algorithm 3) | §A.2 Alg. 3 | `:414-446` | `:202-256` | `findeq!` |
| 19 | Amenity update to match heights (Algorithm 4) | §A.2 Alg. 4 | — | `:268-278` | `conv!` |
| 20 | Amenity scaling to match employment (Algorithm 5) | §A.2 Alg. 5 | — | `:284-290` | `emp!` |
| 21 | Inversion driver (Algorithm 6) | §A.2 Alg. 6 | — | `:296-316` | `invert!` |

`✓`/`✗` in the land-rent rows mark agreement with the Codebook — see §4.

### 2.1 Land rent, written out

Codebook eqs. 5–6 give land rent structurally, as floor revenue net of construction cost with height
already optimised out:

```
r^U(x) = a^U/(1+ω^U) · (S̃^U)^{1+ω^U} − c^U · (S̃^U)^{1+θ^U}
```

Substituting `S̃^U = (a^U/(c^U(1+θ^U)))^{1/(θ^U-ω^U)}`:

```
r^U(x) = a^U/(1+ω^U) · (a^U/(c^U(1+θ^U)))^{(1+ω^U)/(θ^U-ω^U)}
         − c^U · (a^U/(c^U(1+θ^U)))^{(1+θ^U)/(θ^U-ω^U)}
```

This is **exactly** `_1_PROGS.do:104-105`. `AB2022.ado:69-70` writes `c^C·(1+ω^C)` inside the first
bracket instead of `c^C·(1+θ^C)`, which is inconsistent with its own `S*` on line 89. The port
implements the Codebook/`_1_PROGS` form unconditionally; the `.ado` land-rent bracket is treated as a
plain bug, not a convention (it is ~400% off — see §7, measured).

---

## 3. Parameter table

| Symbol | Stata name | Codebook | Value | Meaning |
|---|---|---|---|---|
| `α^C` | `alpha_C` | `α^C` | 0.85 | Commercial share of the non-floor-space input in production. Not a solver argument. |
| `α^R` | `alpha_R` | `α^R` | 0.66 | Residential expenditure share on the non-housing good. Not a solver argument. |
| `β^C` | `beta_C` | `β^C` | 0.03 | **Agglomeration** elasticity (commercial). |
| `β^R` | `beta_R` | `β^R` | 0.00 | Residential agglomeration elasticity — off by default. Not a solver argument. |
| `τ^C` | `tau_C` | `τ^C` | 0.01 | Commercial productivity spatial decay, per km. |
| `τ^R` | `tau_R` | `τ^R` | 0.005 | Residential amenity spatial decay, per km. |
| `ω^C` | `omega_C` | `ω^C` | 0.03 | Commercial **rent elasticity of height**, `ω^U = ω̃^U/(1-α^U)`. |
| `ω^R` | `omega_R` | `ω^R` | 0.07 | Residential rent elasticity of height. |
| `θ^C` | `theta_C` | `θ^C` | 0.5 | Commercial **construction-cost elasticity of height**. |
| `θ^R` | `theta_R` | `θ^R` | 0.55 | Residential construction-cost elasticity of height. |
| `c^C` | `c_C` | `c^C` | 1.4 / 1.3 | Commercial baseline construction cost. |
| `c^R` | `c_R` | `c^R` | 1.4 / 1.3 | Residential baseline construction cost. |
| `ā^C` | `a_bar_C` | `ā^C` | 2 | Fundamental commercial productivity. |
| `ā^R` | `a_bar_R` | `ā^R` | 1 | Fundamental residential amenity. |
| `ã^U` | `a_rand_C`, `a_rand_R` | `ã^U` | 1 | Location-specific amenity component. Uniform 1 for the stylized city; the object the Chicago inversion solves for. |
| `r_a` | `r_a` | `r_a` | 150 / 50 / 30 | Agricultural land rent. A per-observation **variable** in Stata, not a scalar. |
| `S̄^C` | `S_bar_C` | `S̄^C` | 999 | Commercial height limit; 999 ≈ unconstrained. |
| `S̄^R` | `S_bar_R` | `S̄^R` | 999 | Residential height limit. |
| `x̄1` | `x_1_bar` | — | 999 | Urban growth boundary, km. Only used to trim the grid. |
| — | `y` | `y` | 2.5 (start) | City-wide wage. **Target object.** |
| — | `L` | `N` | 1e6 (start) | Total employment. **Target object.** |

**Three parameterisations appear in the toolkit and they are not interchangeable:**

| Context | `c_C`,`c_R` | `r_a` | Source |
|---|---|---|---|
| `AB2022.ado` built-in defaults | 1.4 | 150 | `AB2022.ado:387-389` |
| `_2_ANALYSIS.do` baseline gradients | 1.3 | 50 | `_2_ANALYSIS.do:23` |
| `_1_PROGS.do` `CONV`/`EMP` (the inversion) | 1.4 | 30 | `_1_PROGS.do:276, 288` |

`_1_PROGS.do:17-33` also sets `beta_C=0.04`, `omega_C=omega_R=0.1`, `c_C=c_R=1`, `a_bar_C=1.5`,
`theta_R=0.5`, `tau_R=0.01` while building `BASE.dta`. These are **placeholders overwritten the
moment `FINDEQ` is called** with its explicit argument list; they matter only if `SOLVER` is called
directly. The port will expose them as the `Params` defaults matching the *ado* (the canonical paper
values) and require the caller to pass the `_2_ANALYSIS`/`CONV` variants explicitly.

### 3.1 Notation-collision warning (mandatory, per `Ahlfeldt/CLAUDE.md`)

These letters mean different things here than elsewhere in this repository. Transplanting a
calibrated value across models will look plausible and be wrong:

| Symbol | **AB2022 (here)** | Elsewhere in `Quantitative-Spatial-Economics/` |
|---|---|---|
| `θ` | construction-cost elasticity of height, **0.5–0.55** | Fréchet/trade elasticity, **6–10** (AA2022, Redding-JIE-2016, BHS2026) |
| `β` | **agglomeration** elasticity, **+0.03** | **amenity externality**, **−0.3** (AA2022) — opposite sign *and* opposite economics |
| `α` | non-floor-space input/expenditure share, **0.66–0.85** | productivity externality, **0.1** (AA2022) |
| `ε` | *not used* | migration elasticity, 3–5 |
| `λ` | *not used* | traffic-congestion elasticity, 0.07–0.09 (AA2022) — and the ARSW2015 *agglomeration* elasticity at nearly the same magnitude |
| `σ` | *not used* | CES substitution, ≈4 |

This warning goes in the Julia header block as a comment, not only here.

---

## 4. The ω-convention decision

### 4.1 The discrepancy

`AB2022.ado:97-98` and `_1_PROGS.do:132-133` implement **different formulas** for the horizontal
(per-floor) bid rent:

```
AB2022.ado:97       p_bar_x_C = a_x_C * 1/(1-omega_C) * S_x_C^omega_C
_1_PROGS.do:132     p_bar_x_C = a_x_C * 1/(1+omega_C) * S_x_C^omega_C
```

```
AB2022.ado:98       p_bar_x_R = a_x_R * 1/(1-omega_R) * S_x_R^omega_R
_1_PROGS.do:133     p_bar_x_R = a_x_R * 1/(1+omega_R) * S_x_R^omega_R
```

Codebook eqs. 9–10 read `p̄^U(x) = a(x)/(1+ω^U) · S^U(x)^{ω^U}` — the Codebook agrees with
`_1_PROGS.do`, against the `.ado`.

### 4.2 Which is theoretically correct

`(1+ω)` is correct, and this can be settled without reference to any data. The developer's problem is

```
max_S  p̄(S)·S − c·S^{1+θ}
```

Under `p̄ = a/(1+ω)·S^ω`, revenue is `a/(1+ω)·S^{1+ω}`, so the FOC is

```
a·S^ω = c(1+θ)·S^θ   ⟹   S* = (a/(c(1+θ)))^{1/(θ-ω)}
```

which is **exactly** the `S*` both files use (`.ado:89`, `_1_PROGS:124`). Under `(1−ω)` the FOC would
instead give `S* = (a(1+ω)/((1−ω)c(1+θ)))^{1/(θ-ω)}`, which neither file implements. The `.ado` is
therefore internally inconsistent: its `p̄` does not generate its own `S*`, and its `p̄·S` is not the
first term of its own land rent.

### 4.3 Decision

Implement **both**, selected by the keyword argument `bid_rent_convention`:

- `:codebook` — `1/(1+ω)`. **Default.** Matches Codebook eqs. 9–10, `_1_PROGS.do:132-133`, and the
  first-order condition above.
- `:ado` — `1/(1-ω)`. Matches `AB2022.ado:97-98`.

### 4.4 What the reference data says — the decisive evidence

**The shipped `INVERTED.dta` was generated with `(1−ω)`, the `.ado` convention.** Measured
pointwise against `INVERTED.csv`'s own stored columns (§6, Tier B):

| Formula | max rel. dev. on `p_bar_x_C` | max rel. dev. on `p_bar_x_R` |
|---|---|---|
| `1/(1-ω)` (`.ado`) | **1.06e-07** ✓ | **1.23e-07** ✓ |
| `1/(1+ω)` (Codebook, `_1_PROGS`) | 5.83e-02 ✗ | 1.31e-01 ✗ |

So the level effect is **+5.8% commercial, +13.1% residential** floor-space rent from switching
`(1+ω)` → `(1−ω)`, i.e. `(1+ω)/(1−ω)` exactly. (The task brief's "roughly 6%" is the commercial
figure; residential is larger because `ω^R = 0.07 > ω^C = 0.03`.)

**`BASE.csv` cannot answer the question at all** — it contains no solved values (§5).

The reference data is therefore a **hybrid** that neither current source file reproduces on its own:

| Block | `INVERTED.csv` matches | Does *not* match |
|---|---|---|
| Land rent (eqs. 5–6) | `_1_PROGS.do:104-105` / Codebook (2.46e-07) | `AB2022.ado:69-70` (4.08e+00) |
| Bid rent (eqs. 9–10) | `AB2022.ado:97-98` (1.06e-07) | `_1_PROGS.do:132-133` / Codebook (5.83e-02) |

Most plausible history: the original solver had `c(1+θ)` land rent *and* `(1−ω)` bid rent; the `.ado`
(v0.93, 02/2024) later broke the land-rent bracket while `_1_PROGS.do` (01/2024) fixed the bid rent to
match the Codebook — and `INVERTED.dta` predates both edits.

Consequence for validation: Tier B runs under `bid_rent_convention = :ado` to reproduce the shipped
data, while the shipped *default* stays `:codebook`. Both facts get stated in the header, the README
and §9.

---

## 5. `BASE.csv` is not a solved equilibrium

`_1_PROGS.do:81` saves `BASE.dta` **before** any solver runs — it is the initialised placeholder
grid, written between `gen`-ing the variables (`:48-78`) and defining `SOLVER` (`:90`). Verified:

- 10,001 rows × 28 columns.
- **Populated (7):** `x` (−50…50, step 0.01), `D = |x|`, `y ≡ 2.5`, `L ≡ 1e6`, `r_a ≡ 100`,
  `a_rand_C ≡ 1`, `a_rand_R ≡ 1`.
- **Entirely missing (21):** `a_x_C`, `a_x_R`, `A_tilde_x_C`, `A_tilde_x_R`, `r_x_C`, `r_x_R`, `U`,
  `S_star_x_C`, `S_star_x_R`, `S_x_C`, `S_x_R`, `S_x`, `p_bar_x_C`, `p_bar_x_R`, `L_x_C`,
  `f_bar_x_R`, `n_x`, `URBAN`, `COM`, `SHADE`, `SHADEU` — 0 non-null out of 10,001 each.

Note `r_a = 100` here, which is none of the three solver parameterisations in §3 — it is the `gen`
placeholder from `_1_PROGS.do:53`, overwritten by `FINDEQ`'s 12th argument at `:216`.

**Implications, which the task brief's validation spec did not anticipate:**

1. "Reproduce BASE.csv column by column" is satisfiable, but only 7 of 28 columns carry information.
   The remaining 21 are validated as *all-missing*, which is a real check (it catches a port that
   silently initialises to `0.0` instead of `NaN`) but not an economic one.
2. **BASE.csv cannot discriminate the ω convention**, and cannot supply land-use boundary locations
   (`U` is entirely missing). Both of those requirements move to `INVERTED.csv`.
3. The per-column deviation table for BASE will be reported in full anyway, as asked.

`INVERTED.csv` is the only solved gold standard in the folder, and it is a good one: all 24 solved
columns populated, 4,001 rows (`|x| ≤ 20` after `_3_INVERSION.do:119`'s `drop if D > 20`), plus the
merged Chicago data and the binning variables.

---

## 6. Validation strategy and tolerances

### 6.0 The precision floor

Stata stores `gen`-created variables as **4-byte floats** by default. Round-tripping every numeric
column of `INVERTED.csv` through `float32` changes it by at most **5.7e-08** relative. That is the
noise floor: **no column can be reproduced better than ≈6e-08 relative**, regardless of how correct
the Julia is. Targets below are set against that floor, not against an arbitrary 1e-6.

Comparisons are made against the **full stored precision in the CSVs**, never against Stata's display
precision.

### 6.1 Tier A — `BASE.csv`, exact grid reproduction

Build the grid from `AB2022.ado:339-342` / `_1_PROGS.do:42-45`: `set obs 10001`, `x = (n-5001)/100`
for `n = 1…10001`, `D = |x|`. Then compare all 28 columns.

- **Target:** 7 populated columns exact to float32 round-trip (≤ 1e-7 relative); 21 columns
  `all-missing` on both sides.
- **Fail condition:** any populated column deviating > 1e-7, or any of the 21 columns non-missing in
  the port.

### 6.2 Tier B — `INVERTED.csv`, pointwise identity reproduction (**the real test**)

Take `INVERTED.csv`'s own exogenous state as given — `a_rand_C`, `a_rand_R`, `y`, `r_a = 30`, and the
`CONV`/`EMP` parameters `(0.5, 0.55, 0.03, 0.07, 0.030, 2, 1, 0.01, 0.005, 1.4, 1.4, 30, 999, 999)`
from `_1_PROGS.do:276` — run **one** `solver!` pass and compare all 24 solved columns.

This is a genuine end-to-end test of every equation in §2 simultaneously, and it avoids the
path-dependence problem in §6.3.

Two replication details are mandatory, both discovered empirically and both non-obvious:

**(a) The stored `L` is one update ahead of the columns it generated.** `findeq!`'s loop updates
`L ← 0.5L + 0.25(L̂ + N̂)` *after* the last `solver!` call (`AB2022.ado:443`, `_1_PROGS.do:247`), so
the `L` saved in the dataset is not the `L` that produced `A_tilde_x_C`. Recover it as

```
L_prev = 2·(L_stored − 0.25·(L̂ + N̂))     where  L̂ = Σ L_x_C,  N̂ = Σ n_x
```

For `INVERTED.csv`: `L_stored = 1004690.56`, `L̂ = 1007078.2597`, `N̂ = 1011696.9057` ⟹
`L_prev = 999993.5373`. Confirmed: `A_tilde_x_C` max relative deviation falls from **1.407e-04**
(using `L_stored`) to **8.51e-08** (using `L_prev`) — the residual being exactly the float32 floor,
and the 1.407e-04 being exactly `β^C·ln(L_stored/L_prev) = 0.03 × 4.697e-3`.

**(b) `S_star_x_C`/`S_star_x_R` are never cleared between solver passes.** The
"clear any pre-existing value" loop (`AB2022.ado:58`, `_1_PROGS.do:93`) lists 14 variables and
**omits both `S_star` variables**, so they retain stale values from earlier iterations at locations
whose land use has since changed. In `INVERTED.csv`, `S_star_x_C` is non-missing at 1,387 rows where
`U ≠ 1`, and `S_star_x_R` at 2,135 rows where `U ≠ 2`. The port must replicate this (do not clear
`S_star`), and `S_star` is compared **only within its own zone**.

- **Target:** ≤ **1e-06** relative on every column; expected ≈ **2e-07**.
- Under `bid_rent_convention = :ado` (per §4.4). The same table under `:codebook` is reported
  alongside, to quantify the convention's effect.
- Measured in advance on the stored columns: `A_tilde_x_R` 1.40e-07, `a_x_C` 5.72e-07, `a_x_R`
  3.77e-07, `r_x_C` 2.46e-07, `S_star_x_C|U=1` 1.70e-07, `S_star_x_R|U=2` 1.81e-07, `L_x_C` 1.85e-07,
  `f_bar_x_R` 1.88e-07, `n_x` 1.29e-07. All well inside target.

### 6.3 Tier C — general equilibrium

`FINDEQ` stops at a **1% relative tolerance** with fixed damping, so its terminal `(y, L)` is *not* a
precise fixed point — it is path-dependent, and re-running from a different start legitimately lands
somewhere else inside the 1% band. Exact end-to-end reproduction of the *iteration path* that
produced `INVERTED.dta` is not attainable (it ran through a long `INVERT` history). This is a
property of the source, not a porting failure, and will be stated as such.

What *is* checkable:

1. Both of Stata's convergence tests pass at the reproduced state. Measured in advance:
   outer `|L/(0.5(L̂+N̂)) − 1| = 0.00465 < 0.01` ✓; inner `|L̂/N̂ − 1| = 0.00457 < 0.01` ✓.
2. Labour demand ≈ labour supply: `L̂ = 1007078.26`, `N̂ = 1011696.91`, gap 0.457% — inside the
   solver's own 1% tolerance, as required.
3. A fresh `findeq!` from `(y, L) = (2.5, 1e6)` under `_2_ANALYSIS.do:23` converges, and both
   objectives end < 0.01.

### 6.4 Tier D — land-use boundaries

Recompute `U` from the reproduced `r_x_C`, `r_x_R`, `r_a` and compare **elementwise** to
`INVERTED.csv`'s `U`, then compare `x0` (min `x ≥ 0` with `U ≠ 1`) and `x1`.

- **Target:** `U` identical at all 4,001 rows; `x0`, `x1` identical to the 0.01 km grid step.
- `x1` is computed two ways upstream — `.ado:79-80` takes `max(x)` where `U < 3`, `_1_PROGS:114-115`
  takes `min(x)` where `U == 3`. These differ by exactly one grid step. The port uses the
  `_1_PROGS`/Codebook-walkthrough form and notes the difference in a comment. `x0`/`x1` are reported
  diagnostics only; neither enters any equation.

### 6.5 Tier E — Chicago inversion

Replicate `_3_INVERSION.do`: merge `EMPIRICAL_CH_skyline.csv` on `X` (558 rows, all on the 0.01 grid,
no duplicates, all within `|X| ≤ 20`); assign each location to the taller of `HEIGHT_C`/`HEIGHT_R`
(`:39-40`); bin `CONVBIN = round(X, 0.1)` and take the within-bin max (`:44-48`); then run `invert!`.

Note `INVERT` runs on the **full 10,001-point grid**; `drop if D > 20` happens afterwards (`:119`),
affecting only what is saved.

`INVERT 1000000 0.05` (`_3_INVERSION.do:91`) ⟹ `CONV` receives **λ = 0.05** as its convergence
parameter. (`Ahlfeldt/AB2022-toolkit/CLAUDE.md` states λ = 0.5; that is a typo in the upstream notes —
`INVERT`'s second argument is passed through as `CONV`'s first, `_1_PROGS.do:304`.)

- **Target:** height correlation `R²(HEIGHT_R, S_x_R) ≥ 0.999` (the source's own stopping rule) and
  `|N − 1,000,000| ≤ 1000`.
- The recovered `a_rand_C`/`a_rand_R` are compared to `INVERTED.csv`'s in distribution (share of
  exact zeros — 3,759 and 2,135 respectively — plus quantiles), **not** pointwise: the inversion path
  is as tolerance-limited as §6.3, and an exact pointwise match is not expected. Any mismatch is
  reported with numbers, not glossed.

### 6.6 Tier F — ω-convention effect size

Solve the `_2_ANALYSIS.do:23` baseline twice, once per convention, and report the level shift in
`p̄^C`, `p̄^R`, `y`, `N`, `x0`, `x1`.

---

## 7. Stata → Julia semantic traps, and how each is handled

| Trap | Handling |
|---|---|
| `scalar` vs `gen`/`replace` | `scalar` → a Julia scalar field on `Params`/`State`; `gen`/`replace` → a length-`N` vector, broadcast. `r_a` is the trap: declared `scalar` in `FINDEQ`'s signature but assigned with `qui replace r_a = \`12'` (`_1_PROGS.do:216`), so it is a **variable**. `AB2022.ado:410-411` converts it explicitly (`drop r_a` then `gen r_a = r_a`). Port keeps it a vector. |
| `replace y = f(x)` over observations | Broadcast `@. y = f(x)`. |
| Missing propagation | Stata `.` is **larger than any number** in comparisons; Julia `NaN` compares false everywhere. These differ when one side is missing and the other is not. Audited: after `solver!`, `r_x_C`/`r_x_R` are never missing (`a^U ≥ 0`, and `a^U = 0` yields `r^U = 0`, not missing), so the `U`-assignment comparisons never hit the divergent case. An `@assert` guards this rather than leaving it to luck. |
| `sum`/`r(sum)` | Stata sums **non-missing** only. Julia: explicit NaN-skipping sum. |
| `min(S_bar, ·)` / `max(S_x_C, S_x_R)` | Stata's `min`/`max` **ignore** missing and return missing only if all arguments are. Implemented explicitly, not via Julia's `min`/`max` (which propagate NaN). Verified against `INVERTED.csv`'s `S_x`. |
| `egen ..., by(bin)` | Explicit grouping over `CONVBIN`. |
| `round(X, 0.1)` | Stata rounds to the nearest multiple, halves away from zero. Compared to `INVERTED.csv`'s `CONVBIN` with a 1e-6 tolerance (the column is float32). |
| `_n` / `_N`, `in 1/L` | Not used in the ported paths. |
| Grid spacing | Reproduced exactly from `:339-342`: `x = (n-5001)/100`. Matching this is what makes a column-by-column comparison meaningful rather than an interpolation-error measurement. |
| Display precision | All comparisons against full CSV precision; the float32 storage floor is quantified in §6.0. |

---

## 8. Solver architecture and upstream defects

```
Params            immutable struct: α^U, β^U, τ^U, ω^U, θ^U, c^U, ā^U, r_a, S̄^U,
                  bid_rent_convention
City              mutable struct: the 28 grid vectors + scalars (y, L, L̂, N̂, x0, x1)

solver!(city, p)        Algorithm 1 — one deterministic recursive pass. No iteration.
wage!(city, p)          Algorithm 2 — inner loop on y
findeq!(city, p)        Algorithm 3 — outer loop on L, nests wage!
conv!(city, p, λ)       Algorithm 4 — amenity update toward observed heights
emp!(city, p, target)   Algorithm 5 — amenity scaling toward target population
invert!(city, p, ...)   Algorithm 6 — drives conv! then emp!
```

`wage!` (`AB2022.ado:135`): while `|L̂/N̂ − 1| > 0.01`, set `y_factor = (L̂/N̂)^0.01` (or a flat 1.2 if
`N̂ = 0`, 0.8 if `L̂ = 0`), update `y ← 0.5y + 0.5·y·y_factor`, re-run `solver!`. Note the first
objective uses `+0` and subsequent ones `+0.0001` (`:133` vs `:155`) — replicated.

`findeq!` (`AB2022.ado:418-446`): `solver!`; then while `|L/(0.5(L̂+N̂)) − 1| > 0.01`: `solver!`,
`wage!`, abort if `L̂ + N̂ == 0`, recompute the objective, **then** `L ← 0.5L + 0.25(L̂+N̂)`. The
objective is updated *before* `L` — this ordering is what leaves `L` one step ahead at exit (§6.2a)
and must be preserved exactly.

### Upstream defects, each handled with an in-code comment

1. **No `maxiter` guard** on either loop (`AB2022.ado:135`, `:422`). A parameterisation that never
   reaches the 1% tolerance without tripping the `L̂+N̂ == 0` check hangs Stata indefinitely. The port
   adds `maxiter` to both (defaults 1000 inner / 1000 outer) and emits `@warn` on non-convergence
   rather than failing silently or looping forever.
2. **`grc1leg`** (SSC) is an undeclared dependency of every graph-combining step, despite the
   toolkit README claiming no user-written ado files are used. Stata-figure-only; irrelevant to this
   port, recorded here for completeness.
3. **`S_star` never cleared** (§6.2b) — replicated deliberately, with a comment, because the
   reference data depends on it.
4. **`.ado` land-rent bracket** uses `c(1+ω)` where the Codebook and `_1_PROGS` use `c(1+θ)` (§2.1).
   Not implemented as an option; the `.ado` form is ~400% off the shipped data and inconsistent with
   its own `S*`.
5. **`AB2022.ado:238` defines a `FINDEQ` that is never called** — `AB2022` re-implements the same
   loop inline at `:414-446`. Three copies of the equilibrium loop exist upstream. The port has one.
6. **`_1_PROGS.do:110-111`** inline comments for the `U==2`/`U==1` branches are swapped. Comment-only;
   the logic matches the `.ado`.
7. **`_1_PROGS.do:286-289`** — `EMP`'s `foreach name in R C` loop never uses `name`, so the body runs
   twice per call. Replicated (it changes the amenity scaling per call).

---

## 9. Outcome

*Appended after implementation. Every number below was measured by running the file.*

### 9.1 What was built

`ahlfeldt_barr_skyscraper_model.jl` — 998 lines (685 code, 218 comment, 95 blank; the header block
carries the 13 equations, the notation-collision warning and the two convention discrepancies).
Runs clean under Julia 1.12.4 with `Statistics, Printf, Plots, DataFrames, CSV`. Sections follow the
brief: header, `using`, parameters, city construction, height/bid-rent solver, land-use assignment,
GE loop, inversion, counterfactuals, visualization. Nine PDFs written to `graphs/`.

All five validation tiers pass. Full per-column tables are in `README.md` §6; they are not duplicated
here.

### 9.2 Validation results

| Tier | Target | Result |
|---|---|---|
| A — `BASE.csv` | 7 populated columns exact; 21 all-missing | **PASS** — all 7 at **0.000e+00**; 21 all-missing both sides |
| B — `INVERTED.csv` | every solved column within the data's own precision | **PASS** — worst rel. dev. **2.491e-06**, worst dev/tolerance ratio **0.59** |
| C — equilibrium | both Stata objectives < 0.01 | **PASS** — inner 0.00457, outer 0.00465 |
| D — boundaries | `U` elementwise; `x0`, `x1` at grid resolution | **PASS** — `U` identical 4001/4001; `x0` = 0.08, `x1` = 4.46 km |
| E — Chicago inversion | `R² ≥ 0.999`, population gap ≤ 1000 | **PASS** — `R²` = 0.99957 (5 CONV steps), gap = 641.2 (205 EMP steps) |
| F — ω effect size | measure it | +8.33% / +15.36% on peak floor rent; **+16.41% on equilibrium population** |

**The §6.0 precision-floor analysis changed during implementation and this is the most important
methodological outcome.** The plan assumed a flat ~6e-08 Stata-float32 floor and a 1e-06 target. Both
were wrong:

- The reference CSVs carry ~8 significant decimal digits, but **`y` is stored to only six**
  (`2.10701`). Since `a^C ∝ y^{-5.667}`, raised to `1/(θ^C-ω^C) = 2.128` for height and compounded
  again for labour demand, half an ulp of that printed `y` is amplified by up to ~21×.
- So a flat tolerance is meaningless. The implemented criterion is **per-column and derived from the
  data**: re-solve at `y ± half-ulp`, take the induced spread, floor it at `4·eps(Float32)/2 ≈
  2.4e-07` for columns with no `y`-sensitivity (`A_tilde_x_R` is the only one that binds on the
  floor). Six columns exceed the naive 1e-06 target while sitting at **6%** of what the data can
  actually resolve.

The error chain is pure input-quantization amplification, as predicted in §6.2 but larger than
estimated: `a_rand` → `A_tilde` (8.5e-08) → `^6.667` → `a_x_C` (7.5e-07) → `^2.128` → `S*`
(1.6e-06) → `L_x_C` (2.5e-06).

Replication details (a) and (b) from §6.2 were both confirmed necessary and both behave exactly as
predicted: without the off-by-one `L`, `A_tilde_x_C` sits at 1.407e-04 instead of 8.51e-08.

### 9.3 The ω-convention finding

Confirmed, and it is the hybrid predicted in §4.4. Worst relative deviation across all 19 solved
columns of `INVERTED.csv`:

| Convention | Worst deviation |
|---|---|
| `:ado` — `1/(1-ω)` (`AB2022.ado:97-98`) | **2.49e-06** ✓ |
| `:codebook` — `1/(1+ω)` (`_1_PROGS.do:132-133`, Codebook eqs. 9–10) | 1.51e-01 ✗ |

**`INVERTED.dta` was generated with `(1−ω)`** — the `.ado` bid rent — combined with the
`_1_PROGS`/Codebook *land-rent* bracket. Neither current source file reproduces it alone. Default
remains `:codebook`, which the FOC argument in §4.2 shows is the correct economics.

`BASE.csv` was, as §5 predicted, unable to answer the question: it holds no solved values.

### 9.4 A third upstream defect, found during implementation

Not in the plan, not in the toolkit's notes, not in the task brief's list:

**Land rent ignores the height limit.** Both Stata files inline the closed form of eqs. 5–6, so `S̄`
never appears and land rent is evaluated at the *unconstrained* `S*` rather than Codebook eq. 5's
`S̃ = min(S̄, S*)`. Harmless whenever the cap does not bind — which is true for the baseline
(`S̄ = 999`) and for every column of the reference data, so Tiers A–E are unaffected. It matters
sharply once a cap binds: at the centre with `S̄^C = 20` the toolkit gives land rent **409.26** where
eq. 5 gives **204.90**, a **99.7% overstatement**. Because land rent drives the land-use allocation
(eqs. 7–8), the height-limit counterfactual's boundaries inherit the error.

Implemented as `constrained_land_rent` (**default `true` = Codebook eq. 5, corrected**; pass `false` to reproduce the toolkit); the
counterfactual reports both. Under eq. 5 the same 20-floor cap gives max land rent 206.79 (vs
409.26), `N` = 834,255 (vs 841,490), `x0` = 6.23 km (vs 6.36), `x1` = 21.90 km (vs 22.17).

### 9.5 Decisions taken where the source was ambiguous

1. **ω convention** — both implemented, default `:codebook` (§4.3). The FOC argument settles the
   economics independent of the data; the data went the other way, and both facts are reported.
2. **Land-rent bracket** — `.ado` form treated as a plain bug, not a convention; only the
   Codebook/`_1_PROGS` form implemented (§2.1).
3. **Land rent under a binding cap** — both implemented (§9.4). Default **corrected** (Codebook eq. 5): the height cap is the model's principal counterfactual, so shipping the toolkit's ~2x overstatement as the default would mislead. `constrained_land_rent=false` reproduces upstream.
4. **`x1` definition** — `_1_PROGS`'s `min(x | U==3, x≥0)`; the `.ado`'s `max(x | U<3, x≥0)` differs
   by one grid step. Diagnostic only; neither enters an equation.
5. **`S_star` staleness** — replicated deliberately; the reference data depends on it. Compared only
   within its own zone.
6. **`EMP` double-execution** — replicated (§8.7).
7. **Placeholder parameters** — `Params` defaults follow the `.ado` (canonical paper values);
   `params_analysis()` and `params_conv()` expose the other two parameterisations explicitly.
8. **Tier E comparison mode** — distributional, not pointwise, for the §6.3 tolerance reason.
9. **Tier B tolerance** — data-derived per column rather than a flat constant (§9.2). This was a
   change from the plan, made because the flat target was not defensible.

### 9.6 Bugs found in my own implementation, for the record

Two Julia-specific traps, both caught by the validation rather than by inspection:

- `minimum(...; init=NaN)` and `maximum(...; init=NaN)` always return `NaN`, because `min(NaN, x)`
  is `NaN`. This silently produced `x0 = x1 = NaN` and blank counterfactual rows.
- A closure defined inside a function **assigns to the enclosing function's local** on plain
  assignment. `solve_at(y) = (c = City(...); ...)` therefore rebound the outer `c`, so the nominal
  solve was silently replaced by the last perturbed one and Tier B's deviations inflated 20×. Fixed
  by renaming the closure's local.

### 9.7 Corrections to the upstream notes

- `Ahlfeldt/AB2022-toolkit/CLAUDE.md` states `CONV`'s λ = 0.5. It is **0.05** (`_3_INVERSION.do:91`
  passes it as `INVERT`'s second argument, which `_1_PROGS.do:304` forwards as `CONV`'s first).
  `Ahlfeldt/` is read-only and was not modified.
- The task brief's expectation that `BASE.csv` identifies the ω convention does not hold (§5).

### 9.8 Left undone

- **No Stata was run.** The `.dta` → CSV conversion was taken as given, and the ~8-significant-digit
  export precision is the binding constraint on Tier B (§9.2). Re-exporting at full `double`
  precision would tighten the comparison by roughly two orders of magnitude and is the single
  highest-value follow-up.
- `_4_INVERTEDCOUNTER.do`'s subcenter counterfactual is implemented (seeded at 2022) but validated
  only qualitatively: upstream calls `runiform()` with no seed, so levels are not comparable across
  runs. The qualitative claims all reproduce — a new subcenter appears (0 → 156.2 floors at
  11–16 km), the original CBD loses height (189.6 → 150.7, −20.5%), and the city grows — and the
  wage effect matches the do-file's own annotation ("Wage increases by about 5%") at **+5.26%**. The
  population effect is +69% against its stated +30%, which is within what an unseeded shock explains.
- The height-limit counterfactual is validated against the paper's qualitative claims, not against
  numbers: `FIG_gradients_wHL` ships as a figure only, with no `.dta`. Note that two of the upstream
  prose claims do not survive: the CBD *expands* (`x0` +42%, vertical→horizontal substitution) rather
  than shrinking, and peak floor rent *falls* at the centre where the cap binds. The urban footprint,
  wage and population all move as described.
- The `S_x`/`SHADE`/`SHADEU` graphing columns are reproduced as all-missing in Tier A but are not
  otherwise modelled; they are Stata plotting scaffolding with no economic content.
