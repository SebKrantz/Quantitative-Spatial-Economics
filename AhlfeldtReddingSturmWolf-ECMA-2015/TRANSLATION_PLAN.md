# TRANSLATION_PLAN.md — ARSW (2015) MATLAB toolkit to Julia

Port of Gabriel M. Ahlfeldt's didactic MATLAB toolkit for

> Ahlfeldt, G. M., Redding, S. J., Sturm, D. M., Wolf, N. (2015).
> *The Economics of Density: Evidence from the Berlin Wall.* Econometrica 83(6), 2127–2189.

into a single self-contained Julia file, `ahlfeldt_redding_sturm_wolf_density_model.jl`,
running on **synthetic** Berlin-like data.

Sections 1-12 were written **before** any Julia code was produced. Section 13 ("Post-mortem")
was appended after implementation and records what actually happened, including the places
where the plan turned out to be wrong.

---

## 1. Why synthetic data, and why that makes the test sharper

The toolkit needs three `.mat` files (`prepdata_big_TD.mat`, `prepdata_big_TD86.mat`,
`ttpublic_2006_ren.mat`) that are external downloads from an HU-Berlin file server and are not
in the checkout. They are deliberately **not** downloaded here. Consequently:

* `prepdata_TD.m`, `prepdata_TD86.m`, `prepdata_TD06ttpub.m` are **read for schema only**, never
  ported. They only rebuild those `.mat` files.
* The paper's empirical numbers (λ = 0.0710, ε ≈ 6.69–6.83, ν = 0.07/0.0987, …) are used as
  **parameter inputs**, never as targets. Nothing in this port claims to reproduce Berlin.
* `MAPIT.m` (Berlin shapefile choropleths, Mapping Toolbox) and `modbezirk.m` (hard-coded 23→12
  historic Bezirke crosswalk) are out of scope.

The compensation is that synthetic data enables a **stronger** test than matching Berlin numbers.
Supplement §S.3.1.6 (Proposition S.3) proves a one-to-one mapping

```
{α, β, μ, ε, κ}  +  observed {Q, H_M, H_R, K, τ}   ->   {Ã, B̃, φ̃}
```

and §S.3.2 (Proposition S.4) extends it with {λ, δ, η, ρ} to {ã, b̃}. So if we *plant* known
fundamentals, solve the model *forward* to manufacture "observed" data, and then run the toolkit's
*inversion* on that data, we must get the planted fundamentals back. Any deviation is an
implementation bug, with no data-quality excuse. That is the primary validation in this port.

---

## 2. Notation collisions — read before naming anything

ARSW shares letters with the other models in this repository and means different things by them.
Every symbol below is documented again in the Julia file's header block, with an explicit warning.

| Symbol | **ARSW2015 meaning (this file)** | Same letter elsewhere in this repo |
|---|---|---|
| `lambda` λ | **Agglomeration (production-externality) elasticity**, 0.0710. *Higher λ raises productivity.* | AllenArkolakis-RES-2022: traffic **congestion** elasticity 0.07–0.09. *Higher λ raises transport cost.* Nearly identical magnitude, opposite economics. **Never transplant.** |
| `epsilon` ε | Fréchet shape for bilateral **commuting** preference shocks, 6.83 / 6.694 | Redding-JIE-2016 / AA2022: migration elasticity 3–5 |
| `alpha` α | **Labour** share in production, 0.80 (floor-space share in firm costs is 1−α = 0.20) | AA2022: productivity externality 0.1 |
| `beta` β | Expenditure share on the **tradable numeraire**, 0.75 (floor-space share 1−β = 0.25) | AA2022: amenity externality −0.3 |
| `theta` θ | **Commercial floor-space share** of a block, θ_i ∈ [0,1] | Redding/AA2022: trade / Fréchet elasticity 6–10 |
| `pi_ij` π_ij | **Commuting probability** (live i, work j) | Redding-JIE-2016: bilateral **trade** share |
| `mu` μ | Non-land input share in floor-space production | — |
| `kappa` κ | Iceberg commuting-cost parameter, κ = ν/ε | — |

Note especially that α is the **labour** share, not the floor-space share: the production function
is `Y_j = A_j H_Mj^α (θ_j L_j)^{1-α}` and the paper sets "the share of firm expenditure on
commercial floor space (1−α) equal to 0.20" (p. 2167). The upstream `ARSW2015-toolkit/CLAUDE.md`
parameter table calls α the "floor-space input share", which is the complement of the truth;
`cftualprep_end_TD.m` itself comments `alpha=0.80; % Set input share of labour in production`,
which is correct.

---

## 3. The model, compactly

Index set: `N` blocks, each simultaneously a possible workplace `j` and a residence `i`.

**Commuting / residence choice.** Workers draw i.i.d. Fréchet(ε) shocks over (residence,
workplace) pairs:

```
(4)     φ_ij = exp(-ε κ τ_ij) · B̃_i^ε · Q_i^{-(1-β)ε} · w̃_j^ε ,   π_ij = φ_ij / Φ ,  Φ = Σ_ij φ_ij
(5)     π_Ri = Σ_j π_ij ,      π_Mj = Σ_i π_ij ,        H_Ri = π_Ri H ,  H_Mj = π_Mj H
(6)     π_{j|i} = φ_ij / Σ_j' φ_ij'
(S.20)  E[w|i] = Σ_j π_{j|i} w̃_j          (total worker income vv_i = E[w|i]·H_Ri)
(9)     Ū = γ · Φ^{1/ε},   γ = Γ((ε-1)/ε)
```

**Production, land market.**

```
(10)  Y_j   = A_j · H_Mj^α · (θ_j L_j)^{1-α}
(12)  q_j   = (1-α) (α/w_j)^{α/(1-α)} A_j^{1/(1-α)}        (zero profit; inverted for w or A)
(13)  no-arbitrage θ_j ∈ {0} / [0,1] / {1} as q_j ≶ ξ_j Q_j
(S.30) θ_j L_j       = ((1-α) A_j / q_j)^{1/α} H_Mj        (commercial land clearing = Eq. 18)
(S.29) (1-θ_i) L_i   = (1-β) E[w|i] H_Ri / Q_i             (residential land clearing = Eq. 19)
(S.31) L_i = φ_i K_i^{1-μ}                                 (= Eq. 15, density of development)
```

**Inversion equations actually coded.**

```
(S.44)  H_Mj = Σ_i [ ω_j e^{-ν τ_ij} / Σ_s ω_s e^{-ν τ_is} ] H_Ri ,   ω_j = w̃_j^ε ,  ν = κ ε
(S.46)  W_i  = Σ_s (w̃_s / e^{κ τ_is})^ε           (residential commuting market access, CMA)
(S.47)  B̃_i/B̄ = (H_Ri/H̄_R)^{1/ε} · (Q_i/Q̄)^{1-β} · (W_i/W̄)^{-1/ε}
(S.48)  Ã_j   = (q_j/(1-α))^{1-α} (w̃_j/α)^α
```

**Endogenous agglomeration (Section 7).**

```
(20)  A_j = a_j · Υ_j^λ ,   Υ_j = Σ_s e^{-δ τ_js} (H_Ms / K_s)
(21)  B_i = b_i · Ω_i^η ,   Ω_i = Σ_s e^{-ρ τ_is} (H_Rs / K_s)
```

---

## 4. Algorithm map: codebook number → supplement equation → MATLAB file → Julia function

`matlab_source/` filenames are the flattened, section-prefixed copies in this folder.

| Alg. | Equations | MATLAB file | Julia function | Status |
|---|---|---|---|---|
| 1 | S.44 | `section6_optimepsilon_comegaoptO.m` | `solve_transformed_wages` | CORE |
| 2 | S.64 (GMM moment) | `section6_optimepsilon_cdensityoptren.m` | `epsilon_objective` | STRETCH |
| 3 | — (patternsearch) | `section6_optimepsilon_optimepsilon_TD86.m` | `estimate_epsilon` (Optim.jl) | STRETCH |
| 4 | S.44, S.48 | `section6_calibration_comegaoptC.m` | `solve_wages_and_productivity` | CORE |
| 5 | S.46, S.47 | `section6_calibration_camen.m` | `solve_amenities` | CORE |
| 6 | (12), Φ=H | `section6_calibration_calcal_adj_TD.m` | `rescale_fundamentals` | CORE |
| 7 | S.20 | `section6_calibration_expincome.m` | `expected_income` | CORE |
| 8 | S.29–S.31 | `section6_calibration_cdensity.m` | `solve_density` | CORE |
| 9 | (4),(5),(12) fixed point | `section6_exogcftual_cmodexog.m` | `invert_simultaneous` | CORE |
| 10 | S.29–S.31 (+ levels) | `section6_exogcftual_cdensityE.m` | `solve_density` (`levels=true`) | CORE |
| 11 | (4),(5),(9),(10),(12),(18),(19) | `section6_exogcftual_smodexog.m` | `solve_equilibrium_exog` | CORE |
| 12 | (20) / S.55 | `section7_counterfactual_cprod.m` | `decompose_productivity` | CORE |
| 13 | (21) / S.56 | `section7_counterfactual_cres.m` | `decompose_amenity` | CORE |
| 14 | (9) | `section7_counterfactual_ubar.m` | `reservation_utility` | CORE |
| 15 | as 11 + (20),(21) | `section7_counterfactual_smodendog.m` | `solve_equilibrium_endog` (closed city) | CORE |
| 16 | as 15, Ū fixed, H free | `section7_counterfactual_ussmodendog.m` | `solve_equilibrium_endog` (`closed=false`) | CORE |
| — | drivers | `calcal_TD`, `cftualprep_TD`, `cftualexog_TD`, `cftualprep_end_TD`, `cftualendog_*` | inlined in the `main()` driver | CORE |
| — | data prep | `prepdata_TD*.m` | **read for schema only, not ported** | OUT |
| — | mapping | `MAPIT.m`, `modbezirk.m` | **out of scope** | OUT |

Algorithms 15 and 16 are merged into one Julia function with a `closed::Bool` keyword, because
`ussmodendog.m` is `smodendog.m` plus one extra update line for total employment `H`; keeping two
near-identical 400-line copies is exactly the duplication the port should remove.
Algorithms 8 and 10 likewise merge (`cdensityE.m` = `cdensity.m` returning `LM`, `LR` as well).

---

## 5. Parameter table

| Symbol | Julia name | Section 6 value | Section 7 (GMM) value | Meaning |
|---|---|---|---|---|
| α | `alpha` | 0.80 | 0.80 | **Labour** share in production (1−α = 0.20 floor space in firm costs) |
| β | `beta` | 0.75 | 0.75 | Expenditure share on the tradable good (1−β = 0.25 on residential floor space) |
| μ | `mu` | 0.75 (see note) | 0.75 | Non-land input share in floor-space production; `L = φ K^{1-μ}` |
| ε | `epsilon` | 6.83 | 6.694 | Fréchet shape, commuting |
| ν = κε | `kappaeps` | 0.07 | 0.0987 | Commuting decay per minute |
| κ | `kappa` | derived `= kappaeps/epsilon` | derived | Iceberg commuting cost |
| λ | `lambda` | 0 (no spillovers) | 0.0710 | Density elasticity of **productivity** |
| δ | `delta` | 0 | 0.3617 | **Productivity**-spillover distance decay |
| η | `eta` | 0 | 0.1553 | Density elasticity of residential **amenity** |
| ρ | `rho` | 0 | 0.7595 | **Amenity**-spillover distance decay |
| γ | derived | `Γ((ε-1)/ε)` | same | Scale term in expected utility, Eq. (9) |

**μ note (upstream defect 3).** μ is never a named variable in MATLAB; it appears only as the
literal exponent `0.75` in `cdensity.m:57`, `cdensityE.m:61`, `smodendog.m:85` and
`ussmodendog.m:95`. Read as `K^{1-μ}` that implies μ = 0.25, but the paper (p. 2167) sets "the
share of land in construction costs (1−μ) equal to 0.25", i.e. μ = 0.75 and the exponent should be
**0.25**. In the Julia port `mu` is a named parameter and the exponent is written `1 - mu`.
This is numerically inert everywhere in the toolkit — `cdensity` divides by `K^{1-μ}` and the
solvers immediately multiply by `K^{1-μ}` again — and the port will demonstrate that
numerically rather than assert it.

---

## 6. Synthetic-data design

### 6.1 Geography

* `nx × ny = 15 × 15 = 225` blocks on a 1 km grid (indexing column-major, block (c,r) → i).
  225 keeps every dense `N×N` operation trivial while still giving a genuine spatial structure.
* CBD at the grid centre.
* Block land area `K_i` (km²): smaller near the CBD (finer subdivision), lognormal noise.
  `log K_i = log K0 + g_K · (d_i/d_max) + N(0, σ_K²)`, `g_K > 0`.
* Bilateral travel times `τ_ij` in **minutes** (the unit ν = 0.07 is calibrated in):
  Euclidean distance / 0.5 km·min⁻¹ (30 km/h), with `τ_ii = 1.5` min so own-block commuting is
  not free.

### 6.2 The wall

A vertical cut between columns 7 and 8 splits the city into West (cols 1–7) and East (cols 8–15).

* **Divided regime**: any origin–destination pair straddling the cut has
  `τ_ij += WALL_PENALTY` with `WALL_PENALTY = 300` min. With ν = 0.07 that multiplies φ_ij by
  `e^{-21} ≈ 7.6e-10` — commuting across the wall is severed but the matrix stays strictly
  positive, so no zeros are introduced into the log-based diagnostics.
* **Reunified regime**: penalty 0.

This gives the paper's identifying variation a synthetic analogue: the same fundamentals
{a, b, φ} under two different τ matrices. Used for the comparative-statics checks and for the
"division / reunification" counterfactual.

### 6.3 Planted fundamentals

Drawn once with `Random.seed!(20150615)` (paper's journal year + issue; any fixed seed is fine —
the point is reproducibility).

* `log a_j = g_a · (d_j / d_max) + σ_a ζ_j` — fundamental productivity, *rising* with distance from
  the CBD (`g_a > 0`), matching the paper's finding that *total* A falls with distance while the
  *exogenous* component rises.
* `log b_i = g_b · (d_i / d_max) + σ_b ξ_i` — fundamental amenity, also rising with distance.
* `log φ_i = -g_φ · (d_i / d_max) + σ_φ ψ_i` — density of development falling with distance.
* Total employment `H = 1.0e6`.
* `ζ, ξ, ψ` are independent standard normals, so a and b are not mechanically collinear.

Two block-composition scenarios:

* **`:mixed`** (primary) — every block has `a_j > 0` and `b_i > 0`, so every block is incompletely
  specialised, `q_i = Q_i`, and the single "observed floor price" per block is unambiguous. This is
  the benchmark case of supplement §S.2.6.
* **`:specialised`** (secondary) — a central "office district" (a 3×3 patch) has `b_i = 0`
  (θ = 1, purely commercial, `H_R = 0`), and two "park/housing" patches have `a_j = 0` (θ = 0,
  purely residential, `H_M = 0`). This exercises the `IcsA` / `IcsB` / `Iis` index sets that the
  `:mixed` scenario never touches.

### 6.4 Normalisations planted deliberately

The inversion identifies Ã only up to a normalisation and B̃ only up to a second one, so the
planted values must be put on the same normalisation or the recovery test compares apples to
oranges. Two facts (both verified numerically in the port, not just asserted):

* Scaling all `a_j` by `k` scales `{w, q, Q, Y, vv}` by `k` and leaves `{H_M, H_R, θ, π_ij}`
  **unchanged**.
* Scaling all `b_i` by `c` scales `Φ` by `c^ε` and leaves the entire equilibrium allocation
  **unchanged**.

So the generator does: solve forward once with the raw draws → rescale `a ← a / geomean(A_eq)` and
`b ← b · (H/Φ_eq)^{1/ε}` → re-solve. The second solve reproduces the same allocation with
`geomean(Ã) = 1` and `Φ = H`, which are exactly the normalisations `calcal_adj_TD.m` /
`cmodexog.m` impose on the recovered objects.

### 6.5 Data layer / swapping in the real Berlin data

One clearly marked function is the only thing a user with the three `.mat` files has to replace:

```julia
struct ARSWData
    N     :: Int                 # number of blocks
    Q     :: Vector{Float64}     # observed floor space price, N       (MATLAB: floor06)
    HM    :: Vector{Float64}     # workplace employment, N             (MATLAB: empwpl06)
    HR    :: Vector{Float64}     # residence employment, N             (MATLAB: emprsd06)
    K     :: Vector{Float64}     # geographic area, N                  (MATLAB: area06)
    tau   :: Matrix{Float64}     # bilateral travel time, N x N        (MATLAB: tt06)
    west  :: Vector{Bool}        # former-West flag, N                 (MATLAB: dummywestr)
    coords:: Matrix{Float64}     # N x 2, for plotting only            (MATLAB: Xr/Yr)
end
```

`load_berlin_data(path)` is a documented stub that states the expected `.mat` variable names,
shapes and units, errors out with instructions, and is the single seam between synthetic and real.
Everything downstream takes an `ARSWData`.

---

## 7. Solver architecture

Everything is a hand-written damped fixed point, exactly as upstream — no `fsolve`, no `nlsolve`.
Two departures from MATLAB, both deliberate:

1. **Convergence rule.** Upstream stops on *rounding* equality (`round(x*100) == round(y*100)`,
   i.e. ~1e-2 absolute in levels, or `round(gap*1e4)==0`). That is far too loose for a 1e-8
   recovery test: the recovered fundamentals inherit the forward solve's error. Every solver
   therefore takes `rule = :tol` (default, `max |log(guess) - log(predicted)| < tol`, default
   `tol = 1e-13`) or `rule = :round` (bit-faithful to MATLAB). Both are exercised; the
   `:round` runs are reported alongside so the fidelity claim is checkable.
2. **`Q_i` vs `QT` in `smodexog.m`** — see defect 5 below. The port uses the guess `Q_i`.

Damping follows upstream: fixed 0.5 in `smodexog`/`comegaopt*`/`cmodexog`; adaptive 0.25-when-close
/ 0.5-otherwise in `smodendog`; `H = 0.05·H_up + 0.95·H` in the open city.

Julia functions are pure (no MATLAB globals): parameters travel in an `ARSWParams` struct,
`keep.m` has no analogue and is dropped.

---

## 8. Recovery-test protocol and tolerances

Reported numbers are `max_i |log(recovered_i / planted_i)|` (a scale-free relative error),
alongside the raw `maximum(abs(recovered - planted))`.

| # | Test | Target |
|---|---|---|
| R1 | Forward-solve residual (exogenous, `:mixed`) — `max |log(guess/predicted)|` over `{w,q,Q,θ}` | < 1e-12 |
| R2 | **Sequential inversion** (Alg 4→5→6→7→8) recovers `{Ã, B̃, φ̃}` | < 1e-8 |
| R3 | **Simultaneous inversion** (Alg 9→10) recovers `{Ã, B̃, φ̃}` | < 1e-8 |
| R4 | Sequential vs. simultaneous agree with each other (the cross-check `cftualprep_TD.m` builds in) | < 1e-8 |
| R5 | Same as R2/R3 on the `:specialised` scenario (exercises IcsA/IcsB) | < 1e-8 |
| R6 | **Section 7**: forward-solve with `smodendog` under {λ,δ,η,ρ} > 0, invert, then `cprod`/`cres` recover planted `{ã, b̃}` | < 1e-8 |
| R7 | Endogenous open-city solver (Alg 16) reproduces the closed-city solution when `Ū` is set to the closed-city `Ū` | < 1e-6 (H is updated at weight 0.05, so slower) |

Additional required checks:

| # | Check | Expectation |
|---|---|---|
| C1 | `Σ_j π_{j|i} = 1` for every residence row | `max |Σ - 1|` < 1e-14 |
| C2 | Labour-market clearing: `Σ_i π_{j|i} H_Ri = H_Mj` and `Σ_j π_{i|j} H_Mj = H_Ri` | < 1e-8 |
| C3 | Land-market clearing Eqs. 18 and 19 hold at the solution | < 1e-8 |
| C4 | Zero-shock counterfactual returns changes identically 0 | exactly 0.0 |
| C5 | `λ = δ = η = ρ = 0` → `smodendog` reproduces `smodexog` | < 1e-10 |
| C6 | Comparative statics: raising τ uniformly lowers CMA, lowers wages, lowers Ū | sign check |
| C7 | Scale invariance: `a → k·a` leaves `{H_M,H_R,θ,π}` unchanged and scales `{w,q,Q}` by k; `b → c·b` leaves the whole allocation unchanged | < 1e-10 |
| C8 | μ is inert: changing `mu` changes `φ` but leaves every equilibrium object unchanged | < 1e-12 |
| S1 | *(stretch)* ε recovery: generate with ε known, run Algorithms 1–3, recover ε | reported honestly either way |

For S1 the MATLAB `patternsearch` (Global Optimization Toolbox) is replaced by
`Optim.optimize(..., Brent())` over a bracketed interval plus a `NelderMead` cross-check, since the
objective is one-dimensional. The synthetic analogue of the "12 Bezirke" moment is a partition of
the 225 blocks into 9 contiguous 5×5 "districts"; the moment is the variance of log district wages.

---

## 9. Known upstream defects and how each is handled

1. **δ/η comments swapped** in `cftualprep_end_TD.m:94-96` (`delta` is commented "Density
   elasticity of residential amenity", `eta` "Productivity decay" — both describe the other
   parameter). The *code* is right. The Julia port uses the correct meanings and carries a comment
   saying so; the comment error is not propagated.
2. **`META.m` hard-codes the repo root** behind a `user` switch (`D:/Dropbox/...`). The Julia port
   is self-locating: all paths go through `@__DIR__`, and `graphs/` is created if missing.
3. **μ invisible / wrong** — see §5. Named parameter, exponent written `1 - mu`, inertness
   demonstrated numerically (check C8).
4. **`comegaoptC.m:206-214` sparse-fill loop is broken.** The "New code by GA" writes
   `cprob(Iwpl(idx), Irsd(idx)) = Ecprob(idx)`, indexing an `N×N` matrix with the *logical values*
   `0`/`1` of the mask rather than the positions it selects, and indexing the `nto×nfrom` matrix
   `Ecprob` linearly by a block counter. In MATLAB a `0` subscript is an error. The port uses the
   original, commented-out, correct line `cprob[Iwpl, Irsd] = Ecprob`.
5. **`smodexog.m:152` uses `QT`, not the guess `Q_i`,** inside the commuting probability
   (`EQQ=repmat(QT(Irsd),1,nwpl)`), although its own comment says "We assign **guesses** of
   residence floor space prices". `smodendog.m:167` and `ussmodendog.m:179` correctly use `Q_i`.
   The effect is that in `smodexog` counterfactuals residential floor prices are frozen at their
   baseline values inside the location-choice term, so the price feedback into commuting is
   switched off — and check C5 (λ=δ=η=ρ=0 ⇒ `smodendog` == `smodexog`) cannot hold. The port uses
   `Q_i` in both, and flags this as the one place where it knowingly departs from upstream
   behaviour (a `legacy_QT` keyword reproduces the upstream behaviour for comparison).
6. **`smodendog.m:110` `xtic=toc(xtic./60)`** — divides the timer handle, not the elapsed time.
   Cosmetic; not ported.
7. **`smodendog.m`/`ussmodendog.m` damping comments inverted** — `if maxLDwage < 0.1  % If we are
   far from convergence` is the *close*-to-convergence branch. Code correct, comment wrong; not
   propagated.
8. **`cmodexog.m` computes `Ephi_ij` twice** (an explicit double loop, then the vectorised version
   overwriting it). Only the vectorised version survives; the port computes it once.

*(A ninth item, D9, was found during implementation and is recorded in section 13.5b. It
turned out to be benign, which is why it is not in this pre-implementation list.)*

---

## 10. MATLAB → Julia trap checklist (applied file by file)

* `sum(M)` / `mean` / `max` / `cumsum` / `prod` are **column-wise** in MATLAB, whole-array in Julia
  → every one becomes `sum(M, dims=1)` (or `dims=2`) with an explicit `vec(...)`. This is the
  single biggest source of silent wrong answers; `Ephi=sum(sum(Ephi_ij))` becomes `sum(Ephi_ij)`
  (already the whole-array sum) while `Ephi_j=sum(Ephi_ij)'` becomes `vec(sum(Ephi_ij, dims=1))`.
* `1./x` → `1 ./ x`; `A(i,j)` → `A[i,j]`; `x(x>0)=5` → `x[x .> 0] .= 5`.
* `repmat` → `repeat`; MATLAB's implicit broadcasting → explicit `.` everywhere.
* `diag(M)` extracts in MATLAB; Julia `diagm` *builds* — `diag` is wanted.
* Variables first assigned inside a `while`/`for` body are loop-local in Julia → every quantity
  needed after a solver loop is pre-declared before it.
* MATLAB's shared script workspace has no analogue → structs and explicit arguments; `keep.m`
  dropped.
* `geomean` is not in `StatsBase` for this use → `geomean_pos(x) = exp(mean(log.(x[x .> 0])))`,
  matching MATLAB's behaviour on the masked subvectors the toolkit always passes it.
* `gamma(·)` comes from `SpecialFunctions`.
* The toolkit's `minimize` calls: grepped — the only optimiser call is `patternsearch` in
  `optimepsilon_TD86.m`; there is no function named `minimize` in this toolkit (the five call sites
  mentioned in the brief are `patternsearch`/`fmincon` references, the latter commented out).
  Replaced with `Optim.jl`.

---

## 11. Deliverables and file layout

```
AhlfeldtReddingSturmWolf-ECMA-2015/
├── TRANSLATION_PLAN.md                            (this file)
├── ahlfeldt_redding_sturm_wolf_density_model.jl   (the port; single file, `julia <file>.jl`)
├── README.md                                      (model writeup, synthetic-data caveat up front)
├── graphs/                                        (PDF output, created by the script)
├── matlab_source/  text/  *.pdf                   (staged source material, untouched)
```

Dependencies limited to `LinearAlgebra, Statistics, Printf, Random, StatsBase, SpecialFunctions,
Optim, Plots`. No `Project.toml` (repo convention).

---

## 12. Order of work

1. Skeleton: params, data structs, synthetic geography + wall + planted fundamentals.
2. Shared commuting kernel (`commuting_block`) used by every algorithm.
3. Alg 11 (`smodexog`) first — the forward solve is needed to manufacture data.
4. Algs 4–8, then 9–10. Run R2, R3, R4.
5. Algs 12–16. Run R6, R7, C5.
6. All remaining checks C1–C8.
7. Counterfactuals + figures.
8. Stretch: Algs 1–3, S1.
9. README, then update this file with §13.

---

## 13. Post-mortem — what actually happened

Run log: `julia ahlfeldt_redding_sturm_wolf_density_model.jl`, 225 blocks,
**66 of 66 checks passed**, worst recovery error 2.31e-11, elapsed 18.6 s on this machine.

### 13.1 Scope delivered

All sixteen codebook algorithms, including both stretch items (Algorithms 1-3, the epsilon
estimation). Nothing from the CORE list was dropped. The file is 2,304 lines, single file,
no `Project.toml`, runs with `julia <file>.jl`. Seven PDF figures in `graphs/`.

### 13.2 Recovery tests — planned target vs. actual

All errors are `max_i |log(recovered_i / planted_i)|`.

| # | Test | Target | **Actual** |
|---|---|---|---|
| R1 | Forward-solve residual (COLD start from uniform guesses, 474 iterations) | < 1e-12 | **9.77e-15** |
| R2 | Sequential inversion -> Ã / B̃ / φ̃ / θ | < 1e-8 | **3.33e-15 / 6.66e-16 / 1.55e-15 / 7.11e-15** |
| R3 | Simultaneous inversion -> Ã / B̃ / φ̃ | < 1e-8 | **2.66e-15 / 1.55e-15 / 2.33e-15** |
| R4 | Sequential vs. simultaneous | < 1e-8 | **1.33e-15 / 1.55e-15 / 1.55e-15** |
| R5 | `:specialised` scenario (9 commercial + 10 residential specialised blocks, 206 mixed) | < 1e-8 | **≤ 2.89e-15** |
| R6 | Section 7, {λ,δ,η,ρ} > 0: Ã, B̃, then Algs. 12/13 -> ã, b̃ | < 1e-8 | **≤ 3.11e-15** |
| R7 | Open city started at 65% of H and 20% high prices reproduces the closed city (300 iterations) | < 1e-6 | **2.31e-11** (H: 3.61e-13) |

Propositions S.3/S.4 hold numerically in this port, to machine precision, on both the
`:mixed` and `:specialised` block compositions and with and without the Section 7 externalities.

### 13.3 Checks — planned vs. actual

| # | Check | **Actual** |
|---|---|---|
| C1 | `Σ_j π_{j|i} = 1` per residence / per workplace / overall | 1.78e-15 / 1.11e-15 / 0.0 |
| C2 | Labour-market clearing both directions; `ΣH_M = ΣH_R`; **C2d (added)** the two equivalent `vv` formulas | 8.44e-15 / 9.10e-15 / 3.10e-16 / 9.44e-15 |
| C3 | Land clearing Eqs. 18 / 19 / S.31 | 1.13e-14 / 9.36e-14 / 1.33e-15 |
| C4 | Zero-shock counterfactual, all five objects | **exactly 0.0** |
| C5 | λ=δ=η=ρ=0 ⇒ Alg. 15 == Alg. 11, **under the wall shock** | ≤ 2.33e-15 |
| C6 | τ +25%: CMA −19.0%, w −0.003%, Ū −3.08%; **C6d (added)** open city H −7.53% | all signs correct |
| C7 | The two scale invariances | ≤ 2.33e-15 |
| C8 | μ inert on Ã, B̃, L, θ, solved Q; φ moves by 0.76 | exactly 0.0 |
| C9 | **(added)** `smodexog` vs `smodendog` `Q_e`/`q_e` init: both converge, same allocation | ≤ 1.58e-14 (32 vs 1 iterations) |
| S1 | epsilon: planted 6.830000 -> Brent **6.8300000051** (7.5e-10), NelderMead 6.8299999788 (3.1e-9); Alg. 1 recovers the true wages to 8.9e-16 | PASS |

Two checks had to be **rewritten after the first run because they were vacuous**, which is
worth recording since the plan did not anticipate it:

* **C5** as planned compared the two solvers seeded at the observed data. Both stop at
  iteration 1 there, so they agreed at exactly 0.0 without either solver doing anything. It
  now runs the comparison **under the wall shock**, where each solver travels 374 iterations
  to a genuinely different allocation before being compared.
* **R7** as planned seeded the open-city solver at the closed-city solution, so it converged
  in one iteration and never exercised the `H` update that is the only difference between
  Algorithms 15 and 16. It now starts at 65% of the target `H` with prices 20% off.
* **R1/R5/R6** originally reported the residual of the warm-started final solve (1-2
  iterations). They now report the **cold** solve from uniform guesses (474 / 474 / 470
  iterations).

### 13.4 Where the plan was wrong

1. **The two normalisations are not independent.** Section 6.4 of the plan asserted that
   `a -> a/geomean(A)` and `b -> b·(H/Φ)^{1/ε}` could be applied together. They cannot:
   rescaling `a` by `1/k` scales `Q` and `w` by `1/k`, and `φ_ij` carries
   `Q^{-(1-β)ε} w^ε`, so `Φ` moves by `k^{βε}`. Applying both at once left `B̃` off by
   `k^β` — a 0.161 log error, which R2b caught on the very first run. The generator now
   applies the productivity normalisation first, re-solves, then the amenity one. This is the
   single best argument for the recovery test: the bug was in the *test harness*, was
   economically subtle, and showed up immediately as a number that was not 1e-15.
2. **The wall penalty had to come down from 300 to 45 minutes**, for a reason the plan did
   not foresee. Eq. (S.44) pins transformed wages up to one global scale and the uniqueness
   proof needs the commuting matrix irreducible; a hard wall makes it nearly reducible, so the
   two halves' wage levels are tied only through flows of order `e^{-ν·penalty}` and the
   iteration count scales like `e^{+ν·penalty}`. Measured (`run_wall_conditioning`): 58/69
   iterations at no wall, 94/149 at 15 min, 218/372 at 30, 562/984 at 45, and no convergence
   within 30,000 at 120. The recovery error does not move across those rows, so it is
   conditioning, not correctness — and it is a reason ARSW estimate ε on the West-only sample.
   45 minutes is a 96% barrier and solves in about a second.
3. **Upstream's fixed damping weight of 0.5 does not converge on a large shock.** Not
   anticipated at all. On the wall shock `smodexog.m`'s constant 0.5 settles into a period-2
   limit cycle with the max log gap stuck at 0.26 for 3,000 iterations (0.63 for 20,000 at a
   120-minute penalty). A weight of 0.25 converges to 1e-13 in 374 iterations. Upstream's own
   comment anticipates the problem but only lowers the weight once the gap is already below
   0.1. `DAMP_ROBUST = 0.25` is used wherever a solve crosses a large shock, and the failure
   is reported at runtime rather than hidden.
4. **The plan's §10 claim about `minimize` was wrong**, as suspected: grepping the toolkit
   finds no function called `minimize`. The only optimiser call is `patternsearch` in
   `optimepsilon_TD86.m` (with a commented-out `fmincon` alternative). Replaced with
   `Optim.Brent()` plus a `NelderMead` cross-check.

### 13.5 Two defects found beyond the four the brief listed

* **D5 — `smodexog.m:152` uses the fixed observed `QT` where `smodendog.m:167` correctly uses
  the guess `Q_i`,** inside the commuting probability, against its own comment. This is not
  cosmetic. Under the wall shock, upstream's version differs from the corrected one by
  **0.195 in log H_M, 0.188 in log H_R, 0.157 in log Q**. With the fix, C5 holds to 2.3e-15;
  without it, C5 cannot hold at all. Decisive evidence that it is a bug, not a choice.
  `legacy_QT = true` reproduces upstream for anyone who wants it.
* **D4 — `comegaoptC.m:206-214`'s sparse fill is broken**: it indexes an `N×N` matrix with the
  logical *values* of a mask (a 0 subscript is a MATLAB error) and indexes `Ecprob` linearly
  by a block counter. The original commented-out line is correct and is what the port uses.

### 13.5b One inference that measurement overturned

Reading `smodexog.m:116-117` (`Q_e = q_e = QT`, where `smodendog.m:107-108` has
`Q_e = Q_i`), the obvious conclusion was that upstream's stopping rule -- which compares the
FULL `Q_e` and `Q_i` vectors -- could never be satisfied once any block is completely
specialised, because `Q_e[IcsA] = QT[IcsA]` is never overwritten while `Q_i[IcsA] = 0`. That
would have made `smodexog` unable to report convergence on real Berlin data, which certainly
contains blocks with zero residence or workplace employment.

The measurement (check C9) says otherwise. The damping blend
`q_i <- w*q_e + (1-w)*q_i` drags the structurally unused entries onto their own initial
values too, so the rule IS eventually satisfied: 32 iterations instead of 1 on the
`:specialised` scenario, reaching an allocation identical to 1.6e-14. The only real costs are
the extra iterations and some meaningless nonzero entries in the returned `q[IcsB]`/`Q[IcsA]`
that `Crent` never reads. Recorded as D9 and demoted from "defect" to "benign difference".
Of the defects this port flagged beyond the three named in the brief, the two substantive
ones -- D5 (the frozen `QT` in `smodexog`'s commuting term) and D9 -- were both checked
numerically rather than argued from reading the code. D5 survived and D9 did not.

### 13.6 Decisions where the source was ambiguous

1. **α is the labour share (0.80), not the floor-space share.** The upstream toolkit's
   `CLAUDE.md` parameter table says the opposite. Resolved against the paper (p. 2167: "the
   share of firm expenditure on commercial floor space (1−α) is 0.20"), against
   `cftualprep_end_TD.m`'s own comment ("Set input share of labour in production"), and against
   `comegaoptC.m:222`'s `A = (q/(1-α))^{1-α}(w/α)^α`, which inverts Eq. (12) only under the
   labour-share reading.
2. **μ = 0.75, exponent `(1-mu)` = 0.25** (paper) against the toolkit's hard-coded `K^0.75`
   (which implies μ = 0.25). Shown numerically inert (C8): Ã, B̃, L, θ and the solved Q are
   bit-identical under either convention; only the reported φ changes. `mu = 0.25` reproduces
   upstream's φ exactly.
3. **`smodexog` uses `Q_i`** — see D5 above.
4. **Algorithms 15/16 merged, 8/10 merged**, behind keywords. R7 checks the merge.
5. **`cdensity` uses the single observed floor price for both `LM` and `LR`.** Coherent only
   because a block with `H_M = 0` contributes no `LM`, a block with `H_R = 0` no `LR`, and a
   mixed block has `q = Q`. Checked rather than assumed: R5* reports `max |q - Q| = 0` on the
   incompletely specialised blocks of the `:specialised` scenario.
6. **Convergence by tolerance (1e-13/1e-14) is the default**; `rule = :round` reproduces
   upstream bit-for-bit. Both are run. Upstream's rule gives recovery errors of 4.0e-10 (Ã)
   and 1.2e-11 (B̃) sequentially, 1.3e-10 / 2.1e-10 simultaneously — five orders of magnitude
   worse than the tolerance rule, and enough to mask a real bug. It also stops at a fixed
   ABSOLUTE precision in levels, so its quality depends on the units of the data.
7. **`τ_ii = 1.5` min, not 0** — own-block commuting is not free (≈ crossing a 750 m block at
   30 km/h). Zero would distort the CMA gradient at the centre.
8. **The wall is a travel-time penalty, not a hard zero**, so `τ` stays finite and no `Inf`/`0`
   enters a log or an exponential kernel — and, per 13.4.2, so the commuting system stays
   comfortably irreducible.
9. **The CBD sits at column 8, just inside the East**, deliberately: Berlin's historic centre
   (Mitte) fell on the eastern side of the wall.
10. **The ε-estimation moment** is the variance of log wages over a 3×3 partition of the grid
    into 5×5 districts (the stand-in for `modbezirk.m`'s 12 Bezirke), computed on the
    REUNIFIED city. Upstream uses the divided-era West-only 1986 sample; on the synthetic
    divided city the ω system is ill-conditioned (13.4.2), so the reunified full city is the
    numerically well-posed choice. Stated in the function's docstring.

### 13.7 Left undone

* **The real Berlin data.** `load_berlin_data` is a documented stub: its docstring gives the
  complete `.mat` variable -> field map and the four-line `MAT.jl` implementation. Nothing else
  needs to change, because every algorithm takes an `ARSWData`.
* **`MAPIT.m`** — not ported (Mapping Toolbox + the real index-aligned shapefile). Grid
  heatmaps instead.
* **`modbezirk.m`** — not ported; a literal lookup table for real Berlin districts.
* **The paper's full GMM over {ν, ε, λ, δ, η, ρ}** (supplement §S.4) is not in the toolkit
  either — the toolkit loads the estimates from `roptimis_all_big.mat`. Only the one-step ε
  estimation (Algorithms 1-3) was in scope, and it is done and passes a recovery test.
* **`ξ_i` (the land-use regulatory wedge) and the land price `R_i`/`χ`** are not computed, as
  upstream: everything is in floor-price units.
* **The car-ban counterfactual** (`cftualexog_TD.m`'s first exercise) needs the
  public-transport-only travel-time matrix `ttpub06`, one of the three missing `.mat` files.
  Its structural analogue — re-solving with a different `τ` — is exercised three times over by
  the division, reunification and uniform-τ counterfactuals.
