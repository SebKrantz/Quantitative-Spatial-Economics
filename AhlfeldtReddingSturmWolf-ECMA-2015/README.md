# The Economics of Density — Ahlfeldt, Redding, Sturm & Wolf (2015)

Julia port of Gabriel M. Ahlfeldt's didactic MATLAB toolkit for

> Ahlfeldt, G. M., Redding, S. J., Sturm, D. M., Wolf, N. (2015).
> "The Economics of Density: Evidence from the Berlin Wall." *Econometrica* 83(6), 2127–2189.
> [doi:10.3982/ECTA10876](https://doi.org/10.3982/ECTA10876)

```bash
julia ahlfeldt_redding_sturm_wolf_density_model.jl      # ~19 s, 225 blocks, 66 checks
ARSW_PLOTS=0 julia ahlfeldt_redding_sturm_wolf_density_model.jl   # skip the figures
```

---

## ⚠️ This runs on SYNTHETIC data

**No empirical result of the paper is reproduced or claimed here.** The toolkit needs three
`.mat` files — `prepdata_big_TD.mat` (the 2006 cross-section), `prepdata_big_TD86.mat`
(1986, for the ε estimation) and `ttpublic_2006_ren.mat` (public-transport travel times) —
which are external downloads from an HU-Berlin file server and are not in this repository.
They were deliberately not fetched. So this port:

* generates a **synthetic Berlin-like block geography** with a wall, plants known
  fundamentals, and solves the model forward to manufacture "observed" data;
* uses the paper's estimated parameters (λ = 0.0710, ε = 6.83 / 6.694, ν = 0.07 / 0.0987,
  δ, η, ρ) as **inputs**, never as targets;
* reports gradients, elasticities and counterfactual magnitudes that are properties of the
  synthetic draw, not of Berlin.

**To run it on the real data**, implement the single function `load_berlin_data` (its
docstring gives the complete variable → field map; it is four lines with
[MAT.jl](https://github.com/JuliaIO/MAT.jl)). Nothing else needs to change: every algorithm
in the file takes an `ARSWData` and does not care where it came from.

### Why synthetic data makes the test *stronger*, not weaker

Supplement §S.3.1.6 (Proposition S.3) proves a **one-to-one mapping**

```
{α, β, μ, ε, κ}  +  observed {Q, H_M, H_R, K, τ}   →   {Ã, B̃, φ̃}
```

and §S.3.2 (Proposition S.4) extends it with {λ, δ, η, ρ} to the fundamentals {ã, b̃}. So the
primary validation here is a **recovery test**: plant fundamentals → solve forward → invert →
check they come back. If the mapping theorem holds and the port is right, this *must* work to
machine precision. Matching Berlin's numbers would be a weaker test, because a data problem
and an implementation bug are indistinguishable in it. Here they are not.

**Result: every recovery test lands at 10⁻¹⁵–10⁻¹¹, three to seven orders of magnitude inside
the stated tolerances.** See [Validation](#validation) below.

---

## 1. The model

`N` blocks; each block is simultaneously a candidate **workplace** `j` and a candidate
**residence** `i`. Workers are homogeneous but draw i.i.d. Fréchet(ε) shocks over
(residence, workplace) pairs, so commuting is a discrete choice over `N²` bilateral options.

### 1.1 Commuting and residence choice

| | |
|---|---|
| (4) | `φ_ij = e^{-εκτ_ij} · B̃_i^ε · Q_i^{-(1-β)ε} · w̃_j^ε`, `π_ij = φ_ij / Φ`, `Φ = Σ_ij φ_ij` |
| (5) | `H_Ri = (Σ_j φ_ij / Φ)·H`, `H_Mj = (Σ_i φ_ij / Φ)·H` |
| (6) | `π_{j\|i} = φ_ij / Σ_{j'} φ_{ij'}` |
| (S.20) | `E[w\|i] = Σ_j π_{j\|i} w̃_j`; total worker income `vv_i = E[w\|i]·H_Ri` |
| (9) | `Ū = γ·Φ^{1/ε}`, `γ = Γ((ε-1)/ε)` |

The key structural fact the whole quantification rests on: **the conditional probability
`π_{j|i}` contains neither `B̃_i` nor `Q_i`** — they cancel out of the conditional. So

```
(S.44)   H_Mj = Σ_i [ ω_j e^{-ν τ_ij} / Σ_s ω_s e^{-ν τ_is} ] · H_Ri ,   ω_j = w̃_j^ε, ν = κε
```

determines transformed wages from `{H_M, H_R, τ}` **alone**, with no amenity, floor-price or
productivity information, and — because only `ν` appears — independently of `ε`. Lemmas S.6–S.7
show this system satisfies gross substitution, so the fixed point is unique up to one global
scale (normalised to `geomean(w̃) = 1`).

### 1.2 Production and the floor-space market

| | |
|---|---|
| (10) | `Y_j = A_j · H_Mj^α · (θ_j L_j)^{1-α}` |
| (12) | `q_j = (1-α)(α/w_j)^{α/(1-α)} A_j^{1/(1-α)}` (FOC + zero profit) |
| (13) | no arbitrage: `θ_j = 1` if `q > ξQ`, `∈[0,1]` if `q = ξQ`, `0` if `q < ξQ` |
| (18)/(S.30) | `θ_j L_j = ((1-α)A_j/q_j)^{1/α}·H_Mj` — commercial land clearing |
| (19)/(S.29) | `(1-θ_i)L_i = (1-β)·E[w\|i]H_Ri / Q_i` — residential land clearing |
| (15)/(S.31) | `L_i = φ_i·K_i^{1-μ}` — floor-space supply / density of development |

`θ_i` is the **commercial** floor-space share. Blocks split into three regimes, handled by the
index sets `IcsA` (productive but uninhabited, θ=1), `IcsB` (inhabited but unproductive, θ=0)
and `Iis` (incompletely specialised, `q_i = Q_i` and θ interior).

### 1.3 Inversion

| | |
|---|---|
| (S.46) | `W_i = Σ_s e^{-ντ_is} w̃_s^ε` — residential **commuting market access** (CMA) |
| (S.47) | `B̃_i/B̄ = (H_Ri/H̄_R)^{1/ε}·(Q_i/Q̄)^{1-β}·(W_i/W̄)^{-1/ε}` |
| (S.48) | `Ã_j = (q_j/(1-α))^{1-α}·(w̃_j/α)^α` |

Read S.47 the way the supplement does (p. 30): high residence employment together with high
floor prices must be explained *either* by good commuting market access *or* by attractive
residential amenities. Whatever the data leaves over is `B̃`.

### 1.4 Endogenous agglomeration (Section 7)

```
(20)  A_j = a_j·Υ_j^λ ,   Υ_j = Σ_s e^{-δ τ_js}(H_Ms/K_s)      production externality
(21)  B_i = b_i·Ω_i^η ,   Ω_i = Σ_s e^{-ρ τ_is}(H_Rs/K_s)      residential externality
```

Note the **spatial decay kernel**: the externality is a travel-time-weighted sum of employment
density over *all other* blocks, not a power function of own density. That is richer than the
agglomeration/amenity externalities in this repository's other models (see
[§6](#6-relation-to-the-other-models-in-this-repository)).

---

## 2. What is implemented

All **16 codebook algorithms**, including the two "stretch" ones. `matlab_source/` holds the
flattened, section-prefixed copies of the upstream `.m` files.

| Alg. | Equations | MATLAB | Julia |
|---|---|---|---|
| 1 | S.44 | `comegaoptO.m` | `solve_transformed_wages` |
| 2 | S.64 | `cdensityoptren.m` | `epsilon_objective` |
| 3 | — | `optimepsilon_TD86.m` | `estimate_epsilon` (Optim.jl) |
| 4 | S.44, S.48 | `comegaoptC.m` | `solve_wages_and_productivity` |
| 5 | S.46, S.47 | `camen.m` | `solve_amenities` |
| 6 | (12), Φ=H | `calcal_adj_TD.m` | `rescale_fundamentals` |
| 7 | S.20 | `expincome.m` | `expected_income` |
| 8 | S.29–S.31 | `cdensity.m` | `solve_density` |
| 9 | (4),(5),(12) | `cmodexog.m` | `invert_simultaneous` |
| 10 | S.29–S.31 | `cdensityE.m` | `solve_density` |
| 11 | (4),(5),(9),(10),(12),(18),(19) | `smodexog.m` | `solve_equilibrium_exog` |
| 12 | (20)/S.55 | `cprod.m` | `decompose_productivity` |
| 13 | (21)/S.56 | `cres.m` | `decompose_amenity` |
| 14 | (9) | `ubar.m` | `reservation_utility` |
| 15 | 11 + (20),(21) | `smodendog.m` | `solve_equilibrium_endog(closed=true)` |
| 16 | 15, Ū fixed | `ussmodendog.m` | `solve_equilibrium_endog(closed=false)` |

Algorithms 15/16 and 8/10 are merged into one function apiece with a keyword: upstream's
`ussmodendog.m` is `smodendog.m` plus a single extra line updating total employment, and
`cdensityE.m` is `cdensity.m` returning two more vectors. Keeping four near-identical 400-line
files was the wrong thing to carry across.

**Not ported, deliberately:**

* `prepdata_TD.m`, `prepdata_TD86.m`, `prepdata_TD06ttpub.m` — they only rebuild the missing
  `.mat` files. Read for schema, then stopped.
* `MAPIT.m` — Jenks-classed choropleths on the Berlin shapefile; needs the Mapping Toolbox and
  the real, index-aligned geography. The port draws grid heatmaps instead, which is the honest
  synthetic analogue.
* `modbezirk.m` — a literal 23→12 historic-Bezirke lookup table; meaningless synthetically. Its
  role in the ε estimation is played by a 3×3 partition of the grid into 5×5 districts.
* `keep.m` — a utility for pruning MATLAB's shared script workspace. No Julia analogue.
* `ξ_i` (the land-use regulatory wedge) and the land price `R_i`/`χ` are not computed — as
  upstream, everything is expressed in floor-price units.

---

## 3. The synthetic data

### Geography

15 × 15 = 225 blocks on a 1 km lattice, CBD at the centre, block land area smaller downtown
(finer subdivision) with lognormal noise, travel times in **minutes** — the unit in which
ν = 0.07 is calibrated — at 30 km/h, with `τ_ii = 1.5` min so own-block commuting is not free.

### The wall

A vertical cut between columns 7 and 8; straddling pairs get `τ_ij += 45` minutes, which
multiplies `φ_ij` by `e^{-3.15} ≈ 0.043` — cross-wall commuting falls by 96%. With the
defaults the **CBD lands in column 8, just inside the East** — deliberately, because that is
Berlin's situation: the historic centre (Mitte) fell on the eastern side, which is why
division hit West Berlin's access to the old core hardest.

Running the same fundamentals under the two travel-time matrices gives the paper's identifying
variation a synthetic analogue.

### Planted fundamentals

`log a_j`, `log b_i` rising with distance from the CBD; `log φ_i` falling; independent
lognormal noise on each, so `a` and `b` are not mechanically collinear. Seeded
(`MersenneTwister(20150615)`) so the whole run is reproducible. Two scenarios:

* `:mixed` — every block incompletely specialised, so `q_i = Q_i` and the single observed floor
  price is unambiguous (supplement §S.2.6's benchmark case);
* `:specialised` — a central 3×3 office district with `b_i = 0` (θ=1, `H_R = 0`) and two
  park/housing patches with `a_j = 0` (θ=0, `H_M = 0`), which exercises the `IcsA`/`IcsB` index
  sets that `:mixed` never touches.

### Normalisations

The inversion identifies `Ã` only up to `geomean(Ã) = 1` and `B̃` only up to `Φ = H`, so the
planted values are put on the same footing before comparison. Two exact invariances make that
legitimate, and the port **verifies both numerically** (check C7) rather than asserting them:

* `a → k·a` scales `{w, q, Q, Y, vv}` by `k` and leaves `{H_M, H_R, θ, π_ij}` unchanged;
* `b → c·b` scales `Φ` by `c^ε` and leaves the **entire** allocation unchanged.

They are **not independent**, which is the one subtlety worth flagging: rescaling `a` by `1/k`
scales `Q` and `w`, and `φ_ij` carries `Q^{-(1-β)ε}w^ε`, so `Φ` moves by `k^{βε}`. The amenity
normalisation therefore has to be applied *after* the productivity one. Getting the order wrong
leaves `B̃` off by a factor of `k^β` — a 16% log error, which is exactly what the recovery test
caught on the first run.

---

## 4. Validation

`66 of 66 checks passed`, worst recovery error `2.3e-11`, 19 seconds. All errors are
`max_i |log(recovered_i / planted_i)|`, a scale-free relative error.

### Recovery tests

| # | Test | Target | **Actual** |
|---|---|---|---|
| R1 | Forward-solve residual (cold start, uniform guesses, 474 iterations) | 1e-12 | **9.77e-15** |
| R2 | Sequential inversion (Alg. 4→5→6→7→8) → Ã / B̃ / φ̃ / θ | 1e-8 | **3.33e-15 / 6.66e-16 / 1.55e-15 / 7.11e-15** |
| R3 | Simultaneous inversion (Alg. 9→10) → Ã / B̃ / φ̃ | 1e-8 | **2.66e-15 / 1.55e-15 / 2.33e-15** |
| R4 | Sequential vs. simultaneous (the cross-check `cftualprep_TD.m` builds in) | 1e-8 | **1.33e-15 / 1.55e-15 / 1.55e-15** |
| R5 | `:specialised` scenario, both routes | 1e-8 | **≤ 2.89e-15** |
| R6 | Section 7 (λ,δ,η,ρ > 0): Ã, B̃, then `cprod`/`cres` → ã, b̃ | 1e-8 | **≤ 3.11e-15** |
| R7 | Open city (Alg. 16) started at 65% of H and 20% high prices reproduces the closed city | 1e-6 | **2.31e-11** (H: 3.61e-13) |

### Model checks

| # | Check | **Actual** |
|---|---|---|
| C1 | `Σ_j π_{j\|i} = 1` per residence row / per workplace / overall | 1.78e-15 / 1.11e-15 / 0.0 |
| C2 | Labour-market clearing both directions; `ΣH_M = ΣH_R`; the two equivalent `vv` formulas | 8.44e-15 / 9.10e-15 / 3.10e-16 / 9.44e-15 |
| C3 | Land-market clearing Eqs. 18 / 19 / S.31 | 1.13e-14 / 9.36e-14 / 1.33e-15 |
| C4 | Zero-shock counterfactual (H_M, H_R, Q, θ, Ū) | **exactly 0.0** on all five |
| C5 | λ=δ=η=ρ=0 ⇒ `smodendog` reproduces `smodexog`, **under the wall shock** (374 iterations each) | ≤ 2.33e-15 |
| C6 | τ +25% ⇒ CMA −19.0%, wage −0.003%, Ū −3.08%; open city ⇒ H −7.53% | all signs correct |
| C7 | The two scale invariances (see §3) | ≤ 2.33e-15 |
| C8 | μ is inert: Ã, B̃, L, θ and the solved Q all unchanged; only φ moves | exactly 0.0 (φ moves by 0.76) |
| C9 | `smodexog`'s vs `smodendog`'s `Q_e`/`q_e` initialisation: both converge, same allocation | ≤ 1.58e-14 (32 vs 1 iterations) |
| S1 | **ε recovered**: planted 6.830000 → Brent **6.8300000051** (rel. err 7.5e-10), NelderMead 6.8299999788 | PASS |

**On C5 and C6.** Both would be vacuous if seeded at the answer — with the observed data as the
starting guess every solver stops at iteration 1 and any two of them "agree". So C5 is run under
a real shock (the wall going up, 374 iterations for each solver) and C6's wage prediction is
followed up with the open-city version, where the population margin the closed city cannot use
makes the response first-order (H falls 7.5%).

**On S1 (ε recovery).** This is a genuine test, not a tautology: Algorithm 1 sees only
employment and travel times — **no wage data at all** — and recovers `ω` to 8.9e-16 of the true
wages up to scale; the single moment (the variance of log district wages) then pins ε. Upstream
uses MATLAB's `patternsearch`; here the objective is one-dimensional and smooth, so
`Optim.Brent()` over ε ∈ [2, 24] is used, cross-checked against `NelderMead` from upstream's
starting value of 4.

### Figures (`graphs/`, PDF)

`recovery_fundamentals` (planted vs. inverted, side by side on the grid) ·
`recovery_scatter` (the 45° line) · `convergence_smodexog` · `wall_division` ·
`wall_gradient_break` · `section7_decomposition` · `commuting_market_access`

---

## 5. Counterfactuals

Seeded from the simultaneous quantification, as upstream's `cftualexog_TD.m` and
`cftualendog_*_TD.m` do.

**1. Division** — the wall goes up, same fundamentals. Ū falls 6.15%. Blocks within 1.6 km of
the wall lose 5.78% of residents and 4.19% of floor-price value, against a city-wide floor-price
change of −0.18%: the gradient breaks at the wall, which is the paper's signature.

**2. Reunification** — quantify the model on **divided** data (inversion recovers the planted
fundamentals to 2.3e-14 even there), then take the wall away. Ū rises 6.55%; residence
employment near the former wall rises 7.85% against −1.89% elsewhere, floor prices +4.54%
against −1.13%. This is the direction the paper actually exploits, and the sign pattern matches
its reduced form.

**3. Eastern renewal** (+10% productivity in the East, upstream's second counterfactual):

| | Ū | H_M East | H_M West | A East |
|---|---|---|---|---|
| exogenous fundamentals (Alg. 11) | +4.43% | +21.96% | −26.31% | +10% by construction |
| endogenous agglomeration (Alg. 15) | +4.51% | +24.46% | −29.94% | **+10.90%** |

The shock is +10% to the *fundamental* ã; the gap between +10% and +10.90% is the production
externality of Eq. (20) amplifying it, and the larger employment reallocation is that
amplification feeding back into location choice.

---

## 6. Relation to the other models in this repository

| Mechanism | ARSW2015 (here) | Elsewhere in this repo |
|---|---|---|
| Bilateral **commuting** (residence × workplace Fréchet), `H_M ≠ H_R` | ✔ core | **absent everywhere** |
| Explicit **floor-space market**, endogenous commercial/residential split θ | ✔ core | absent (no model here has a non-tradable housing good) |
| Endogenous **transport cost** (route choice, congestion) | ✘ — `τ_ij` is exogenous data | AllenArkolakis-RES-2022, FuchsFoongWong-MMN-2026, Santamaria-2022, Optimal Transport Networks |
| Multi-region CES **goods trade** | ✘ — a single freely-traded numeraire | all of them |
| Agglomeration externality | spatially decayed, `Σ_s e^{-δτ_js}(H_Ms/K_s)` | QSE-ARE-2017 / AA2022: a local power of own density |

* **QSE-ARE-2017** (`Helpman_E.jl`) is the same family — CES consumption, Cobb-Douglas floor
  space, agglomeration and amenity externalities — but has no bilateral commuting choice: one
  location is a single workplace-and-residence unit, so no `H_M ≠ H_R`, no `π_ij`, no ε.
* **Redding-JIE-2016** is the closest *code* analogue: `thetaepsopt.jl` ↔ `optimepsilon_TD86.m`,
  `solveab.jl`/`solveHab.jl` ↔ `cprod.m`/`cres.m`, and the same sequential-inversion shape.
  **But its `π_ij` is a goods-trade share, not a commuting probability.** The code looks nearly
  isomorphic; the economics is not.
* **AllenArkolakis-RES-2022** endogenises what ARSW takes as given (transport cost via route
  choice and congestion) and takes as given what ARSW endogenises (it has no floor-space market
  and no separate workplace/residence margin). Complements, not nests.

### ⚠️ Notation collisions — the λ row is genuinely dangerous

| Symbol | **ARSW2015 (this folder)** | Elsewhere in this repo |
|---|---|---|
| **λ** | **Agglomeration elasticity of productivity, 0.0710. Higher λ RAISES productivity.** | **AllenArkolakis-RES-2022: traffic CONGESTION elasticity, 0.07–0.09. Higher λ RAISES transport cost.** Same letter, near-identical magnitude, **opposite** economic content. A value transplanted between them looks plausible and is wrong. |
| ε | Fréchet shape for **commuting**, 6.83 / 6.694 | migration elasticity, 3–5 |
| α | **Labour** share in production, 0.80 (floor space in firm costs = 1−α = 0.20) | AA2022: productivity externality, 0.1 |
| β | Expenditure share on the tradable good, 0.75 (floor space = 1−β = 0.25) | AA2022: amenity externality, −0.3 |
| θ | **Commercial floor-space share** of a block, ∈ [0,1] | trade / Fréchet elasticity, 6–10 |
| π_ij | **Commuting probability** (live i, work j) | Redding-JIE-2016: bilateral **trade** share |

---

## 7. Upstream defects, and what was done about each

**D1 — δ/η comments swapped.** `cftualprep_end_TD.m:94-96` comments `delta` as "Density
elasticity of residential amenity" and `eta` as "Productivity decay"; each describes the other's
parameter. The *code* and the codebook are correct. Not propagated.

**D2 — hard-coded repository root.** `META.m:25-34` `cd`s into an author-specific absolute path
behind a `user` integer switch. The Julia port is self-locating through `@__DIR__`, and creates
`graphs/` if it is missing.

**D3 — μ is invisible, and probably wrong.** μ is never a named variable upstream; it appears
only as the literal exponent `0.75` in four files (`cdensity.m:57`, `cdensityE.m:61`,
`smodendog.m:85`, `ussmodendog.m:95`). Read as `K^{1-μ}` that implies μ = 0.25, but the paper
(p. 2167) sets "the share of land in construction costs (1−μ) equal to 0.25", i.e. μ = 0.75 and
the exponent should be **0.25**. Here `mu` is a named parameter and the exponent is written
`(1 - mu)`. The disagreement is **numerically inert** — `cdensity` divides by `K^{1-μ}` and the
solvers immediately multiply by it again — and the port demonstrates that (check C8: Ã, B̃, L, θ
and the solved Q are all bit-identical under μ = 0.25 vs 0.75; only the reported φ moves).

**D4 — `comegaoptC.m`'s sparse fill is broken.** The "New code by GA" block (lines 206–214)
writes `cprob(Iwpl(idx), Irsd(idx)) = Ecprob(idx)`, indexing an `N×N` matrix with the *logical
values* `0`/`1` of the mask rather than the positions it selects, and indexing the `nto×nfrom`
matrix `Ecprob` linearly by a block counter. A `0` subscript is a MATLAB error. The port uses
the original, correct, commented-out line `cprob[Iwpl, Irsd] = Ecprob`.

**D5 — `smodexog.m:152` freezes residential floor prices inside the commuting term.** It builds
`EQQ=repmat(QT(Irsd),1,nwpl)` from the **fixed observed** `QT`, although its own comment says
"We assign *guesses* of residence floor space prices"; `smodendog.m:167` and
`ussmodendog.m:179` correctly use the current guess `Q_i`. The consequence is that in a
`smodexog` counterfactual, residential price changes do not feed back into location choice.

This one is not cosmetic, and the port measures it. Under the wall shock:

| | H_M | H_R | Q | Ū |
|---|---|---|---|---|
| upstream `QT` vs. corrected `Q_i` (max log diff) | **0.195** | **0.188** | **0.157** | 8.3e-4 |

So a 19% difference in workplace employment. With the fix, check C5 (λ=δ=η=ρ=0 ⇒ Alg. 15
reproduces Alg. 11) holds to 2.3e-15; without it, it cannot hold at all. That is decisive: it is
a bug, not a modelling choice. `legacy_QT = true` restores the upstream behaviour for anyone who
wants to reproduce it.

**D9 — `smodexog.m:116-117` initialises `Q_e = q_e = QT`** (the full observed vector) where
`smodendog.m:107-108` uses `Q_e = Q_i`, `q_e = q_i`. Only `Q_e[IcsB]`, `Q_e[Iis]`,
`q_e[IcsA]` and `q_e[Iis]` are written inside the loop, so upstream leaves
`Q_e[IcsA] = QT[IcsA]` against `Q_i[IcsA] = 0` while its stopping rule compares the *full*
vectors — which looks like it must fail forever on any completely specialised block.

**It does not, and this is worth recording as a correction.** The damping blend
`q_i ← w·q_e + (1-w)·q_i` drags the structurally unused entries onto their own initial values
as well, so the rule is eventually satisfied: on the `:specialised` scenario upstream's
initialisation converges in 32 iterations against 1, and reaches an allocation identical to
1.6e-14 (check C9). The only real costs are those extra iterations and some meaningless
nonzero values in the returned `q[IcsB]`/`Q[IcsA]`, which `Crent` never reads. The first
reading of this code was that it broke convergence outright; measuring it said otherwise.
`legacy_Qe_init = true` reproduces upstream.

**D6 — `smodendog.m:110`** `xtic=toc(xtic./60)` divides the timer handle rather than the elapsed
time. Cosmetic; not ported.

**D7 — damping comments inverted.** `smodendog.m:293-304`'s `if maxLDwage < 0.1 % If we are far
from convergence` is the *close*-to-convergence branch. Code right, comment wrong; not
propagated.

**D8 — `cmodexog.m` computes `Ephi_ij` twice**, once with an explicit double loop and then again
vectorised, overwriting it. The port computes it once.

---

## 8. Two findings about the solvers

Both came out of the validation, and both are reported by the script at runtime.

### Upstream's fixed damping weight of 0.5 does not converge on a large shock

`smodexog.m:286` fixes the weight on the predicted value at 0.5 (its `if maxLD < 0.5` branch
assigns 0.5 in both arms, so it is a constant). On the wall shock that iteration settles into a
period-2 limit cycle with the max log gap **stuck at 0.26 for 3,000 iterations**. A weight of
0.25 converges to 1e-13 in **374** iterations. Upstream's own comment anticipates the problem —
"it can be beneficial to choose a smaller weight if the algorithm is bouncing too much" — it
just never lowers the weight while the gap is still large (`smodendog.m` only drops to 0.25
once the gap is already below 0.1). The port exposes `damp` and uses 0.25 wherever a solve has
to cross a large shock.

### A wall degrades the *conditioning* of the inversion, not its correctness

Eq. (S.44) pins transformed wages only up to one global scale, and the uniqueness proof
(Lemmas S.6–S.7, gross substitution) needs the commuting matrix to be irreducible. A hard wall
makes it nearly reducible: the two halves' wage *levels* are then tied together only through
flows of order `e^{-ν·penalty}`, so the slow mode of the damped iteration has modulus
`1 − O(e^{-ν·penalty})` and the iteration count scales like `e^{+ν·penalty}`. Measured:

| wall penalty (min) | `e^{-ν·penalty}` | Alg. 4 iterations | Alg. 9 iterations | Alg. 9 Ã error |
|---|---|---|---|---|
| 0 | 1.000 | 58 | 69 | 2.7e-15 |
| 15 | 0.350 | 94 | 149 | 3.1e-15 |
| 30 | 0.123 | 218 | 372 | 6.9e-15 |
| 45 | 0.043 | 562 | 984 | 1.9e-14 |
| 120 | 2.2e-4 | not converged in 30,000 (gap 2.5e-8) | not converged in 30,000 (Ã error 2.5e-4) | — |

The iteration count tracks `1/e^{-ν·penalty}` almost exactly while the recovery error does not
move — so this is conditioning, not correctness. It is also, incidentally, a reason ARSW estimate
ε on the **West-only** 1986 sample rather than on the divided city as a whole. The default wall
penalty here is 45 minutes, which is a 96% barrier and still solves in about a second.

### Stopping rule

Upstream stops on *rounding* equality (`round(x·100) == round(y·100)` for the equilibrium
solvers, `round(gap·1e4) == 0` for the inversions) rather than on a tolerance. That is fine for
mapping and for the paper's purposes, but it stops at a fixed **absolute** precision in levels,
so how good it is depends on the units of the data, and it is not sharp enough to distinguish a
correct inversion from a subtly wrong one:

| rule | sequential Ã / B̃ | simultaneous Ã / B̃ |
|---|---|---|
| upstream `round(gap·1e4)==0` | 4.0e-10 / 1.2e-11 (39 it) | 1.3e-10 / 2.1e-10 (43 it) |
| tolerance 1e-14 | 3.3e-15 / 6.7e-16 (58 it) | 2.7e-15 / 1.6e-15 (69 it) |

Every solver takes `rule = :tol` (default) or `rule = :round` (bit-faithful to upstream); both
are run and reported.

---

## 9. Contents of this folder

```
ahlfeldt_redding_sturm_wolf_density_model.jl   the port (2,304 lines, single file)
TRANSLATION_PLAN.md                            algorithm map, synthetic design, test protocol,
                                               and a post-mortem of what actually happened
README.md                                      this file
graphs/                                        7 PDF figures, written by the script
matlab_source/                                 all 29 upstream .m files, section-prefixed
text/                                          pdftotext -layout of the paper, supplement, codebook
ARSW2015-paper.pdf  ARSW2015-supplement.pdf  ARSW2015-codebook.pdf
```

Dependencies (all already in this repository's global Julia environment, and no `Project.toml`,
per repo convention): `LinearAlgebra`, `Statistics`, `Printf`, `Random`, `SpecialFunctions`,
`Optim`, `Plots`. Julia 1.12.4.

The upstream checkout lives at `../Ahlfeldt/ARSW2015-toolkit/` (read-only; its `CLAUDE.md` maps
all 16 codebook algorithms to files with equation references and was the starting point for this
port).

---

## 10. Key references

* Ahlfeldt, Redding, Sturm & Wolf (2015), "The Economics of Density: Evidence from the Berlin
  Wall", *Econometrica* 83(6), 2127–2189 — and its **Supplemental Material**, which is where the
  model derivations, the identification propositions S.3/S.4 and the computational appendix live.
* Ahlfeldt, G. M., *ARSW2015 toolkit codebook* — 16 numbered algorithms mapping onto the MATLAB
  files, the map this port follows.
* Davis & Ortalo-Magné (2011) — source of 1−β = 0.25.
* Valentinyi & Herrendorf (2008) — source of 1−α = 0.20.
* Combes, Duranton & Gobillon (2014); Epple, Gordon & Sieg (2010) — source of 1−μ = 0.25.
