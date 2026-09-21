# Implementation Plan: Brockhaus, Hinz & Serfaty (2026) — Chokepoint Disruptions with Endogenous Freight Costs and Monopoly Tolls

> **Status: implemented.** See [`brockhaus_hinz_serfaty_chokepoint_model.jl`](./brockhaus_hinz_serfaty_chokepoint_model.jl) and README §8 for results. All 12 steps of §12 are done and all verification items of §8 pass. Five things came out differently than planned:
>
> 1. **The expenditure system had to be non-homogeneous.** §5 planned to solve the levels system as written and read wages off it. That system is scale-free, so relative country sizes are fixed by $(\pi, \gamma, \alpha)$ alone and the implied wage map has derivative of order $\theta$ — it diverges. The solver instead takes labour income $\hat w_i VA_i$ as given (making the Leontief map a strict contraction) and updates wages from an excess-demand function, with an adaptive step. `solve_levels` survives, but only for building the baseline, where the eigenvector *is* what determines country sizes.
> 2. **The transport block needs adaptive damping, not a fixed step.** The congestion-relief loop gain $\eta^R\chi_s\lambda_s(1-\rho_r)$ exceeds one for high-$\chi$ sectors at the calibrated elasticities, so a fixed step oscillates and produces a spurious branch on which large fees look profitable. This surfaced as a non-monotone calibration curve.
> 3. **The fee search needs a grid scan before golden section.** With a large $\eta^R$ the revenue function is a narrow spike just above $\phi = 1$ on a long flat floor; pure golden section wanders on the flat part.
> 4. **The $\chi = 0$ limiting check in §8.17 was wrong.** A route *closure*'s cost is $\chi$-independent — dropping a route from the CES raises the composite by $\rho_{open}^{-1/\eta^R}$ regardless — because the baseline cost gaps are inferred from observed shares. The implemented check uses a *fee* instead, which does vanish at $\chi = 0$.
> 5. **$\eta^R$ calibrates to 140, not 350**, and $\chi_{energy}$ was set to 0.10 rather than 0.06. Both are consequences of identifying off this baseline rather than GTAP; 0.10 is also closer to the ≈0.12 implied by the paper's own Hormuz revenue numbers.

## Context

Implement the quantitative block of Brockhaus, Hinz & Serfaty (2026), *"Navigating Shocks: The Ripple Effects of Shipping Route Closures"* (Banque de France WP 1057), as a self-contained Julia script. See [`README.md`](./README.md) for the full model summary and equation references.

The model is a multi-country, multi-sector Caliendo–Parro economy with input–output linkages in which the bilateral trade wedge is **endogenous** through a nested mode-over-route transport block featuring:

1. **Route congestion** on individual passages, $\lambda_s$
2. A **global shipping capacity constraint**, $\gamma^G_s$
3. **Monopoly toll-setting** by canal authorities, $\phi_q$, whose revenue is rebated to the collector country

**Target file:** `BrockhausHinzSerfaty-2026/brockhaus_hinz_serfaty_chokepoint_model.jl`

**Reference implementations for conventions:** `FuchsFoongWong-MMN-2026/fuchs_foong_wong_multimodal_model.jl` (nested mode-over-route CES, congestion, exact-hat solver), `AllenArkolakis-RES-2022/allen_arkolakis_traffic_model.jl` (hat-algebra counterfactual structure), `QRE-HoRaUE-2025/workhorse_qre_framework.jl` (nested fixed-point style).

---

## 0. The Data Problem and the Two-Track Answer

The paper's baseline is **GTAP 11** (160 countries × 65 sectors, ref. year 2017), plus proprietary AIS trajectories and Panjiva customs micro-data. None of that is in this repository, and GTAP is licensed.

Following the repo's existing pattern (`fuchs_foong_wong_multimodal_model.jl` on a synthetic US-inspired network, `ffw_thailand_calibration.jl` on real exported data), the implementation is split:

- **Track A (this plan, the deliverable):** one self-contained script with a **stylized 24-country × 8-sector calibration** that carries the paper's real geography of chokepoints, its calibrated elasticities, and its actual Gulf-exposure shares. It reproduces the *mechanisms* and the *qualitative incidence pattern* — not the third decimal of Table 8.
- **Track B (documented hooks, deferred):** a `load_gtap()` / `load_ais_route_shares()` interface plus a `export_searoutes.py` helper, so a user with GTAP 11 and AIS access can swap in the real baseline without touching the solver.

Everything below the data-assembly layer is **identical across tracks** — the solver never sees where the numbers came from.

---

## 1. Model Dimensions and Data Structures

```julia
N  = 24    # countries
S  = 8     # sectors
M  = 3     # modes: 1 = sea, 2 = air, 3 = other (land/pipeline/services)
R  = 8     # sea routes (see §2.2); air and other have a single route each
```

### Baseline economy (a `Baseline` struct / NamedTuple)

| Field | Dims | Meaning |
|---|---|---|
| `π[i,j,s]` | N×N×S | share of $j$'s sector-$s$ expenditure sourced from $i$ (sums to 1 over `i`) |
| `X[j,s]` | N×S | total sector-$s$ expenditure (intermediate + final) in $j$ |
| `Y[i,s]` | N×S | gross output of $(i,s)$, **net of toll take** |
| `γ_VA[i,s]` | N×S | value-added share of gross output |
| `γ_IO[i,k,s]` | N×S×S | share of input $k$ in gross output of $(i,s)$; `γ_VA + sum(γ_IO, dims=k) = 1` |
| `α_fin[j,s]` | N×S | final-demand shares, sums to 1 over `s` |
| `I[j]`, `VA[j]`, `D[j]` | N | income, value added, fixed nominal deficit |
| `μ[i,j,s,m]` | N×N×S×M | mode shares |
| `ρ[i,j,s,m,r]` | N×N×S×M×R | route shares within mode |
| `δ[i,j,s,m,r]` | N×N×S×M×R | ton-km capacity weights |
| `Ξ̄[m,r]`, `Ψ̄m[m]`, `Ψ̄` | — | baseline congestion / capacity aggregates |
| `φ_base[q]` | n_auth | baseline toll **levels** by authority |
| `Π[j]` | N | baseline toll revenue accruing to $j$ |

### Parameters

```julia
θ[s]      # trade elasticity (Fontagné et al. 2022, σ_s − 1)  — PLACEHOLDERS, see §3
ηR        # route elasticity μ^d, CALIBRATED (paper: 350)
ηM[s]     # mode elasticity (Ko et al. 2025, Table 4; pooled fallback 2.44)
λ[s]      # route congestion (0.135 = 0.38 × 0.35, uniform)
γG[s]     # global capacity (0.44, uniform)
γM[s]     # mode capacity (0.0 — switched off in the paper's quantification)
χ[s]      # freight cost share (0.01 electronics/pharma … 0.30 bulk; trade-wtd ≈ 0.09)
```

Memory sanity check: `ρ` at (24, 24, 8, 3, 8) is 110k `Float64` — trivial. At the paper's (160, 160, 65, 3, 8) it is 400M entries ≈ 3.2 GB, so **Track B must exploit that route shares are common across sectors within a pair** and store `ρ[i,j,m,r]` (2.5 MB), reconstructing sectoral variation on the fly. Note this in the code.

---

## 2. Geography: Countries, Passages, Authorities

### 2.1 Countries

Chosen to span every incidence group the paper reports:

| Group | Countries |
|---|---|
| Canal collectors | Egypt, Panama |
| Hormuz collectors | Iran ($g=0.95$), Oman ($g=0$) |
| Gulf exposed | Qatar, Kuwait, Iraq, Bahrain ($g=1$); Saudi Arabia ($g=0.80$); UAE ($g=0.75$) |
| Europe | Germany, France, Netherlands, Italy, Malta |
| Asia / Americas | China, Japan, Korea, India, Singapore, United States |
| Rival hydrocarbon exporters | Brunei, Equatorial Guinea |
| Residual | Rest of World |

Hormuz exposure, exactly as in Section 5.1:

$$e_{ij} = g_i(1-g_j) + g_j(1-g_i)$$

so intra-Gulf trade is unexposed.

### 2.2 The eight sea routes

Base passages (from AIS in the paper; hand-assigned from geography in Track A): `Suez`, `Panama`, `Cape`, `Direct`. Hormuz exposure creates a **composite** of each:

```
r = 1  Suez             r = 5  Hormuz–Suez
r = 2  Panama           r = 6  Hormuz–Panama
r = 3  Cape             r = 7  Hormuz–Cape
r = 4  Direct           r = 8  Hormuz–Direct
```

Baseline route shares split as `ρ[.., r]      = (1 − e_ij) * ρ_geo[.., r]` and
`ρ[.., r+4] = e_ij * ρ_geo[.., r]`, where `ρ_geo` is the passage share ignoring Hormuz.

**Two structural rules from Section 5.1:**
- A Hormuz transit **precedes** the choice of onward passage, so tolls on the strait and the canal **stack multiplicatively**: $\Phi_{q} = \phi_{Hormuz}\phi_{Suez}$ on route 5.
- Closing an onward passage **also closes its Hormuz composite** (Red Sea closure kills routes 1 *and* 5).

### 2.3 Authority → route map

```julia
authorities = (
  Suez    = (routes = [1, 5], collector = :EGY),
  Panama  = (routes = [2, 6], collector = :PAN),
  Hormuz  = (routes = [5, 6, 7, 8], collectors = [:IRN => 0.5, :OMN => 0.5]),
)
```

Combined wedge on route `q`: `Φ[q] = prod(φ[a] for a in authorities_of(q))`. Toll take share on `(q, s)`:

$$t_{q,s} = \frac{\Phi_q^{\chi_s}-1}{\Phi_q^{\chi_s}}, \qquad
\Pi_a = \sum_{q:\,a\in A(q)}\sum_{i,j,s}\frac{\ln\phi_a}{\ln\Phi_q}\, t_{q,s}\, X_{ij,sq} \tag{13}$$

Guard the $\ln\phi_a/\ln\Phi_q$ split for $\Phi_q \to 1$: fall back to equal shares (or an $\varepsilon$-regularized limit) when `log(Φ) < 1e-12`.

### 2.4 Track A route-share assignment

Hand-code `ρ_geo` from region pairs (a small lookup table), e.g.:

| Pair type | Suez | Panama | Cape | Direct |
|---|---|---|---|---|
| Europe/Med ↔ East & South Asia | 0.85 | 0.02 | 0.13 | 0 |
| Europe/Med ↔ Gulf | 0.90 | 0 | 0.10 | 0 |
| East Asia ↔ US East Coast | 0.20 | 0.55 | 0.25 | 0 |
| East Asia ↔ US West Coast | 0 | 0 | 0 | 1 |
| Europe ↔ West/Southern Africa | 0 | 0 | 0.30 | 0.70 |
| Intra-region / all other | 0 | 0 | 0 | 1 |

> **Important:** route shares are **data**, not model output. With $\eta^R = 350$ the implied baseline cost differences are tiny ($\bar d_r \propto \rho_r^{-1/\eta^R}$; a 0.85/0.15 share ratio implies a 0.5% cost gap), so trying to *generate* shares from a distance-based logit would be numerically meaningless. Only *changes* in shares are model-driven — which is exactly what hat algebra needs.

---

## 3. Sectors and Elasticities

| # | Sector | $\chi_s$ | $\eta^M_s$ (Table 4) | $\theta_s$ | Sea-mode share (typ.) |
|---|---|---|---|---|---|
| 1 | Agriculture & food | 0.20 | 4.87 (AFP) | 8 † | 0.95 |
| 2 | Energy / hydrocarbons | 0.30 | 11.08 (ENG) | 15 † | 0.98 |
| 3 | Minerals & bulk | 0.30 | 2.68 (NMM) | 5 † | 0.98 |
| 4 | Chemicals | 0.08 | 3.13 (CHM) | 5 † | 0.90 |
| 5 | Metals | 0.08 | 5.21 (I_S) | 6 † | 0.95 |
| 6 | Machinery & vehicles | 0.03 | 3.17 (OME) | 4 † | 0.85 |
| 7 | Electronics & pharma | 0.01 | 1.67 (ELE) | 5 † | 0.60 |
| 8 | Services & other | 0.00 | 2.44 (pooled) | 3 † | 0.00 (mode 3 only) |

† **$\theta_s$ are placeholders.** The paper uses Fontagné et al. (2022), $\theta_s = \sigma_s - 1$, aggregated to GTAP sectors; those values are not reproduced in the paper. Flag them explicitly in the code with a `# TODO: replace with Fontagné et al. (2022)` comment and expose them as a top-level constant so they are trivial to swap.

$\lambda_s = 0.135$, $\gamma^G_s = 0.44$, $\gamma^M_s = 0$ uniformly across sectors, per Section 5.2.

Also record for the code header: the **effective fee elasticity is $\eta^R \chi_s$** (≈31 at the trade-weighted $\chi$), which is what actually disciplines the toll.

---

## 4. Building an Internally Consistent Baseline

The baseline must satisfy market clearing *exactly*, or the "all hats = 1" test in §8 fails. Build it by solving the levels system once rather than by assembling numbers that almost balance.

```
build_baseline(params; seed = 42)
```

1. **Country sizes** `VA[i]`: assign roughly realistic value added, scaled so world VA hits a realistic level (≈ \$80–100tn) — this is what makes dollar revenue targets ($10.25bn Suez) meaningful.
2. **IO structure**: draw `γ_VA[i,s]` from sector-typical ranges (services high, manufacturing low), and `γ_IO[i,k,s]` from a plausible sector-to-sector matrix with a within-sector diagonal bump. Normalize so `γ_VA + Σ_k γ_IO = 1`.
3. **Final demand shares** `α_fin[j,s]`: sector-typical, higher on services in rich countries.
4. **Trade shares** `π[i,j,s]`: gravity, $\pi_{ij,s} \propto T_{i,s}(w_i \bar d_{ij,s})^{-\theta_s}$, with great-circle distance and a home bias, then column-normalized. Give hydrocarbon exporters high `T` in sector 2 so that the Gulf's Hormuz exposure actually bites.
5. **Mode and route shares** from §2 and the table in §3.
6. **`δ` weights**: `δ[i,j,s,m,r] = route_distance[i,j,m,r] * tons_per_dollar[s]` (bulk-heavy sectors weigh more per dollar). Route distances: hand-coded per passage type in Track A; `searoute-py` output in Track B.
7. **Solve the levels system** for `(Y, X, I)` given `π`, `γ`, `α_fin`, `D` (start with `D = 0`, `Π = 0`):

$$Y_{i,s} = \sum_j \tilde\pi_{ij,s} X_{j,s}, \qquad
X_{j,s} = \sum_k \gamma^{IO}_{j,s,k} Y_{j,k} + \alpha^{fin}_{j,s} I_j, \qquad
I_j = \sum_s \gamma^{VA}_{j,s} Y_{j,s} + \Pi_j + D_j$$

where $\tilde\pi_{ij,s} = \pi_{ij,s}\big(1 - \sum_{m,r}\mu_{ij,sm}\rho_{ij,smr} t_{q(m,r),s}\big)$ is the **net-of-toll** sourcing share. This is linear in `(Y, X)` given `Π`; iterate the outer `Π` loop 3–5 times (it converges immediately, since `Π` is a small share of world income).
8. **Set `VA[i]` to the implied** $\sum_s \gamma^{VA}_{i,s} Y_{i,s}$ (or rescale in step 1 and re-solve) so the baseline is a fixed point of the solver by construction.
9. Compute `Ξ̄[m,r]`, `Ψ̄m[m]`, `Ψ̄` from baseline flows.

---

## 5. The Solver

Three nested loops. All conventions match the repo: damped updates, `round(·, digits=6)` convergence checks, NamedTuple returns, `println` diagnostics.

```julia
solve_counterfactual(base, shock; maxiter=2000, tol=1e-8, damp=0.3)
```

`shock` carries: `d̄_hat[i,j,s,m,r]` (default 1), `route_open[r]::BitVector`, and `φ_new[a]` **levels**. Hats on tolls are `φ̂_a = φ_new[a] / base.φ_base[a]`.

### L3 — Transport block (innermost), given bilateral flows `Xijs`

```
solve_transport_block(base, shock, Xijs) → (d̂[i,j,s], μ', ρ', Ξ̂, Ψ̂m, Ψ̂, Φ, t_share)
```

Iterate to a fixed point on `(Ξ̂, Ψ̂m, Ψ̂)`:

1. Route-level: $\ln \hat d_{smr} = \ln \hat{\bar d} + \chi_s\lambda^R_s\ln\hat\Xi_{mr} + \chi_s\gamma^M_s\ln\hat\Psi_m + \chi_s\gamma^G_s\ln\hat\Psi + \chi_s\ln\hat\phi_{mr}$
2. Within-mode composite and route shares (14)–(15)
3. Across-mode composite and mode shares (16)–(17)
4. Flows `X_smr = Xijs .* μ' .* ρ'`; update `Ξ̂[m,r] = Σδ X_smr / Ξ̄[m,r]`, similarly `Ψ̂`
5. Damp and check

> **Numerical guardrail (critical).** With $\eta^R = 350$, forming `d̂^(-ηR)` directly overflows: `1.5^{-350} ≈ 1e-62`, `0.8^{-350} ≈ 1e34`. **Do the CES in logs with log-sum-exp:**
> ```julia
> ln_d_sm = -logsumexp(log.(ρ0_r) .- ηR .* ln_d_smr) / ηR
> ln_ρ_new = log.(ρ0_r) .- ηR .* (ln_d_smr .- ln_d_sm)
> ```
> Same for the mode nest. Never exponentiate before the LSE.

> **Closed routes.** Do not represent closure as `d̂ = Inf` (it produces `Inf * 0 = NaN` downstream). Use the `route_open` mask: drop closed routes from the LSE and set their share to 0. If *all* routes in a mode close, set that mode's composite to `Inf` in logs and drop it from the mode nest; if all modes close for a pair, the pair's trade goes to zero (guard `θ`-powers against `0^0`).

### L2 — Price block, given `ŵ` and `d̂`

Standard Caliendo–Parro inner loop:

$$\hat c_{i,s} = \hat w_i^{\gamma^{VA}_{i,s}}\prod_k \hat P_{i,k}^{\gamma^{IO}_{i,s,k}}, \qquad
\hat P_{j,s} = \Big[\sum_i \pi_{ij,s}\big(\hat c_{i,s}\hat d_{ij,s}\big)^{-\theta_s}\Big]^{-1/\theta_s}$$

Iterate on `P̂` (contraction; 100–300 iterations at `tol = 1e-10`). Then

$$\pi'_{ij,s} = \pi_{ij,s}\Big(\frac{\hat c_{i,s}\hat d_{ij,s}}{\hat P_{j,s}}\Big)^{-\theta_s}$$

### L1 — Outer loop on wages

Per iteration:

1. `d̂, μ', ρ', Φ, t_share = solve_transport_block(base, shock, Xijs)` using the current `Xijs`
2. `P̂, ĉ, π' = solve_prices(ŵ, d̂, base)`
3. Net-of-toll shares `π̃' = π' .* (1 .- Σ_{m,r} μ' .* ρ' .* t_share)`
4. Solve the expenditure system for `(Y', X', I')` exactly as in §4.7, with `Π'` from (13) evaluated at **levels** `Φ` on counterfactual flows `X'_{ij,s,q}`
5. Labor market clearing gives the wage target directly:
   `ŵ_target[i] = Σ_s γ_VA[i,s] * Y'[i,s] / VA[i]`
6. `ŵ = damp * ŵ_target + (1-damp) * ŵ`; renormalize to the numeraire `Σ_i ŵ_i VA_i = Σ_i VA_i`
7. Converged when `round.(ŵ, digits=6) == round.(ŵ_prev, digits=6)` **and** the goods-market residual `maxabs(Y' - Σ_j π̃' X') < tol`

**Closure:** nominal deficits `D[j]` fixed at baseline; world value added as numeraire (Section 4.4).

### Outputs

```julia
(; ŵ, P̂, d̂, π′, μ′, ρ′, Y′, X′, I′, Π′, Ξ̂, Ψ̂,
   W_hat,        # real income change: Î_j / ∏_s P̂_{j,s}^{α_fin[j,s]}
   ΔW,           # 100*(W_hat - 1), percent — Table 8
   Δτ,           # 100*(Π′_j - Π_j)/I_j, direct rent channel — Table 10
   route_flow_change,   # % change in Σ_{ijs} X_{ij,s,r} by route — Table 9
   toll_revenue)        # by authority, in dollars
```

---

## 6. Toll Calibration (Section 5.1)

This is the part with no analogue in the other repo implementations, so give it its own section in the file.

```julia
best_response_fee(base, authority; φ_others, bracket=(1.0, 3.0), ge=true)
```
Golden-section / Brent maximization of $\Pi_a(\phi_a; \phi_{-a})$ from (13). Two evaluation modes:
- `ge = false` — **partial equilibrium**: hold `Xijs` at baseline and run only the transport block. Cheap; matches the toy model of Section 4.1. Use it to get a starting bracket.
- `ge = true` — full `solve_counterfactual`. Use for the final answer.

```julia
calibrate_route_elasticity(base; target = 10.25e9, bracket = (50.0, 2000.0))
```
Root-find on `ηR` such that the revenue-maximizing Suez fee, **priced in isolation** (all other wedges at their reference level $\phi = 1$), delivers the observed Suez toll revenue. The paper's answer: `ηR = 350`, `φ_Suez = 1.59`, model revenue \$10.5bn (within 2.5% of the \$10.25bn target).

The mapping `ηR → revenue` is monotone (higher $\eta^R$ ⇒ more elastic route choice ⇒ lower optimal fee and revenue), so a simple bisection or `Roots.find_zero` on a bracket suffices. Log the trace so the user can see the calibration curve.

**Then re-baseline.** Apply `(φ_Suez, φ_Panama)` as a counterfactual to the raw ($\phi = 1$) baseline and **overwrite the baseline** with the resulting equilibrium — updating `π, μ, ρ, X, Y, I, Π, Ξ̄, Ψ̄m, Ψ̄` and setting `φ_base`. From then on:
- Counterfactual fees enter route costs as `φ̂ = φ_new/φ_base`
- Toll revenue (13) is always evaluated at the **level** `φ_new`

Keep the raw `φ = 1` baseline in a separate object — Section 6.3 needs it.

**Untargeted validation:** Panama's revenue-maximizing fee at the calibrated `ηR` should land near \$3.6bn with `φ_Panama ≈ 1.16` (observed \$3.18bn, drought-depressed). Print this as a check, not a target.

---

## 7. Scenarios (Section 6)

```julia
scenario_red_sea_closure(base)
```
`route_open[[1, 5]] .= false` — Suez *and* the Hormuz–Suez composite. Cargo reroutes to Cape and Panama, congestion rises on both, global capacity tightens, and Suez rents collapse to zero.

```julia
scenario_hormuz_fee(base, φ; split = [:IRN => 0.5, :OMN => 0.5])
```
`φ_new[:Hormuz] = φ` on routes 5–8, stacking multiplicatively with the Suez/Panama wedges. Run for `φ ∈ {1.05, 1.35, 1.50}`.

```julia
scenario_no_fee_world(raw_base)
```
Re-run **both** scenarios from the raw `φ = 1` baseline (Section 6.3). This is the paper's key comparison and reverses the incidence.

Also worth adding, cheaply:
- **Laffer sweep**: `φ_Hormuz ∈ 1.0:0.05:2.0` → revenue and Hormuz throughput. The paper reports strong convexity (−19%, −60%, −67% throughput) and revenue flattening between 1.35 and 1.50, which the sweep should reproduce as a single-peaked curve.
- **Nash vs. isolation**: solve the two-authority fee game to a fixed point and compare with pricing in isolation. The paper argues they nearly coincide when strategic interaction is weak — a cheap and interesting check of its own claim.
- **Mechanism sensitivity**: $\chi \to 0$ (no pass-through), $\lambda = \gamma^G = 0$ (no congestion feedback), $\eta^R \to \infty$ (cost equalization across routes), mode nest off ($|\mathcal N| = 1$, recovering Section 4.3).

---

## 8. Verification

Run with `julia brockhaus_hinz_serfaty_chokepoint_model.jl`. Each item below should be an explicit printed check.

**Identity / consistency**
1. **Null shock**: `solve_counterfactual(base, no_shock)` returns `ŵ = P̂ = d̂ = 1` and `ΔW = 0` to `1e-8`.
2. **Baseline self-reproduction** (Section 4.4): re-solving the calibrated baseline *at its own fees* leaves route costs and rents unchanged — this is the test that the levels-vs-changes handling of `φ` is right.
3. **Market clearing**: `Y[i,s] ≈ Σ_j π̃[i,j,s] X[j,s]` and `Σ_i (I_i - VA_i - Π_i - D_i) ≈ 0`.
4. **Numeraire**: `Σ_i ŵ_i VA_i = Σ_i VA_i`.
5. **Shares**: `sum(π, dims=1) ≈ 1`, `sum(μ, dims=4) ≈ 1`, `sum(ρ over open routes) ≈ 1`, all non-negative, before and after every scenario.
6. **Toll accounting**: total take under stacked wedges equals the markup at the combined wedge — i.e. the log-share split in (13) sums to $(\Phi^{\chi}-1)/\Phi^{\chi} \cdot X$, and is *strictly less* than charging each authority separately.

**Analytic**
7. **Table 7 reproduction** (exact, no solver): $d = \phi^\chi$ vs. $d = 1+\chi(\phi-1)$ at $\phi = 1.59$ for $\chi \in \{0.01, 0.03, 0.08, 0.30, 1.00\}$ → `1.0046/1.0059`, `1.0139/1.0176`, `1.0376/1.0468`, `1.1483/1.1756`, `1.5855/1.5855`.
8. **Toy monopolist** (Section 4.1): with $d^S = 0$ and a single route pair, the numerical optimum matches $F^\star = d^L(\mu-1)^{-1/\mu}$, and $F^\star \to d^L$ as $\mu \to \infty$.
9. **Additive vs. multiplicative**: run one scenario under both `freight_decomposition = :multiplicative` and `:additive`; welfare effects should differ by a second-order amount and preserve the country ranking.

**Qualitative (the paper's actual claims)**
10. **Incidence with fees**: Egypt is the largest loser under a Red Sea closure, and the rent component `Δτ` accounts for the overwhelming majority of it.
11. **Incidence without fees**: Egypt's loss shrinks by an order of magnitude, while route-lengthened economies (Qatar, UAE, Malta, Gabon-type exporters) deepen. **The two rankings should cross.**
12. **World cost is similar with and without fees** — the toll moves incidence, not the total.
13. **Rerouting**: Suez −100%, Cape and Panama both up (paper: +18.4%, +19.8%).
14. **Hormuz**: losses concentrate on the Gulf; large traders lose < 0.1%; **Oman ≫ Iran** despite the equal revenue split (Oman's own trade bypasses the strait, Iran's does not); rival hydrocarbon exporters gain.
15. **Cross-chokepoint spillover**: under the Hormuz fee, the Hormuz–Suez composite collapses and **Egypt loses** — a fee at one chokepoint taxes the other.
16. **Laffer**: Hormuz revenue is single-peaked in `φ` and near-flat between 1.35 and 1.50.
17. **Limiting cases**: `χ = 0` ⇒ zero welfare effect of any route shock; `λ = γG = 0` ⇒ smaller spillovers to unaffected pairs; `ηR → ∞` ⇒ delivered costs equalize across open routes.

Items 10–17 are the bar for "working". Matching Table 8 to the decimal is **not** achievable without GTAP 11 and the AIS route shares, and the script should say so in its header.

---

## 9. Plots (`graphs/*.pdf`)

1. **Calibration curve** — Suez toll revenue vs. $\phi_{Suez}$, with the \$10.25bn target line and the calibrated optimum marked (this *is* the identification argument, so it should be the first figure).
2. **Hormuz Laffer curve** — revenue and Hormuz throughput vs. $\phi$, with the three scenario fees marked.
3. **Five largest gains/losses per scenario** — horizontal bars (Figures 10, 13, 16).
4. **Route reallocation** — grouped bars by route × scenario (Figure 14, Table 9).
5. **With vs. without fees** — paired bars for the Table 6 countries; the crossing pattern is the paper's punchline.
6. **Rent decomposition** — stacked bars of $\Delta\tau$ and the residual for Egypt, Panama, Iran, Oman, World (Table 10).
7. **Mechanism scatter** — country welfare loss against (a) Suez route exposure and (b) trade-weighted $\chi$, showing that exposure alone does not predict incidence — rent collection does.
8. **Country × scenario heatmap** of $\Delta W$.

Choropleths (Figures 8, 11, 15) need country geometries and are out of scope for a self-contained script; note the omission.

---

## 10. File Layout

```
# Header: paper reference, mechanisms, key equations, what is and is not reproducible
using LinearAlgebra, Statistics, Plots, Random, StatsBase, Roots, Printf

# 1. Parameters and elasticities            (§3)
# 2. Geography: countries, g_i, routes, authorities   (§2)
# 3. Baseline construction                  (§4)
# 4. Transport block (L3)                   (§5)
# 5. Price block (L2) and expenditure system
# 6. Counterfactual solver (L1)             (§5)
# 7. Toll revenue, best responses, calibration   (§6)
# 8. Scenarios                              (§7)
# 9. Verification suite                     (§8)
# 10. Main execution + plots                (§9)
```

New dependency relative to the rest of the repo: **`Roots.jl`** (or hand-rolled bisection — preferable, to keep the dependency list at the repo's usual five packages). Golden-section search for the fee optimum is ~20 lines; write it inline rather than adding `Optim.jl`.

---

## 11. Track B Hooks (documented, not implemented now)

- `load_gtap(path)` → `(π, X, γ_VA, γ_IO, α_fin, D)` from GTAP 11 flat files, aggregated to a user-chosen country/sector set.
- `load_ais_route_shares(path)` → `ρ_geo[i,j,r]` from a CSV of passage shares by directed country pair.
- `export_searoutes.py` — mirrors `FuchsFoongWong-MMN-2026/export_thai_network.R`: uses `searoute-py` (Halili 2026) to compute, for each country pair, port-weighted maritime distance with each passage open and closed, plus the passage classification of the shortest route. Output: one CSV consumed by `load_ais_route_shares`.
- Everything above the loaders is dimension-agnostic; the only code change needed for (160, 65) is the `ρ[i,j,m,r]` storage collapse noted in §1.

---

## 12. Suggested Sequence

| Step | Deliverable | Verified by |
|---|---|---|
| 1 | Parameters, geography, authority map | shares sum to 1; `e_ij` symmetric; intra-Gulf unexposed |
| 2 | Baseline construction | checks 3–5 |
| 3 | Transport block with LSE + closure masks | check 17; route shares respond correctly to a hand-set `φ` |
| 4 | Price + expenditure blocks, full solver | checks 1, 3, 4 |
| 5 | Toll revenue (13) and best responses | checks 6, 8; Table 7 (check 7) |
| 6 | `ηR` calibration and re-baselining | check 2; `φ_Suez ≈ 1.59`; Panama untargeted check |
| 7 | Scenarios + no-fee world | checks 10–15 |
| 8 | Laffer sweep, Nash-vs-isolation, sensitivity | checks 16, 17, 9 |
| 9 | Plots and printed summary tables | visual inspection against Figures 10, 13, 14, 16 |

Steps 1–4 are the bulk of the work and are structurally close to `fuchs_foong_wong_multimodal_model.jl`. Steps 5–6 are the genuinely new part — the monopoly toll layer has no precedent in this repository.

---

## Conventions (following existing implementations)

- **Solver pattern:** nested fixed points with damped updates (0.25–0.5), `round.(x, digits=6)` convergence tests, `maxiter` 2000–3000
- **Returns:** NamedTuples throughout
- **Diagnostics:** `println` for convergence banners (`>>>> Transport block converged <<<<`), `@printf` for result tables
- **Self-contained:** one file, no `include()`s, no external data required to run
- **Reproducible:** `Random.seed!(1)` for the synthetic baseline
- **Comments:** equation numbers (`# eq 13`) throughout, traceable to `README.md` and the paper
- **Header block:** state up front that this is a stylized calibration, that Track A reproduces mechanisms rather than GTAP-exact magnitudes, and list the paper's own four caveats (README §5.5) so they are not lost in translation
