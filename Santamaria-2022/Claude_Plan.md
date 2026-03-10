# Implementation Plan: Santamaria (JEEA 2026) — Spatial Model with Endogenous Infrastructure

## Context

Implement the quantitative spatial trade model from Santamaria (2026) "Reshaping Infrastructure: Evidence from the Division of Germany" as a self-contained Julia script. The model extends Redding (2016) / Allen-Arkolakis (2014) with a benevolent government choosing infrastructure investments across locations on a transport network. Key novelty: transport costs depend on the entire network via least-cost paths, and infrastructure investments are location-level with spillovers to all routes transiting a location.

**File:** `/Users/sebastiankrantz/Documents/IfW Kiel/Quantitative-Spatial-Economics/Santamaria-JEEA-2026/santamaria_infrastructure_model.jl`

**Reference implementations for coding conventions:** `QRE-HoRaUE-2025/simple_qre_model.jl`, `QRE-HoRaUE-2025/workhorse_qre_framework.jl`, `Redding-JIE-2016/GTFM_JuliaPrograms/regions/functions/solveLw.jl`

---

## Model Summary

**Agents:** L mobile workers choosing among N locations with heterogeneous (Fréchet) preferences. Firms produce differentiated varieties under monopolistic competition (Krugman 1980). A government allocates infrastructure investments to maximize aggregate welfare.

**Key parameters:**
| Symbol | Description | Value |
|--------|-------------|-------|
| σ | Elasticity of substitution | 7 |
| α | Expenditure share on tradables | 0.7 |
| ε | Fréchet shape (migration elasticity) | 3 |
| γ | Returns to infrastructure investment | 0.84 |
| F | Fixed cost per firm | 1 (normalized) |
| Z | Government budget | calibrated |
| c | Marginal construction cost | 1 (uniform) |

**Uniqueness condition (eq 15):** `1/(σ-1) - 1/(α*ε) - (1-α)/α ≤ 0`
With default params: `1/6 - 1/2.1 - 3/7 ≈ 0.167 - 0.476 - 0.429 = -0.738 ≤ 0` ✓

---

## File Structure

### 1. Packages
```julia
using Random, Statistics, Plots, LinearAlgebra, StatsBase, Graphs, Optim
```
- `Graphs.jl` for network representation and Dijkstra's shortest path algorithm
- `Optim.jl` for the government's infrastructure optimization (interior-point or L-BFGS-B with bounds)

### 2. Parameters (top of file)
```
N_side = 5          # Grid side length (N = N_side^2 = 25 locations)
N = N_side^2        # Total locations
sigma = 7.0         # Elasticity of substitution
alpha = 0.7         # Expenditure share on tradables
epsilon = 3.0       # Fréchet shape parameter
gamma = 0.84        # Infrastructure elasticity (decreasing returns)
F = 1.0             # Fixed cost (normalized)
Lbar = 1000.0       # Total labor supply
c = 1.0             # Marginal construction cost (uniform)
Z = ...             # Government budget (calibrated)
```

### 3. Geography: Network Construction

Generate a regular grid graph with 4-nearest-neighbor connectivity (like Figure 3 in paper):

```
function build_grid_network(N_side)
```

- Place N = N_side² locations on unit square grid
- Connect each to 4 nearest neighbors (up/down/left/right)
- Compute Euclidean distances for each edge
- Return: `coords` (Nx2), `graph` (SimpleGraph), `edge_distances` (Dict mapping edges to distances)

Also generate:
- `H` (land endowments): log-normal, normalized (or set = area ∝ 1/N)
- `A` (productivities): log-normal, geo-mean normalized to 1

### 4. Transport Cost Functions

#### `compute_edge_weights(graph, edge_distances, phi, gamma)`
Compute the ad-valorem shipping cost along each edge (eq 7):
```
w(x,y) = 0.5 * (d_{x,y}/phi[r(x)]^gamma + d_{x,y}/phi[r(y)]^gamma)
```
where `r(x)` maps vertex x to its location/region.

Returns a weighted adjacency / edge-weight structure for shortest path computation.

#### `compute_transport_costs(graph, edge_weights, locations)`
Use Dijkstra's algorithm on the weighted graph to compute least-cost paths between all location pairs. Transport cost (eq 10):
```
T[n,i] = 1 + sum of w(x,y) along least-cost path from n to i
T[n,n] = 1  (normalized)
```
Returns: NxN transport cost matrix `T`, and optionally the path indicator matrices `I_{n,i}^{x,y}` (needed for FOCs / gradient).

### 5. Spatial Equilibrium Solver

#### `solve_equilibrium(A, H, T, sigma, alpha, epsilon, F, Lbar; maxiter=2000, tol=1e-6)`

Solve eqs 11-14 for `{w, P, r, L}` given infrastructure (transport costs).

**Algorithm** (nested fixed-point, following solveLw.jl pattern):

**Outer loop** (population, max 2000 iterations):
1. Given `L`, compute number of varieties: `M_i = L_i / (σF)`

**Inner loop** (wages, max 2000 iterations):
2. Compute delivered prices: `p_{n,i} = (σ/(σ-1)) * (w_i/A_i) * T_{n,i}`
3. Price index (eq 12): `P_n^{1-σ} = Σ_i (L_i/(σF)) * p_{n,i}^{1-σ}`
4. Trade shares: `π_{n,i} = (L_i/(σF)) * p_{n,i}^{1-σ} / P_n^{1-σ}`  (share of n's expenditure on i's goods = share of i exporting to n)
   - Note: In the paper, `X_{n,i}` is expenditure of *n* on goods from *i*. The trade share `π_{n,i} = X_{n,i} / (α * v_n * L_n)`.
   - But for the balanced trade condition (eq 11), we need: `w_n * L_n = Σ_i π_{i,n} * (w_i * L_i / α)` — i.e., income of n = sum over all destinations i of n's share in i's expenditure × i's total expenditure.

   Actually, let's be more careful. Define:
   - Total income in n: `v_n * L_n = w_n * L_n / α`  (since land rents are redistributed)
   - Total expenditure on tradables in n: `α * v_n * L_n = w_n * L_n`
   - Trade share (share of n's tradable expenditure on goods from i): `π_{n,i} = (M_i * p_{n,i}^{1-σ}) / P_n^{1-σ}`
   - Goods market clearing (eq 11): `w_n * L_n = Σ_i π_{i,n} * w_i * L_i`
     (labor income of n = total spending by all locations on n's goods)

5. Check balanced trade: `income_n = w_n * L_n`, `expend_n = Σ_i π_{i,n} * w_i * L_i`
6. Update wages: `w_e = w .* (expend ./ income) .^ (1/σ)`, damp: `w = 0.25*w_e + 0.75*w`
7. Normalize: `w ./= geomean(w)`
8. Convergence check: `round.(income, digits=6) == round.(expend, digits=6)`

**Back to outer loop:**
9. Compute rents (eq 13): `r_n = ((1-α)/α) * w_n * L_n / H_n`
10. Compute real income / indirect utility: `U_n = v_n / (P_n^α * r_n^{1-α})`  where `v_n = w_n/α`
11. Update population (eq 14): `L_e_n = (U_n^ε / Σ_k U_k^ε) * Lbar`
12. Damp: `L = 0.25*L_e + 0.75*L`, or `L = L .* (L_e ./ L) .^ (1/ε)` then damp
13. Convergence check on L

**Returns:** NamedTuple `(; w, L, P, r, tradesh, welfare, U_tilde)` where:
- `U_tilde` (eq 6): `δ * [Σ_n (v_n / P_n^α * r_n^{1-α})^ε]^{1/ε}` with `δ = Γ(ε/(ε-1))`

### 6. Government Optimization

#### `compute_welfare(phi, graph, edge_distances, A, H, ...)`
Wrapper: given `phi` vector → compute edge weights → transport costs → solve equilibrium → return `-U_tilde` (negative for minimization).

#### `optimize_infrastructure(graph, edge_distances, A, H, params...; Z, phi_init, phi_lb=ones(N))`

Solve the government's problem (eq 24):
```
max_{phi} U_eq(phi)   s.t.  Σ c_n * phi_n ≤ Z,  phi_n ≥ 1 ∀n
```

**Algorithm** (following Section III.H):
1. Start from `phi_0 = ones(N)` (or `phi_0 = Z/N * ones(N)` ensuring budget binds) or provided `phi_init`
2. Use `Optim.jl` with:
   - Objective: `-U_tilde(phi)` (minimize negative welfare)
   - Constraints: box constraints `phi_n ≥ 1` (can use Fminbox), budget handled via penalty or reparametrization
   - **Budget constraint handling:** Reparametrize as `phi = 1 .+ (Z - N) * softmax(theta)` so that `Σ(phi-1) = Z - N` automatically (or equivalently `Σ phi = Z` if c=1). This converts the constrained problem into an unconstrained one.
   - Alternatively: use `NelderMead` or `LBFGS` with a penalty term for budget violation, or use the `IPNewton` method from Optim.jl.
   - **Practical approach:** Use projected gradient or simply: allocate budget proportional to marginal welfare gains, iteratively. Since the paper acknowledges finding local optima, a simpler approach is acceptable.

   **Simplest robust approach:**
   - Reparametrize: `phi_n = 1 + s_n * (Z - N*c) / (c * N)` where `s_n = exp(theta_n) / Σ exp(theta_k)` are softmax shares of the *excess* budget. This ensures `phi_n ≥ 1` and `Σ c * phi_n = Z` automatically.
   - Optimize over unconstrained `theta` using `Optim.optimize(..., NelderMead())` or `LBFGS()`.
   - Each evaluation: theta → phi → transport costs → equilibrium → welfare.

**Returns:** optimal `phi`, equilibrium at optimum, welfare level.

### 7. Counterfactual Analysis

#### `counterfactual_infrastructure(phi_baseline, phi_counterfactual, ...)`
Compare welfare under two different infrastructure vectors:
- Solve equilibrium under each
- Compute welfare change: `ΔW = (U_tilde_cf - U_tilde_base) / U_tilde_base * 100`

Implement key counterfactuals from the paper:
1. **No infrastructure** (`phi = ones(N)`) vs. optimal → total gains from infrastructure
2. **Uniform infrastructure** vs. optimal → gains from spatial targeting
3. **Constrained optimal** (some phi fixed at pre-existing levels) vs. **unconstrained optimal** → cost of path dependence
4. **Shock simulation**: Remove some locations from the trade network (like division), re-optimize → infrastructure reshaping

### 8. Diagnostics and Derived Quantities

#### `compute_market_access(tradesh, w, L, P, T, sigma)`
- Trade flows: `X_{n,i} = π_{n,i} * w_n * L_n`  (eq 3)
- Market access (outward, eq from Donaldson-Hornbeck): `MA_n = Σ_i T_{n,i}^{1-σ} * E_i / MA_i` (iterative)
- Centrality: for each location k, count how many shortest paths transit k (betweenness centrality from Graphs.jl, or from path indicators)
- Network statistics: average path length, diameter

#### `invert_productivity(L_data, H, T, sigma, alpha, epsilon, F, Lbar)`
Given observed population distribution, solve for the productivity vector `A` that rationalizes it in equilibrium (used in calibration). This is the inverse problem: find `A` such that `solve_equilibrium(A, ...).L ≈ L_data`.

**Algorithm:** Iterative. Start with `A = ones(N)`. Solve equilibrium → get `L_model`. Update `A` using the excess ratio `A_new = A .* (L_data ./ L_model) .^ step`, iterate until `L_model ≈ L_data`.

### 9. Main Execution

```julia
# 1. Build network, set parameters, generate geography
# 2. Print uniqueness condition check
# 3. Solve equilibrium with uniform infrastructure (phi=1)
#    Print: convergence, welfare, trade balance check
# 4. Solve for optimal infrastructure (government problem)
#    Print: optimal phi, welfare gain, infrastructure allocation pattern
# 5. Counterfactuals:
#    a. Compare optimal vs uniform infrastructure → gains from targeting
#    b. Simulate "division" shock (remove top row of grid from trade)
#       Re-solve optimal infrastructure → show reshaping
#    c. Constrained vs unconstrained re-optimization → cost of path dependence
# 6. Print summary statistics: population Gini, trade openness, welfare gains
# 7. Plots
```

### 10. Plots (saved to PDF)

1. **Network visualization**: Grid with edge thickness ∝ infrastructure quality, node size ∝ population
2. **Optimal infrastructure**: Heatmap of `phi` across grid locations (cf. Figure 4 in paper)
3. **Population distribution**: Bar chart baseline vs. counterfactual
4. **Division counterfactual**: Side-by-side heatmaps of optimal phi before vs. after "division"
5. **Welfare comparison**: Bar chart of welfare levels under different infrastructure scenarios
6. **Exports vs centrality**: Scatter of location-level exports and betweenness centrality vs. infrastructure investment (illustrating eq 23 intuition)

---

## Key Equations Reference

| Eq | Description | Formula |
|----|-------------|---------|
| 1 | Utility | `U_n(ω) = b_n(ω) * (C_n/α)^α * (H_n/(1-α))^{1-α}` |
| 2 | Pricing | `p_i = (σ/(σ-1)) * w_i/A_i` |
| 3 | Trade demand | `X_{n,i} = (L_i/(σF)) * p_{n,i}^{1-σ} / P_n^{1-σ} * α * X_n` |
| 4 | Population | `L_n/L = (v_n / P_n^α r_n^{1-α})^ε / Σ(...)^ε` |
| 6 | Expected utility | `Ũ = δ [Σ (v_n/P_n^α r_n^{1-α})^ε]^{1/ε}` |
| 7 | Edge shipping cost | `w(x,y) = 0.5*(d_{x,y}/phi_{r(x)}^γ + d_{x,y}/phi_{r(y)}^γ)` |
| 10 | Transport cost | `T_{n,i} = 1 + Σ_ℓ I_{n,i}^{x,y} * w(x,y)` |
| 11 | Balanced trade | `w_n L_n = Σ_i (L_i/(σF))(σ/(σ-1) w_i/A_i T_{n,i})^{1-σ} P_i^{σ-1} w_i L_i` |
| 12 | Price index | `P_n^{1-σ} = Σ_i (L_i/(σF))(σ/(σ-1) w_i/A_i T_{n,i})^{1-σ}` |
| 13 | Land rent | `r_n = ((1-α)/α) * w_n L_n / H_n` |
| 14 | Residential choice | `L_n/L = (v_n/P_n^α r_n^{1-α})^ε / Σ(...)^ε` |
| 15 | Uniqueness | `1/(σ-1) - 1/(αε) - (1-α)/α ≤ 0` |
| 16 | Budget | `Σ c_n phi_n ≤ Z` |
| 23 | Infrastructure FOC | `phi_k^{γ+1} ∝ exports(k) + centrality(k)` |

---

## Conventions (following existing implementations)

- **Solver pattern**: Nested fixed-point iteration, damped updates (0.25/0.75), convergence via `round.(x, digits=6)` equality, max 2000 iterations
- **Normalization**: Wages by geometric mean, population sums to Lbar
- **Returns**: NamedTuples for all function outputs
- **Diagnostics**: `println()` for convergence messages, `@show` for key results
- **Self-contained**: All functions defined within the single file, no external includes
- **Reproducible**: `Random.seed!(1)` for geography generation
- **Comments**: Equation references (eq X) throughout for traceability to paper

---

## Verification

1. Run with `julia santamaria_infrastructure_model.jl`
2. **Uniqueness condition**: Verify `1/(σ-1) - 1/(αε) - (1-α)/α < 0` with default parameters
3. **Trade balance**: `Σ_i X_{n,i} ≈ w_n * L_n` for all n (balanced trade)
4. **Welfare equalization**: All workers have same expected utility (by construction of Fréchet model)
5. **Optimal > uniform**: Welfare under optimal infrastructure > welfare under uniform phi
6. **Budget binds**: `Σ c_n * phi_n ≈ Z` at optimum
7. **Division reshaping**: After removing locations, optimal infrastructure shifts away from removed boundary (qualitative check matching paper's findings)
8. **Path dependence cost**: Unconstrained optimal > constrained optimal (welfare)
9. **Gravity**: Trade flows decrease with distance (log-linear, negative coefficient)
10. **Infrastructure correlates**: phi correlates positively with exports and centrality (eq 23 intuition)
