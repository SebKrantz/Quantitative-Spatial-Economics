# Implementation Plan: Allen & Arkolakis (RES 2022) – Welfare Effects of Transportation Infrastructure Improvements

## Context

The user wants a Julia implementation of the Allen & Arkolakis (2022) model from "The Welfare Effects of Transportation Infrastructure Improvements." This paper develops a quantitative spatial framework with **endogenous transportation costs and traffic congestion**, yielding analytical expressions for transport costs (via Leontief inverse of the adjacency matrix), traffic flows (gravity equation), and equilibrium spatial distribution of economic activity. The key innovation vs. standard spatial models is that transport costs are endogenous: agents choose routes optimally through a network, and traffic congestion raises costs on heavily-used links.

The implementation follows the coding conventions established in:
- `/Santamaria-JEEA-2026/santamaria_infrastructure_model.jl` (single file, params at top, detailed comments, returns vectors/matrices)
- `/QRE-HoRaUE-2025/simple_qre_model.jl` and `workhorse_qre_framework.jl`

## Output File

`/AllenArkolakis-RES-2022/allen_arkolakis_traffic_model.jl` (~900–1100 lines)

## Model Summary

Two model variants sharing identical transport cost / traffic / congestion structure:

1. **Economic Geography Model**: Agents choose location, trade goods. Unknowns: {y_i, l_i} (income share, labor share). Eqs (10)-(11), with congestion: (28)-(29). Exact-hat: (36)-(37).
2. **Urban Model**: Agents choose residence + workplace + route. Unknowns: {l^R_i, l^F_i}. Eqs (19)-(20), with congestion: (30)-(31). Exact-hat: (38)-(39).

**Key equations**:
- Transport costs: τ_{ij} = b_{ij}^{-1/θ} where B = (I - A)^{-1}, A = [t_{kl}^{-θ}] (eq 21, Leontief inverse)
- Traffic gravity: Ξ_{kl} = t_{kl}^{-θ} × P_k^{-θ} × Π_l^{-θ} (eq 24)
- Congestion: t_{kl} = t̄_{kl} × Ξ_{kl}^λ (eq 25)
- Combined: t_{kl} = t̄_{kl}^{1/(1+θλ)} × P_k^{-θλ/(1+θλ)} × Π_l^{-θλ/(1+θλ)} (eq 26)
- Welfare: W̄ = χ^{-1/θ}
- Infrastructure improvement: t̂̄_{kl} = lanes_kl^{-λ} (eq 42)

## Implementation Structure

### 1. Header & Packages (~20 lines)
```julia
using LinearAlgebra, Statistics, Plots, Random
```
No Graphs.jl needed (transport costs computed via matrix Leontief inverse, not Dijkstra).

### 2. Parameters (~30 lines)

**Economic Geography (US Highway):**
- θ = 8 (trade elasticity / Fréchet shape)
- α = 0.1 (productivity externality)
- β = -0.3 (amenity externality / housing)
- λ_congestion = 0.092 (traffic congestion parameter)
- L̄ = 100 (aggregate labor, for synthetic example)

**Urban (Seattle):**
- θ_urban = 6.83
- α_urban = -0.12
- β_urban = -0.1
- λ_urban = 0.071

**Grid:**
- N_side = 5, N = 25
- t̄_kl = 1.5 for connected links (as in Figure 1)

### 3. Network Construction (~40 lines)
`build_grid_network(N_side; t_bar=1.5)` → Returns:
- `adjacency`: N×N matrix of t̄_{kl} (∞ for unconnected)
- `coords`: N×2 position matrix
- `connected`: N×N Bool matrix

### 4. Transport Cost Functions (~60 lines)

`compute_transport_costs(t_bar_matrix, theta)`
- Compute A = [t_{kl}^{-θ}] (setting disconnected to 0)
- Check spectral radius of A < 1
- B = (I - A)^{-1} (Leontief inverse)
- τ_{ij} = b_{ij}^{-1/θ} (eq 21)
- Returns: τ matrix (N×N), B matrix, A matrix

`compute_transport_costs_with_congestion(t_bar_matrix, P_inv, Pi_inv, theta, lambda)`
- First compute endogenous t_{kl} from eq (26): t_{kl} = t̄_{kl}^{1/(1+θλ)} × P_k^{-θλ/(1+θλ)} × Π_l^{-θλ/(1+θλ)}
- Then compute τ_{ij} via Leontief inverse
- Returns: τ matrix, traffic Ξ matrix, t matrix

### 5. Equilibrium Solver – Economic Geography (~120 lines)

`solve_econ_geo_equilibrium(A_bar, u_bar, t_bar_matrix, theta, alpha, beta, lambda, Lbar; ...)`

Solve eqs (28)-(29) for {y_i, l_i, χ}:
- These are 2N equations in 2N unknowns (y_i, l_i) plus scalar χ
- χ = (L̄^{α+β} / W̄)^θ determined by normalization Σ y_i = 1 and Σ l_i = 1
- Use damped fixed-point iteration (following Santamaria pattern):
  1. Initialize y, l uniformly
  2. Compute RHS of (28)-(29) given current y, l
  3. Update with damping
  4. Enforce normalization
  5. Repeat until convergence

Returns: NamedTuple{y, l, chi, tau, traffic, market_access_P, market_access_Pi, welfare}

### 6. Equilibrium Solver – Urban Model (~120 lines)

`solve_urban_equilibrium(A_bar, u_bar, t_bar_matrix, theta, alpha, beta, lambda, Lbar; ...)`

Solve eqs (30)-(31) for {l^R_i, l^F_i, χ}:
- Analogous structure to economic geography
- Uses same damped fixed-point

Returns: NamedTuple{l_R, l_F, chi, tau, traffic, welfare}

### 7. Traffic and Trade/Commuting Prediction (~60 lines)

`predict_trade_from_traffic(Xi, Y, E)` (eq 34)
- D^X = diag(½(Y+E) + ½(Ξ1 + Ξ'1))
- C^X = (D^X - Ξ)^{-1}
- X_{ij} = c^X_{ij} × Y_i × E_j

`predict_commuting_from_traffic(Xi, L_R, L_F)` (eq 35)
- D^L = diag(½(L^R+L^F) + ½(Ξ1 + Ξ'1))
- C^L = (D^L - Ξ)^{-1}
- L_{ij} = c^L_{ij} × L^R_i × L^F_j

### 8. Link Intensity (~30 lines)

`compute_link_intensity(tau, t, theta)` (eq 23)
- π^{kl}_{ij} = (τ_{ij} / (τ_{ik} × t_{kl} × τ_{lj}))^θ
- Returns 4D array or function

### 9. Exact-Hat Counterfactual Solver – Economic Geography (~100 lines)

`counterfactual_econ_geo(Xi, Y, E, t_hat_bar, theta, alpha, beta, lambda; ...)`

Solve eqs (36)-(37) for {ŷ_i, l̂_i, χ̂}:
- Input: observed traffic Ξ_{ij}, income Y_i = E_i, infrastructure change t̂̄_{kl}
- Weight shares: s^out_{ij} = Ξ_{ij}/(E_i + Σ_k Ξ_{ik}), s^in_{ij} = Ξ_{ji}/(Y_i + Σ_k Ξ_{ki})
- Own shares: s^own_out_i = E_i/(E_i + Σ_k Ξ_{ik}), s^own_in_i = Y_i/(Y_i + Σ_k Ξ_{ki})
- Fixed-point iteration on {ŷ, l̂} with normalization Σ ŷ_i l̂_i = 1 (income shares × labor shares sum correctly)
- Welfare change: Ŵ = χ̂^{-1/θ}

Returns: NamedTuple{y_hat, l_hat, chi_hat, welfare_change}

### 10. Exact-Hat Counterfactual Solver – Urban (~100 lines)

`counterfactual_urban(Xi, L_R, L_F, t_hat_bar, theta, alpha, beta, lambda; ...)`

Solve eqs (38)-(39) for {l̂^R_i, l̂^F_i, χ̂}:
- Analogous structure to economic geography
- Weight shares use L^F_i and L^R_i instead of E_i and Y_i

Returns: NamedTuple{l_R_hat, l_F_hat, chi_hat, welfare_change}

### 11. Welfare Elasticity Computation (~60 lines)

`compute_welfare_elasticities(Xi, Y_or_LR, E_or_LF, t_bar_matrix, theta, alpha, beta, lambda; model=:econ_geo, delta=0.01)`
- For each link (k,l) in the network:
  - Construct t̂̄ = I everywhere except link (k,l) and (l,k) where t̂̄ = 1 - δ
  - Solve counterfactual
  - Record ∂ln W̄ / ∂ln t̄_{kl} ≈ (Ŵ - 1) / δ
- Returns: matrix of welfare elasticities per link

### 12. Main Execution Block (~200 lines)

**Part A: Synthetic Grid Example (replicating Figure 1)**
1. Build 5×5 grid, identical fundamentals (Ā_i = ū_i = 1)
2. Solve urban equilibrium with λ=0 (no congestion) → Figure 1(a)
3. Solve urban equilibrium with λ=0.05 → Figure 1(b)
4. Solve with λ=0.05, L̄=1000 → Figure 1(c)
5. Solve with λ=0.05, L̄=10000 → Figure 1(d)
6. Plot: heatmaps of population distribution + traffic on edges for each case
7. Verify: congestion pushes activity to edges, scale dependence

**Part B: Economic Geography Counterfactuals**
1. Solve equilibrium in levels on the grid
2. Compute traffic flows Ξ
3. Run counterfactual: 1% improvement to central link
4. Compute welfare elasticities for all links
5. Compare with/without congestion

**Part C: Urban Counterfactuals**
1. Solve equilibrium in levels on the grid
2. Run counterfactual: 1% improvement to central link
3. Compute welfare elasticities for all links
4. Compare with/without congestion

**Part D: Diagnostics**
1. Uniqueness check: α+β ≤ 0 (econ geo), α,β ≤ ½(1/θ - λ) (urban)
2. Predict trade/commuting from traffic → correlation check
3. Spectral radius check on adjacency matrix
4. Print summary statistics

### 13. Plots (~60 lines)
- Heatmap of equilibrium population (grid layout)
- Traffic on edges (colored by intensity)
- Welfare elasticity map
- Scatter: welfare elasticity with vs. without congestion
- Save to `/AllenArkolakis-RES-2022/graphs/`

## Key Implementation Notes

1. **No Dijkstra needed**: Unlike Santamaria, this paper uses the Leontief inverse (I-A)^{-1} to compute transport costs analytically. This is a key feature – the stochastic route choice with Fréchet shocks yields a closed-form matrix inverse.

2. **Adjacency matrix A**: Set a_{ij} = t_{ij}^{-θ} for connected links, 0 for unconnected. Must check spectral radius < 1. For the grid with t̄=1.5 and θ=4-8, this should be fine since 1.5^{-4} = 0.198 and each row has at most 4 nonzero entries.

3. **Fixed-point with congestion**: The full equilibrium with congestion (eqs 28-31) involves a complex feedback loop. The key insight is that despite the feedback, the system has the same dimensionality (2N equations, 2N unknowns). Use damped iteration with appropriate starting values.

4. **Exact-hat algebra**: The counterfactual system (eqs 36-39) has the same mathematical structure as the levels system but uses observed data (traffic Ξ, income/population) as weights. This is preferable for empirical work since it doesn't require knowledge of fundamentals {Ā_i, ū_i}.

5. **Scale dependence**: With λ > 0, equilibrium depends on L̄ (unlike standard spatial models). This is a key result demonstrated in the synthetic example.

6. **Welfare**: W̄ = χ^{-1/θ} where χ = (L̄^{α+β}/W̄)^θ. In counterfactuals, welfare change Ŵ = χ̂^{-1/θ}.

## Verification

1. **Run**: `julia AllenArkolakis-RES-2022/allen_arkolakis_traffic_model.jl`
2. **Check**:
   - Uniqueness conditions satisfied for all parameter constellations used
   - Spectral radius of A < 1
   - Equilibrium shares sum to 1 (Σ y_i = Σ l_i = 1)
   - No congestion case: scale-invariant (changing L̄ doesn't change distribution)
   - With congestion: center locations lose activity as L̄ increases
   - Welfare elasticities: positive for all links, highest for central/hub links
   - Congestion reduces welfare gains (elasticities lower with λ>0)
   - Trade/commuting prediction from traffic: positive correlation
   - PDF plots generated in `graphs/` folder
