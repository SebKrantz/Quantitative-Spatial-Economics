# Spatial Model with Endogenous Infrastructure Investment
# Based on Santamaria (JEEA 2026) "Reshaping Infrastructure: Evidence from the Division of Germany"
#
# Implements a multi-region quantitative trade model with:
# - Monopolistic competition à la Krugman (1980)
# - Mobile workers with Fréchet heterogeneous preferences
# - Network-based transport costs with least-cost path routing
# - Endogenous infrastructure investment by a benevolent government
#
# The model extends Redding (2016) / Allen-Arkolakis (2014) by introducing
# location-level infrastructure investments that affect the entire transport
# cost matrix through shortest-path routing on a network graph.

using Random
using Statistics
using Plots
using LinearAlgebra
using StatsBase: geomean
using Graphs
using Optim
using SpecialFunctions: gamma as Γ

# ============================================================
# Parameters (Table D.3 in paper)
# ============================================================

N_side = 5                  # Grid side length
N = N_side^2                # Number of locations (25)
sigma = 7.0                 # Elasticity of substitution across varieties (σ > 1)
alpha = 0.7                 # Expenditure share on tradable goods
epsilon = 3.0               # Fréchet shape parameter (migration elasticity)
gamma_infra = 0.84          # Returns to infrastructure investment (γ < 1 = decreasing returns)
F = 1.0                     # Fixed cost per firm (normalized)
Lbar = 1000.0               # Total labor supply
c_cost = 1.0                # Marginal construction cost (uniform across locations)
Z_budget = 2.5 * N          # Government infrastructure budget

# Uniqueness condition (eq 15): 1/(σ-1) - 1/(α*ε) - (1-α)/α ≤ 0
# Agglomeration force λ = 1/(σ-1), Dispersion force η = -1/(αε) - (1-α)/α
lambda_agg = 1 / (sigma - 1)
eta_disp = -1 / (alpha * epsilon) - (1 - alpha) / alpha
uniqueness_cond = lambda_agg + eta_disp
println("Uniqueness condition (eq 15): λ + η = $(round(lambda_agg, digits=4)) + $(round(eta_disp, digits=4)) = $(round(uniqueness_cond, digits=4))")
uniqueness_cond <= 0 ? println("  Unique stable equilibrium guaranteed.") : println("  WARNING: Uniqueness not guaranteed!")

# Fréchet constant δ = Γ(ε/(ε-1)) for expected utility (eq 6)
delta_frechet = Γ(epsilon / (epsilon - 1))

# ============================================================
# Network Construction (Section III.E)
# ============================================================

"""
    build_grid_network(N_side; spacing=1.0)

Build a regular grid network with 4-nearest-neighbor connectivity.
Each of the N = N_side² locations is connected to its up/down/left/right neighbors.
This matches the example geography in Figure 3 of the paper.

Returns NamedTuple with:
- graph: SimpleGraph with N vertices
- coords: Nx2 matrix of (x,y) coordinates
- edge_distances: Dict mapping Edge => Euclidean distance
"""
function build_grid_network(N_side; spacing=1.0)
    N = N_side^2
    g = SimpleGraph(N)
    coords = zeros(N, 2)

    # Place locations on a regular grid
    for i in 1:N_side, j in 1:N_side
        idx = (i - 1) * N_side + j
        coords[idx, :] = [(j - 1) * spacing, (N_side - i) * spacing]
    end

    # Connect 4-nearest neighbors (up, down, left, right)
    for i in 1:N_side, j in 1:N_side
        idx = (i - 1) * N_side + j
        # Right neighbor
        if j < N_side
            add_edge!(g, idx, idx + 1)
        end
        # Down neighbor
        if i < N_side
            add_edge!(g, idx, idx + N_side)
        end
    end

    # Compute edge distances
    edge_distances = Dict{Edge{Int}, Float64}()
    for e in edges(g)
        d = sqrt(sum((coords[src(e), :] .- coords[dst(e), :]) .^ 2))
        edge_distances[e] = d
        edge_distances[Edge(dst(e), src(e))] = d  # symmetric
    end

    return (; graph=g, coords, edge_distances)
end

# ============================================================
# Transport Cost Functions (Section III.E, eqs 7, 10)
# ============================================================

"""
    compute_edge_weights(graph, edge_distances, phi, gamma_infra)

Compute the ad-valorem shipping cost along each edge (eq 7):
    w(x,y) = 0.5 * (d_{x,y} / φ_{r(x)}^γ + d_{x,y} / φ_{r(y)}^γ)

In this grid network, each vertex IS a location, so r(x) = x.
Higher infrastructure φ reduces shipping costs. γ < 1 implies decreasing returns.

Returns Dict mapping Edge => shipping cost.
"""
function compute_edge_weights(graph, edge_distances, phi, gamma_infra)
    weights = Dict{Edge{Int}, Float64}()
    for e in edges(graph)
        i, j = src(e), dst(e)
        d = edge_distances[e]
        # Eq 7: average cost across the two endpoint infrastructure levels
        w_ij = 0.5 * (d / phi[i]^gamma_infra + d / phi[j]^gamma_infra)
        weights[e] = w_ij
        weights[Edge(j, i)] = w_ij  # symmetric
    end
    return weights
end

"""
    compute_transport_costs(graph, edge_weights, N)

Compute the NxN transport cost matrix using Dijkstra's shortest path algorithm.
Transport cost between locations n and i (eq 10):
    T_{n,i} = 1 + Σ_ℓ I_{n,i}^{x,y} * w(x,y)
where the sum is along the least-cost path from n to i.
T_{n,n} = 1 (normalized, free intra-location trade).

Returns NxN matrix T.
"""
function compute_transport_costs(graph, edge_weights, N)
    # Build a weight matrix for Dijkstra
    n_vertices = nv(graph)
    wt_matrix = zeros(n_vertices, n_vertices)
    for e in edges(graph)
        i, j = src(e), dst(e)
        w = edge_weights[e]
        wt_matrix[i, j] = w
        wt_matrix[j, i] = w
    end

    # Compute shortest paths from each location
    T = ones(N, N)  # T[n,n] = 1
    for n in 1:N
        ds = dijkstra_shortest_paths(graph, n, wt_matrix)
        for i in 1:N
            if i != n
                T[n, i] = 1.0 + ds.dists[i]  # eq 10: T = 1 + least-cost path distance
            end
        end
    end

    return T
end

"""
    compute_transport_costs_and_paths(graph, edge_weights, N)

Same as compute_transport_costs but also returns the shortest path predecessor
matrices (for computing betweenness centrality and infrastructure gradients).
"""
function compute_transport_costs_and_paths(graph, edge_weights, N)
    n_vertices = nv(graph)
    wt_matrix = zeros(n_vertices, n_vertices)
    for e in edges(graph)
        i, j = src(e), dst(e)
        w = edge_weights[e]
        wt_matrix[i, j] = w
        wt_matrix[j, i] = w
    end

    T = ones(N, N)
    parents = zeros(Int, N, N)  # parents[n, i] = predecessor of i on shortest path from n

    for n in 1:N
        ds = dijkstra_shortest_paths(graph, n, wt_matrix)
        for i in 1:N
            if i != n
                T[n, i] = 1.0 + ds.dists[i]
            end
        end
        parents[n, :] = ds.parents
    end

    return T, parents
end

# ============================================================
# Spatial Equilibrium Solver (Section III.F, eqs 11-14)
# ============================================================

"""
    solve_equilibrium(A, H, T, sigma, alpha, epsilon, F, Lbar;
                      maxiter=2000, tol=1e-6, verbose=true)

Solve the spatial equilibrium for given fundamentals and transport costs.
Finds the vector {w, L, P, r} satisfying:
  - Goods market clearing / balanced trade (eq 11)
  - Price index (eq 12)
  - Land market clearing (eq 13)
  - Residential choice / labor allocation (eq 14)

Uses nested fixed-point iteration following the Redding (2016) pattern:
  Outer loop: update population L (residential choice)
  Inner loop: update wages w (balanced trade)

Returns NamedTuple with equilibrium quantities.
"""
function solve_equilibrium(A, H, T, sigma, alpha, epsilon, F, Lbar;
                           maxiter=2000, tol=1e-6, verbose=true)
    N = length(A)
    markup = sigma / (sigma - 1)  # constant markup (eq 2)

    # Initialize
    L = fill(Lbar / N, N)       # uniform population
    w = ones(N)                  # uniform wages

    # Precompute T^(1-σ)
    T_1ms = T .^ (1 - sigma)    # NxN matrix

    converged_L = false
    tradesh = zeros(N, N)

    if verbose
        println(">>>> Start Spatial Equilibrium Solver <<<<")
    end

    for outer in 1:maxiter
        # Number of varieties in each location (eq 34): M_i = L_i / (σF)
        M = L ./ (sigma * F)

        # ---- Inner loop: wage convergence (balanced trade, eq 11) ----
        converged_w = false
        for inner in 1:maxiter
            # Delivered price: p_{n,i} = markup * (w_i / A_i) * T_{n,i}  (eq 2 + transport)
            # p_{n,i}^{1-σ} = markup^{1-σ} * (w_i/A_i)^{1-σ} * T_{n,i}^{1-σ}
            unit_cost_1ms = (markup .* w ./ A) .^ (1 - sigma)  # N-vector: (markup * w_i/A_i)^{1-σ}

            # Numerator of trade shares: M_i * p_{n,i}^{1-σ} = (L_i/(σF)) * markup^{1-σ} * (w_i/A_i)^{1-σ} * T_{n,i}^{1-σ}
            # This is an NxN matrix where [n,i] = export from i to n
            nummat = T_1ms .* (ones(N) * (M .* unit_cost_1ms)')  # [n,i]

            # Price index: P_n^{1-σ} = Σ_i nummat[n,i]  (eq 12)
            P_1ms = sum(nummat, dims=2)[:]  # N-vector

            # Trade shares: π_{n,i} = nummat[n,i] / P_n^{1-σ}
            # π_{n,i} = share of n's tradable expenditure on goods from i
            tradesh = nummat ./ P_1ms  # NxN matrix [n,i]

            # Tradable expenditure in location n = α * v_n * L_n = w_n * L_n
            # (since v_n = w_n/α and tradable expenditure = α * v_n * L_n)
            income = w .* L  # labor income = tradable expenditure

            # Revenue of location n from selling to all destinations:
            # expend_n = Σ_i π_{i,n} * (w_i * L_i)  (income of n from all markets)
            expend = tradesh' * income  # [n] = Σ_i tradesh[i,n] * income[i]

            # Check convergence
            income_r = round.(income .* 1e6)
            expend_r = round.(expend .* 1e6)

            if income_r == expend_r
                converged_w = true
                break
            else
                # Update wages (excess demand raises wages)
                w_e = w .* (expend ./ income) .^ (1 / sigma)
                w = 0.25 .* w_e .+ 0.75 .* w
                w ./= geomean(w)  # normalize
            end
        end

        if !converged_w && verbose && outer > 20
            println("  WARNING: Wage loop did not converge at outer iteration $outer")
        end

        # ---- Outer loop: population convergence (residential choice, eq 14) ----

        # Price index: P_n = P_n^{1-σ}^{1/(1-σ)}
        P_1ms = sum(tradesh .* (ones(N) * P_1ms'), dims=2)  # recompute from tradesh
        # Actually, we already have P_1ms from the inner loop. Let's recompute properly.
        unit_cost_1ms = (markup .* w ./ A) .^ (1 - sigma)
        M = L ./ (sigma * F)
        nummat_final = T_1ms .* (ones(N) * (M .* unit_cost_1ms)')
        P_1ms = sum(nummat_final, dims=2)[:]
        P = P_1ms .^ (1 / (1 - sigma))

        # Land rent (eq 13): r_n = ((1-α)/α) * w_n * L_n / H_n
        r = ((1 - alpha) / alpha) .* w .* L ./ H

        # Total income per worker: v_n = w_n / α (land rents redistributed lump-sum)
        v = w ./ alpha

        # Indirect utility: U_n = v_n / (P_n^α * r_n^{1-α})  (eq 5)
        U = v ./ (P .^ alpha .* r .^ (1 - alpha))

        # Population shares (eq 14): L_n/L = U_n^ε / Σ_k U_k^ε
        U_eps = U .^ epsilon
        L_e = (U_eps ./ sum(U_eps)) .* Lbar

        # Check convergence
        L_r = round.(L .* 1e6)
        L_e_r = round.(L_e .* 1e6)

        if L_r == L_e_r
            converged_L = true
            if verbose
                println(">>>> Spatial Equilibrium Converged (outer=$outer) <<<<")
            end
            break
        else
            # Damped update
            L_update = L .* (L_e ./ L) .^ (1 / (epsilon * (1 - alpha)))
            L = 0.25 .* L_update .+ 0.75 .* L
        end
    end

    if !converged_L && verbose
        println("  WARNING: Population loop did not converge")
    end

    # Final computation of equilibrium quantities
    M = L ./ (sigma * F)
    unit_cost_1ms = (markup .* w ./ A) .^ (1 - sigma)
    nummat_final = T_1ms .* (ones(N) * (M .* unit_cost_1ms)')
    P_1ms = sum(nummat_final, dims=2)[:]
    P = P_1ms .^ (1 / (1 - sigma))
    r = ((1 - alpha) / alpha) .* w .* L ./ H
    v = w ./ alpha
    U = v ./ (P .^ alpha .* r .^ (1 - alpha))
    tradesh = nummat_final ./ P_1ms

    # Trade flows: X_{n,i} = π_{n,i} * w_n * L_n  (expenditure of n on goods from i)
    tradeflows = tradesh .* (w .* L)  # NxN matrix

    # Expected utility / aggregate welfare (eq 6):
    # Ũ = δ * [Σ_n (v_n / P_n^α * r_n^{1-α})^ε]^{1/ε}
    U_tilde = delta_frechet * sum(U .^ epsilon) ^ (1 / epsilon)

    # Domestic trade share (share of own goods in expenditure)
    dom_share = diag(tradesh)

    return (; w, L, P, r, v, U, U_tilde, tradesh, tradeflows, M, dom_share,
              converged=converged_L)
end

# ============================================================
# Productivity Inversion (Section IV.A.2)
# ============================================================

"""
    invert_productivity(L_target, H, T, sigma, alpha, epsilon, F, Lbar;
                        maxiter=200, tol=1e-4, verbose=true)

Given observed population distribution L_target, find the productivity vector A
that rationalizes it in equilibrium. Uses iterative updating:
    A_new = A .* (L_target ./ L_model)^step

This is the calibration procedure from Section IV.A.2 of the paper.
"""
function invert_productivity(L_target, H, T, sigma, alpha, epsilon, F, Lbar;
                             maxiter=200, tol=1e-4, verbose=true)
    N = length(L_target)
    A = ones(N)

    if verbose
        println(">>>> Inverting Productivity Parameters <<<<")
    end

    for iter in 1:maxiter
        result = solve_equilibrium(A, H, T, sigma, alpha, epsilon, F, Lbar;
                                   maxiter=2000, verbose=false)
        L_model = result.L

        # Check convergence
        max_dev = maximum(abs.(L_model .- L_target) ./ L_target)
        if verbose && iter % 20 == 0
            println("  Inversion iter $iter: max relative deviation = $(round(max_dev, digits=6))")
        end

        if max_dev < tol
            if verbose
                println(">>>> Productivity Inversion Converged (iter=$iter) <<<<")
            end
            return A
        end

        # Update productivity: raise A where L_model < L_target (need to attract more workers)
        # The elasticity from eq 43 guides the step size
        step = alpha / (1 / epsilon + (1 - alpha) - alpha / (sigma - 1))
        A = A .* (L_target ./ L_model) .^ step
        A ./= geomean(A)  # normalize
    end

    if verbose
        println("  WARNING: Productivity inversion did not converge")
    end
    return A
end

# ============================================================
# Government Infrastructure Optimization (Section III.G-H)
# ============================================================

"""
    theta_to_phi(theta, Z_budget, c_cost, N)

Convert unconstrained parameter vector θ to infrastructure investment vector φ
using softmax reparametrization that automatically satisfies:
  - φ_n ≥ 1 for all n (minimum infrastructure level)
  - Σ c_n * φ_n = Z (budget constraint binds)

The excess budget (Z - N*c) is allocated proportional to softmax(θ).
"""
function theta_to_phi(theta, Z_budget, c_cost, N)
    excess_budget = Z_budget - N * c_cost  # budget available beyond minimum φ=1
    if excess_budget <= 0
        return ones(N)  # no budget for improvements
    end
    # Softmax shares of excess budget
    theta_centered = theta .- maximum(theta)  # numerical stability
    exp_theta = exp.(theta_centered)
    shares = exp_theta ./ sum(exp_theta)
    # φ_n = 1 + s_n * (Z - N*c) / c
    phi = 1.0 .+ shares .* (excess_budget / c_cost)
    return phi
end

"""
    welfare_objective(theta, graph, edge_distances, A, H, sigma, alpha, epsilon,
                      gamma_infra, F, Lbar, Z_budget, c_cost, N)

Compute negative aggregate welfare for given infrastructure allocation θ.
Used as objective function for the government's optimization problem.
    theta → φ (via softmax) → edge weights → transport costs → equilibrium → -Ũ
"""
function welfare_objective(theta, graph, edge_distances, A, H, sigma, alpha, epsilon,
                           gamma_infra, F, Lbar, Z_budget, c_cost, N)
    phi = theta_to_phi(theta, Z_budget, c_cost, N)
    edge_wts = compute_edge_weights(graph, edge_distances, phi, gamma_infra)
    T = compute_transport_costs(graph, edge_wts, N)
    result = solve_equilibrium(A, H, T, sigma, alpha, epsilon, F, Lbar;
                               maxiter=1000, verbose=false)
    return -result.U_tilde  # minimize negative welfare
end

"""
    optimize_infrastructure(graph, edge_distances, A, H, sigma, alpha, epsilon,
                           gamma_infra, F, Lbar, Z_budget, c_cost;
                           phi_init=nothing, maxiter=200, verbose=true)

Solve the government's problem (eq 24):
    max_{φ} Ũ(φ)   s.t.  Σ c_n φ_n ≤ Z,  φ_n ≥ 1 ∀n

Uses Nelder-Mead optimization with softmax reparametrization to handle constraints.
The paper's solution algorithm (Section III.H) is:
  1. Start from initial investment vector
  2. Solve spatial equilibrium for given infrastructure
  3. Take utility-maximizing step
  4. Repeat until convergence

Returns NamedTuple with optimal phi, equilibrium, and welfare.
"""
function optimize_infrastructure(graph, edge_distances, A, H, sigma, alpha, epsilon,
                                 gamma_infra, F, Lbar, Z_budget, c_cost;
                                 phi_init=nothing, maxiter=200, verbose=true)
    N = length(A)

    if verbose
        println(">>>> Starting Infrastructure Optimization <<<<")
    end

    # Initialize theta (unconstrained parameters)
    if phi_init === nothing
        theta0 = zeros(N)  # uniform allocation
    else
        # Invert softmax: theta ∝ log(phi - 1) approximately
        excess = Z_budget - N * c_cost
        shares = (phi_init .- 1.0) .* c_cost ./ excess
        shares = max.(shares, 1e-10)  # avoid log(0)
        theta0 = log.(shares)
    end

    # Objective function
    obj = theta -> welfare_objective(theta, graph, edge_distances, A, H, sigma, alpha, epsilon,
                                     gamma_infra, F, Lbar, Z_budget, c_cost, N)

    # Optimize using Nelder-Mead (derivative-free, robust for noisy objectives)
    result = Optim.optimize(obj, theta0, NelderMead(),
                            Optim.Options(iterations=maxiter, show_trace=verbose,
                                         show_every=20, g_tol=1e-8))

    theta_opt = Optim.minimizer(result)
    phi_opt = theta_to_phi(theta_opt, Z_budget, c_cost, N)

    if verbose
        println(">>>> Infrastructure Optimization $(Optim.converged(result) ? "Converged" : "Completed") <<<<")
        println("  Iterations: $(Optim.iterations(result))")
        println("  Budget used: $(round(sum(c_cost .* phi_opt), digits=2)) / $Z_budget")
    end

    # Solve final equilibrium at optimum
    edge_wts = compute_edge_weights(graph, edge_distances, phi_opt, gamma_infra)
    T_opt = compute_transport_costs(graph, edge_wts, N)
    eq_opt = solve_equilibrium(A, H, T_opt, sigma, alpha, epsilon, F, Lbar;
                               verbose=false)

    return (; phi=phi_opt, T=T_opt, eq=eq_opt, optim_result=result)
end

"""
    optimize_infrastructure_constrained(graph, edge_distances, A, H, sigma, alpha, epsilon,
                                        gamma_infra, F, Lbar, Z_budget, c_cost, phi_lower;
                                        maxiter=200, verbose=true)

Solve the constrained government's problem where some infrastructure is pre-existing:
    max_{φ} Ũ(φ)   s.t.  Σ c_n φ_n ≤ Z,  φ_n ≥ φ_n^{lower} ∀n

This models path-dependence: the government cannot reduce existing infrastructure (eq 26).
The excess budget above the lower bound is allocated via softmax.
"""
function optimize_infrastructure_constrained(graph, edge_distances, A, H, sigma, alpha, epsilon,
                                              gamma_infra, F, Lbar, Z_budget, c_cost, phi_lower;
                                              maxiter=200, verbose=true)
    N = length(A)
    lower_cost = sum(c_cost .* phi_lower)

    if lower_cost >= Z_budget
        if verbose
            println("  Budget fully consumed by pre-existing infrastructure")
        end
        edge_wts = compute_edge_weights(graph, edge_distances, phi_lower, gamma_infra)
        T = compute_transport_costs(graph, edge_wts, N)
        eq = solve_equilibrium(A, H, T, sigma, alpha, epsilon, F, Lbar; verbose=false)
        return (; phi=phi_lower, T, eq, optim_result=nothing)
    end

    # Effective budget for additional investment
    Z_excess = Z_budget - lower_cost

    if verbose
        println(">>>> Starting Constrained Infrastructure Optimization <<<<")
        println("  Pre-existing cost: $(round(lower_cost, digits=2)), excess budget: $(round(Z_excess, digits=2))")
    end

    theta0 = zeros(N)

    function obj_constrained(theta)
        # Additional investment on top of lower bound
        theta_c = theta .- maximum(theta)
        exp_t = exp.(theta_c)
        shares = exp_t ./ sum(exp_t)
        phi_add = shares .* (Z_excess / c_cost)
        phi = phi_lower .+ phi_add

        edge_wts = compute_edge_weights(graph, edge_distances, phi, gamma_infra)
        T = compute_transport_costs(graph, edge_wts, N)
        result = solve_equilibrium(A, H, T, sigma, alpha, epsilon, F, Lbar;
                                   maxiter=1000, verbose=false)
        return -result.U_tilde
    end

    result = Optim.optimize(obj_constrained, theta0, NelderMead(),
                            Optim.Options(iterations=maxiter, show_trace=verbose,
                                         show_every=20))

    theta_opt = Optim.minimizer(result)
    theta_c = theta_opt .- maximum(theta_opt)
    exp_t = exp.(theta_c)
    shares = exp_t ./ sum(exp_t)
    phi_opt = phi_lower .+ shares .* (Z_excess / c_cost)

    if verbose
        println(">>>> Constrained Optimization $(Optim.converged(result) ? "Converged" : "Completed") <<<<")
    end

    edge_wts = compute_edge_weights(graph, edge_distances, phi_opt, gamma_infra)
    T_opt = compute_transport_costs(graph, edge_wts, N)
    eq_opt = solve_equilibrium(A, H, T_opt, sigma, alpha, epsilon, F, Lbar; verbose=false)

    return (; phi=phi_opt, T=T_opt, eq=eq_opt, optim_result=result)
end

# ============================================================
# Market Access and Centrality (Section IV.B)
# ============================================================

"""
    compute_market_access(T, w, L, sigma)

Compute market access for each location, following Donaldson and Hornbeck (2016):
    MA_n = Σ_i T_{n,i}^{1-σ} * E_i / MA_i

where E_i = w_i * L_i is expenditure. Solved iteratively.
"""
function compute_market_access(T, w, L, sigma)
    N = length(w)
    T_1ms = T .^ (1 - sigma)
    E = w .* L  # expenditure

    MA = ones(N)
    for iter in 1:500
        MA_new = T_1ms * (E ./ MA)
        if maximum(abs.(MA_new .- MA) ./ MA) < 1e-8
            break
        end
        MA = MA_new
    end
    return MA
end

"""
    compute_betweenness_centrality(graph, edge_weights, N)

Compute weighted betweenness centrality for each location.
Centrality measures how many shortest paths transit through a location,
which determines infrastructure investment incentives (eq 23).
"""
function compute_betweenness_centrality(graph, edge_weights, N)
    # Use path information to count transits
    n_vertices = nv(graph)
    wt_matrix = zeros(n_vertices, n_vertices)
    for e in edges(graph)
        i, j = src(e), dst(e)
        w = edge_weights[e]
        wt_matrix[i, j] = w
        wt_matrix[j, i] = w
    end

    centrality = zeros(N)
    for n in 1:N
        ds = dijkstra_shortest_paths(graph, n, wt_matrix)
        for i in 1:N
            if i == n
                continue
            end
            # Trace back the shortest path from i to n
            current = i
            while current != n && current != 0
                parent = ds.parents[current]
                if parent != n && parent != i && parent != 0
                    centrality[parent] += 1.0
                end
                current = parent
            end
        end
    end
    # Normalize
    centrality ./= (N * (N - 1))
    return centrality
end

# ============================================================
# Counterfactual Analysis (Section V)
# ============================================================

"""
    compare_welfare(phi1, phi2, graph, edge_distances, A, H, params...)

Compare aggregate welfare under two infrastructure allocations.
Returns percentage welfare change: ΔW = (Ũ₂ - Ũ₁) / Ũ₁ * 100
"""
function compare_welfare(phi1, phi2, graph, edge_distances, A, H,
                         sigma, alpha, epsilon, gamma_infra, F, Lbar)
    ew1 = compute_edge_weights(graph, edge_distances, phi1, gamma_infra)
    T1 = compute_transport_costs(graph, ew1, length(A))
    eq1 = solve_equilibrium(A, H, T1, sigma, alpha, epsilon, F, Lbar; verbose=false)

    ew2 = compute_edge_weights(graph, edge_distances, phi2, gamma_infra)
    T2 = compute_transport_costs(graph, ew2, length(A))
    eq2 = solve_equilibrium(A, H, T2, sigma, alpha, epsilon, F, Lbar; verbose=false)

    dW = (eq2.U_tilde - eq1.U_tilde) / eq1.U_tilde * 100
    return (; eq1, eq2, dW)
end

"""
    simulate_division(graph, edge_distances, A, H, coords, N_side,
                      sigma, alpha, epsilon, gamma_infra, F, Lbar, Z_budget, c_cost;
                      removed_rows=1)

Simulate a "division" shock by removing locations from the top row(s) of the grid.
This is analogous to the division of Germany removing East German districts.
Re-solves for optimal infrastructure in the remaining locations.

Returns results for the full economy and the post-division economy.
"""
function simulate_division(graph, edge_distances, A, H, coords, N_side,
                           sigma, alpha, epsilon, gamma_infra, F, Lbar, Z_budget, c_cost;
                           removed_rows=1)
    N = N_side^2
    # Locations to keep (remove top row(s))
    kept = Int[]
    for i in 1:N
        row = div(i - 1, N_side) + 1
        if row > removed_rows
            push!(kept, i)
        end
    end
    N_kept = length(kept)

    println("\n>>>> Simulating Division: removing $(N - N_kept) locations <<<<")

    # Build sub-network for remaining locations
    # Create mapping from old indices to new indices
    old_to_new = Dict{Int, Int}()
    for (new_idx, old_idx) in enumerate(kept)
        old_to_new[old_idx] = new_idx
    end

    # Build new graph with only kept locations
    g_new = SimpleGraph(N_kept)
    edge_dist_new = Dict{Edge{Int}, Float64}()
    for e in edges(graph)
        i, j = src(e), dst(e)
        if haskey(old_to_new, i) && haskey(old_to_new, j)
            ni, nj = old_to_new[i], old_to_new[j]
            add_edge!(g_new, ni, nj)
            edge_dist_new[Edge(ni, nj)] = edge_distances[e]
            edge_dist_new[Edge(nj, ni)] = edge_distances[e]
        end
    end

    A_new = A[kept]
    H_new = H[kept]
    coords_new = coords[kept, :]
    Lbar_new = Lbar * N_kept / N  # proportional labor force

    # Adjust budget: same total budget for post-division economy
    # (the government reallocates the full budget to remaining locations)
    Z_new = Z_budget

    return (; graph=g_new, edge_distances=edge_dist_new, A=A_new, H=H_new,
              coords=coords_new, kept, N=N_kept, Lbar=Lbar_new, Z=Z_new)
end

# ============================================================
# Main Execution
# ============================================================

println("\n" * "="^60)
println("Santamaria (JEEA 2026): Spatial Model with Endogenous Infrastructure")
println("="^60)

# ---- 1. Build network and generate geography ----

Random.seed!(1)

net = build_grid_network(N_side)
println("\nNetwork: $(N) locations on $(N_side)x$(N_side) grid, $(ne(net.graph)) edges")

# Land endowments: proportional to area (uniform in this grid)
H = ones(N)

# Productivities: log-normal, geo-mean normalized
A = exp.(0.3 .* randn(N))
A ./= geomean(A)

println("\nProductivities (A):")
@show round.(A, digits=3)

# ---- 2. Baseline equilibrium with NO infrastructure (φ = 1) ----

println("\n" * "="^60)
println("CASE 1: Baseline equilibrium (no infrastructure, φ = 1)")
println("="^60)

phi_base = ones(N)
ew_base = compute_edge_weights(net.graph, net.edge_distances, phi_base, gamma_infra)
T_base = compute_transport_costs(net.graph, ew_base, N)

eq_base = solve_equilibrium(A, H, T_base, sigma, alpha, epsilon, F, Lbar)

println("\nBaseline results:")
println("  Welfare (Ũ): $(round(eq_base.U_tilde, digits=4))")
println("  Max trade balance deviation: $(round(maximum(abs.(sum(eq_base.tradeflows, dims=2)[:] .- eq_base.w .* eq_base.L)), digits=8))")
println("  Population range: [$(round(minimum(eq_base.L), digits=2)), $(round(maximum(eq_base.L), digits=2))]")
println("  Wage range: [$(round(minimum(eq_base.w), digits=4)), $(round(maximum(eq_base.w), digits=4))]")
println("  Domestic trade share range: [$(round(minimum(eq_base.dom_share), digits=4)), $(round(maximum(eq_base.dom_share), digits=4))]")

# ---- 3. Uniform infrastructure ----

println("\n" * "="^60)
println("CASE 2: Uniform infrastructure (φ = Z/N for all)")
println("="^60)

phi_uniform = fill(Z_budget / N, N)
ew_uniform = compute_edge_weights(net.graph, net.edge_distances, phi_uniform, gamma_infra)
T_uniform = compute_transport_costs(net.graph, ew_uniform, N)

eq_uniform = solve_equilibrium(A, H, T_uniform, sigma, alpha, epsilon, F, Lbar)

dW_uniform = (eq_uniform.U_tilde - eq_base.U_tilde) / eq_base.U_tilde * 100
println("\nUniform infrastructure results:")
println("  Welfare (Ũ): $(round(eq_uniform.U_tilde, digits=4))")
println("  Welfare gain vs no infrastructure: $(round(dW_uniform, digits=2))%")

# ---- 4. Optimal infrastructure ----

println("\n" * "="^60)
println("CASE 3: Optimal infrastructure (government problem)")
println("="^60)

opt_result = optimize_infrastructure(net.graph, net.edge_distances, A, H,
                                     sigma, alpha, epsilon, gamma_infra, F, Lbar,
                                     Z_budget, c_cost;
                                     maxiter=300, verbose=true)

dW_opt_vs_base = (opt_result.eq.U_tilde - eq_base.U_tilde) / eq_base.U_tilde * 100
dW_opt_vs_uniform = (opt_result.eq.U_tilde - eq_uniform.U_tilde) / eq_uniform.U_tilde * 100

println("\nOptimal infrastructure results:")
println("  Welfare (Ũ): $(round(opt_result.eq.U_tilde, digits=4))")
println("  Welfare gain vs no infrastructure: $(round(dW_opt_vs_base, digits=2))%")
println("  Welfare gain vs uniform: $(round(dW_opt_vs_uniform, digits=2))%")
println("  Infrastructure allocation (φ):")
@show round.(opt_result.phi, digits=3)
println("  Budget used: $(round(sum(c_cost .* opt_result.phi), digits=2)) / $Z_budget")

# ---- 5. Division shock and infrastructure reshaping ----

println("\n" * "="^60)
println("CASE 4: Division shock (remove top row of grid)")
println("="^60)

# First, solve optimal infrastructure for the full pre-division economy
# (already done in Case 3, this is opt_result)

# Pre-division optimal phi restricted to remaining locations
div_data = simulate_division(net.graph, net.edge_distances, A, H, net.coords, N_side,
                             sigma, alpha, epsilon, gamma_infra, F, Lbar, Z_budget, c_cost;
                             removed_rows=1)

# Solve post-division optimal (unconstrained)
opt_post_div = optimize_infrastructure(div_data.graph, div_data.edge_distances,
                                       div_data.A, div_data.H,
                                       sigma, alpha, epsilon, gamma_infra, F,
                                       div_data.Lbar, div_data.Z, c_cost;
                                       maxiter=300, verbose=true)

println("\nPost-division optimal infrastructure:")
@show round.(opt_post_div.phi, digits=3)

# Compare: pre-division optimal (restricted to kept locations) vs post-division optimal
phi_pre_kept = opt_result.phi[div_data.kept]

# Equilibrium under pre-division infrastructure in post-division world
ew_pre_kept = compute_edge_weights(div_data.graph, div_data.edge_distances, phi_pre_kept, gamma_infra)
T_pre_kept = compute_transport_costs(div_data.graph, ew_pre_kept, div_data.N)
eq_pre_kept = solve_equilibrium(div_data.A, div_data.H, T_pre_kept,
                                sigma, alpha, epsilon, F, div_data.Lbar; verbose=false)

dW_reshape = (opt_post_div.eq.U_tilde - eq_pre_kept.U_tilde) / eq_pre_kept.U_tilde * 100
println("\nGains from infrastructure reshaping after division: $(round(dW_reshape, digits=2))%")

# ---- 6. Cost of path dependence ----

println("\n" * "="^60)
println("CASE 5: Cost of path dependence")
println("="^60)

# The constrained problem: the government keeps pre-existing infrastructure
# and allocates additional budget on top (eq 26: φ^post ≥ φ^pre)
# We give the government additional budget = 50% of pre-existing investment
phi_lower = phi_pre_kept
pre_existing_cost = sum(c_cost .* phi_lower)
Z_constrained = pre_existing_cost + 0.5 * pre_existing_cost  # 50% additional budget
Z_unconstrained = Z_constrained  # same total budget, but free to allocate from scratch

println("  Pre-existing infrastructure cost: $(round(pre_existing_cost, digits=2))")
println("  Total budget (constrained): $(round(Z_constrained, digits=2))")

# Constrained optimization: pre-division infrastructure as lower bound, plus extra budget
opt_constrained = optimize_infrastructure_constrained(
    div_data.graph, div_data.edge_distances, div_data.A, div_data.H,
    sigma, alpha, epsilon, gamma_infra, F, div_data.Lbar, Z_constrained, c_cost, phi_lower;
    maxiter=300, verbose=true)

# Unconstrained optimization with same total budget but no lower bounds
opt_unconstrained = optimize_infrastructure(
    div_data.graph, div_data.edge_distances, div_data.A, div_data.H,
    sigma, alpha, epsilon, gamma_infra, F, div_data.Lbar, Z_unconstrained, c_cost;
    maxiter=300, verbose=true)

# Baseline: equilibrium under only the pre-existing infrastructure (no additional investment)
dW_constrained = (opt_constrained.eq.U_tilde - eq_pre_kept.U_tilde) / eq_pre_kept.U_tilde * 100
dW_unconstrained = (opt_unconstrained.eq.U_tilde - eq_pre_kept.U_tilde) / eq_pre_kept.U_tilde * 100
dW_path_dep = dW_unconstrained - dW_constrained

println("\nPath dependence results (relative to pre-existing infra only):")
println("  Constrained optimal welfare gain: $(round(dW_constrained, digits=2))%")
println("  Unconstrained optimal welfare gain: $(round(dW_unconstrained, digits=2))%")
println("  Cost of path dependence: $(round(dW_path_dep, digits=2)) percentage points")

# ---- 7. Market access and centrality analysis ----

println("\n" * "="^60)
println("Diagnostics: Market Access and Centrality")
println("="^60)

MA = compute_market_access(T_base, eq_base.w, eq_base.L, sigma)
centrality = compute_betweenness_centrality(net.graph, ew_base, N)

# Correlation between optimal infrastructure and exports/centrality (eq 23 intuition)
exports = sum(eq_base.tradeflows, dims=1)[:]  # total exports from each location
println("\nCorrelation between optimal φ and:")
println("  Total exports: $(round(cor(opt_result.phi, exports), digits=3))")
println("  Betweenness centrality: $(round(cor(opt_result.phi, centrality), digits=3))")
println("  Market access: $(round(cor(opt_result.phi, MA), digits=3))")

# Gravity check: trade decreases with distance
dists = [sqrt(sum((net.coords[i,:] .- net.coords[j,:]).^2)) for i in 1:N, j in 1:N]
mask = [i != j for i in 1:N, j in 1:N]
log_trade = log.(eq_base.tradeflows[mask])
log_dist = log.(dists[mask])
# Simple OLS: log(trade) = a + β * log(dist)
X_ols = hcat(ones(length(log_dist)), log_dist)
beta_gravity = (X_ols' * X_ols) \ (X_ols' * log_trade)
println("\nGravity equation: log(trade) = $(round(beta_gravity[1], digits=3)) + $(round(beta_gravity[2], digits=3)) * log(distance)")
println("  Distance elasticity: $(round(beta_gravity[2], digits=3)) (expected: negative, paper finds -2.8)")

# ============================================================
# Plots
# ============================================================

println("\n>>>> Generating Plots <<<<")

# Plot 1: Infrastructure heatmaps (optimal vs uniform)
phi_grid_opt = reshape(opt_result.phi, N_side, N_side)'  # reshape for heatmap (row=y, col=x)
phi_grid_uniform = reshape(phi_uniform, N_side, N_side)'

p1 = heatmap(1:N_side, 1:N_side, phi_grid_opt,
    title="Optimal Infrastructure (φ)", xlabel="x", ylabel="y",
    color=:YlOrRd, aspect_ratio=1, clims=(1, maximum(opt_result.phi)))

p2 = heatmap(1:N_side, 1:N_side, phi_grid_uniform,
    title="Uniform Infrastructure (φ)", xlabel="x", ylabel="y",
    color=:YlOrRd, aspect_ratio=1, clims=(1, maximum(opt_result.phi)))

# Plot 2: Population distribution (baseline vs optimal infrastructure)
L_grid_base = reshape(eq_base.L, N_side, N_side)'
L_grid_opt = reshape(opt_result.eq.L, N_side, N_side)'

p3 = heatmap(1:N_side, 1:N_side, L_grid_base,
    title="Population (no infrastructure)", xlabel="x", ylabel="y",
    color=:Blues, aspect_ratio=1)

p4 = heatmap(1:N_side, 1:N_side, L_grid_opt,
    title="Population (optimal infrastructure)", xlabel="x", ylabel="y",
    color=:Blues, aspect_ratio=1)

# Plot 3: Welfare comparison bar chart
welfare_levels = [eq_base.U_tilde, eq_uniform.U_tilde, opt_result.eq.U_tilde]
welfare_labels = ["No infra", "Uniform", "Optimal"]
p5 = bar(welfare_labels, welfare_levels, label="Welfare (Ũ)",
    title="Aggregate Welfare Comparison", ylabel="Expected Utility",
    color=[:gray, :steelblue, :darkgreen], alpha=0.8)

# Plot 4: Division - infrastructure reshaping
if length(div_data.kept) > 0
    N_post = div_data.N
    N_side_post = N_side
    N_rows_post = N_side - 1

    phi_pre_grid = fill(NaN, N_side, N_side)
    phi_post_grid = fill(NaN, N_side, N_side)
    for (new_idx, old_idx) in enumerate(div_data.kept)
        row = div(old_idx - 1, N_side) + 1
        col = mod(old_idx - 1, N_side) + 1
        phi_pre_grid[row, col] = phi_pre_kept[new_idx]
        phi_post_grid[row, col] = opt_post_div.phi[new_idx]
    end

    p6 = heatmap(1:N_side, 1:N_side, phi_pre_grid',
        title="Pre-division φ (kept)", xlabel="x", ylabel="y",
        color=:YlOrRd, aspect_ratio=1)

    p7 = heatmap(1:N_side, 1:N_side, phi_post_grid',
        title="Post-division optimal φ", xlabel="x", ylabel="y",
        color=:YlOrRd, aspect_ratio=1)
end

# Plot 5: Infrastructure vs exports and centrality (eq 23 intuition)
p8 = scatter(exports, opt_result.phi,
    xlabel="Total exports", ylabel="Optimal φ",
    title="Infrastructure vs Exports (eq 23)",
    label="", markersize=5, color=:steelblue)

p9 = scatter(centrality, opt_result.phi,
    xlabel="Betweenness centrality", ylabel="Optimal φ",
    title="Infrastructure vs Centrality (eq 23)",
    label="", markersize=5, color=:darkred)

# Combine and save
p_main = plot(p1, p2, p3, p4, layout=(2, 2), size=(900, 800))
savefig(p_main, joinpath(@__DIR__, "graphs", "infrastructure_equilibrium.pdf"))

p_welfare = plot(p5, layout=(1, 1), size=(500, 400))
savefig(p_welfare, joinpath(@__DIR__, "graphs", "welfare_comparison.pdf"))

p_division = plot(p6, p7, layout=(1, 2), size=(900, 400))
savefig(p_division, joinpath(@__DIR__, "graphs", "division_reshaping.pdf"))

p_mechanism = plot(p8, p9, layout=(1, 2), size=(900, 400))
savefig(p_mechanism, joinpath(@__DIR__, "graphs", "infrastructure_mechanism.pdf"))

println("\nPlots saved to $(joinpath(@__DIR__, "graphs"))/")

println("\n" * "="^60)
println("Done.")
println("="^60)
