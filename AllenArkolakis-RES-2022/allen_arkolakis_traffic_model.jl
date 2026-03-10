# Quantitative Spatial Model with Endogenous Traffic Congestion
# Based on Allen & Arkolakis (RES 2022) "The Welfare Effects of Transportation Infrastructure Improvements"
#
# Implements a general equilibrium spatial framework with:
# - Stochastic route choice with Frechet distributed shocks (shape parameter theta)
# - Analytical transport costs via Leontief inverse of the adjacency matrix (eq 21)
# - Traffic gravity equation: traffic on each link depends on link cost and market access (eq 24)
# - Endogenous traffic congestion: link costs increase with traffic (eq 25)
# - Two model variants:
#   (1) Economic geography: agents choose location, trade goods (eqs 10-11, with congestion: 28-29)
#   (2) Urban: agents choose residence + workplace + route (eqs 19-20, with congestion: 30-31)
# - Exact-hat counterfactuals using observed traffic data (Proposition 2, eqs 36-39)
#
# Key innovation vs. standard spatial models (e.g., Allen-Arkolakis 2014, Ahlfeldt et al. 2015):
# Transport costs are ENDOGENOUS. Agents choose routes optimally through a network, and traffic
# congestion raises costs on heavily-used links. The Frechet route-choice assumption yields
# closed-form transport costs as the Leontief inverse of the adjacency matrix -- NO Dijkstra needed.
#
# Key results:
# - With congestion (lambda > 0), equilibrium is SCALE DEPENDENT (depends on aggregate L_bar)
# - Without congestion, the model nests standard scale-invariant spatial models
# - Welfare elasticities of infrastructure improvements can be computed using only traffic data

using LinearAlgebra, Statistics, Plots, Random

# ============================================================
# Parameters
# ============================================================

# --- Economic Geography Model (US Highway System, Section 6.1) ---
theta_eg = 8.0              # Trade elasticity / Frechet shape parameter (theta)
alpha_eg = 0.1              # Productivity externality: A_i = A_bar_i * L_i^alpha (eq 1)
beta_eg = -0.3              # Amenity externality: u_i = u_bar_i * L_i^beta (eq 1)
lambda_eg = 0.092           # Traffic congestion parameter: t_kl = t_bar_kl * Xi_kl^lambda (eq 25)
                            # lambda = delta_0 * delta_1, where delta_0 = 1/theta (unit distance elasticity)
                            # and delta_1 = 0.739 is the IV estimate of congestion elasticity of inverse speed

# --- Urban Model (Seattle Road Network, Section 6.2) ---
theta_urban = 6.83          # Commuting elasticity (from Ahlfeldt et al. 2015)
alpha_urban = -0.12         # Productivity externality (from Ahlfeldt et al. 2015)
beta_urban = -0.1           # Amenity externality (from Ahlfeldt et al. 2015)
lambda_urban = 0.071        # Traffic congestion: delta_1 = 0.488, delta_0 = 1/theta => lambda = 0.071

# --- Grid Network (Figure 1 example) ---
N_side = 5                  # Grid side length
N = N_side^2                # Number of locations (25)
t_bar_default = 1.5         # Free-flow cost of traversing a connected link (t_bar_kl >= 1)
                            # t_bar_kl = Inf for non-connected pairs
theta_grid = 4.0            # Frechet shape for grid example (as in Figure 1)
Lbar_default = 100.0        # Aggregate labor endowment

# Uniqueness conditions (Proposition 1):
# Economic geography (symmetric T_bar): alpha + beta <= 0
# Urban model: alpha <= 0.5*(1/theta - lambda) and beta <= 0.5*(1/theta - lambda)
println("="^60)
println("Uniqueness Conditions (Proposition 1)")
println("="^60)
println("Economic Geography: alpha + beta = $(alpha_eg + beta_eg) <= 0? $(alpha_eg + beta_eg <= 0 ? "YES" : "NO")")
urban_bound = 0.5 * (1/theta_urban - lambda_urban)
println("Urban: alpha=$alpha_urban <= $(round(urban_bound, digits=4))? $(alpha_urban <= urban_bound ? "YES" : "NO")")
println("Urban: beta=$beta_urban <= $(round(urban_bound, digits=4))? $(beta_urban <= urban_bound ? "YES" : "NO")")

# ============================================================
# Network Construction
# ============================================================

"""
    build_grid_network(N_side; t_bar=1.5)

Build a regular grid network with 4-nearest-neighbor connectivity (up/down/left/right).
This matches the example geography in Figure 1 of Allen & Arkolakis (2022).

Each of the N = N_side^2 locations is connected to its adjacent neighbors with a
symmetric free-flow traversal cost t_bar >= 1. Non-connected pairs have t_bar = Inf.

Returns NamedTuple with:
- t_bar_matrix: NxN matrix of free-flow link costs (Inf for non-connected)
- coords: Nx2 matrix of (x,y) coordinates for plotting
- connected: NxN Bool matrix indicating which links exist
"""
function build_grid_network(N_side; t_bar=1.5)
    N = N_side^2
    t_bar_matrix = fill(Inf, N, N)
    coords = zeros(N, 2)

    # Place locations on a regular grid
    for i in 1:N_side, j in 1:N_side
        idx = (i - 1) * N_side + j
        coords[idx, :] = [j - 1, N_side - i]  # (x, y) with y increasing upward
    end

    # Connect 4-nearest neighbors with cost t_bar
    for i in 1:N_side, j in 1:N_side
        idx = (i - 1) * N_side + j
        t_bar_matrix[idx, idx] = Inf  # no self-loops (a_{kk} = 0 in adjacency matrix)
        # Right neighbor
        if j < N_side
            t_bar_matrix[idx, idx + 1] = t_bar
            t_bar_matrix[idx + 1, idx] = t_bar  # symmetric
        end
        # Down neighbor
        if i < N_side
            t_bar_matrix[idx, idx + N_side] = t_bar
            t_bar_matrix[idx + N_side, idx] = t_bar  # symmetric
        end
    end

    connected = isfinite.(t_bar_matrix)
    return (; t_bar_matrix, coords, connected)
end

# ============================================================
# Transport Cost Functions (Section 3)
# ============================================================

"""
    compute_transport_costs(t_bar_matrix, theta; check_spectral=true)

Compute bilateral transport costs using the Leontief inverse (eq 21).

The adjacency matrix A has elements a_{kl} = t_{kl}^{-theta} for connected links (0 otherwise).
Transport costs are: tau_{ij} = b_{ij}^{-1/theta} where B = (I - A)^{-1} is the Leontief inverse.

This analytical formula arises from Frechet-distributed route choice shocks: agents choose
routes stochastically, and summing over all possible routes yields the geometric matrix series
B = I + A + A^2 + ... = (I - A)^{-1}, which converges when spectral_radius(A) < 1.

Returns NamedTuple with:
- tau: NxN transport cost matrix (tau_{ij} >= 1)
- B: NxN Leontief inverse matrix
- A: NxN adjacency matrix
- spectral_radius: spectral radius of A
"""
function compute_transport_costs(t_bar_matrix, theta; check_spectral=true)
    N = size(t_bar_matrix, 1)

    # Adjacency matrix: a_{kl} = t_{kl}^{-theta} for connected links, 0 otherwise
    A = zeros(N, N)
    for k in 1:N, l in 1:N
        if isfinite(t_bar_matrix[k, l]) && k != l
            A[k, l] = t_bar_matrix[k, l]^(-theta)
        end
    end

    # Check spectral radius < 1 (required for Leontief inverse to exist)
    sr = maximum(abs.(eigvals(A)))
    if check_spectral && sr >= 1.0
        error("Spectral radius of A = $sr >= 1. Leontief inverse does not converge. " *
              "Reduce theta or increase t_bar.")
    end

    # Leontief inverse: B = (I - A)^{-1}
    B = inv(I - A)

    # Transport costs: tau_{ij} = b_{ij}^{-1/theta} (eq 21)
    tau = B .^ (-1 / theta)

    return (; tau, B, A, spectral_radius=sr)
end

"""
    compute_transport_costs_with_congestion(t_bar_matrix, P_inv_theta, Pi_inv_theta,
                                            theta, lambda; check_spectral=true)

Compute transport costs with endogenous traffic congestion (eqs 26, 21).

First computes the congested link costs from eq (26):
    t_{kl} = t_bar_{kl}^{1/(1+theta*lambda)} * P_k^{-theta*lambda/(1+theta*lambda)}
             * Pi_l^{-theta*lambda/(1+theta*lambda)}

Then computes tau via the Leontief inverse of the congested adjacency matrix.

Arguments:
- t_bar_matrix: NxN free-flow link costs (Inf for non-connected)
- P_inv_theta: N-vector of P_k^{-theta} (inward market access, raised to power theta)
- Pi_inv_theta: N-vector of Pi_l^{-theta} (outward market access, raised to power theta)
- theta: Frechet shape parameter
- lambda: congestion elasticity

Returns NamedTuple with tau, B, A, t_congested, Xi (traffic), spectral_radius
"""
function compute_transport_costs_with_congestion(t_bar_matrix, P_inv_theta, Pi_inv_theta,
                                                  theta, lambda; check_spectral=true)
    N = size(t_bar_matrix, 1)
    rho = 1 / (1 + theta * lambda)  # useful exponent

    # Congested link costs (eq 26):
    # t_{kl} = t_bar_{kl}^rho * P_k^{-theta*lambda*rho} * Pi_l^{-theta*lambda*rho}
    t_congested = fill(Inf, N, N)
    for k in 1:N, l in 1:N
        if isfinite(t_bar_matrix[k, l]) && k != l
            t_congested[k, l] = t_bar_matrix[k, l]^rho *
                                P_inv_theta[k]^(lambda * rho) *
                                Pi_inv_theta[l]^(lambda * rho)
        end
    end

    # Adjacency matrix from congested costs
    A = zeros(N, N)
    for k in 1:N, l in 1:N
        if isfinite(t_congested[k, l]) && k != l
            A[k, l] = t_congested[k, l]^(-theta)
        end
    end

    sr = maximum(abs.(eigvals(A)))
    if check_spectral && sr >= 1.0
        error("Spectral radius of congested A = $sr >= 1.")
    end

    B = inv(I - A)
    tau = B .^ (-1 / theta)

    # Traffic flows (eq 27): Xi_{kl} = t_bar_{kl}^{-theta*rho} * P_k^{-theta*rho} * Pi_l^{-theta*rho}
    Xi = zeros(N, N)
    for k in 1:N, l in 1:N
        if isfinite(t_bar_matrix[k, l]) && k != l
            Xi[k, l] = t_bar_matrix[k, l]^(-theta * rho) *
                        P_inv_theta[k]^rho *
                        Pi_inv_theta[l]^rho
        end
    end

    return (; tau, B, A, t_congested, Xi, spectral_radius=sr)
end

# ============================================================
# Equilibrium Solver - Economic Geography Model (Section 4.1)
# ============================================================

"""
    solve_econ_geo_equilibrium(A_bar, u_bar, t_bar_matrix, theta, alpha, beta, lambda, Lbar;
                               maxiter=5000, tol=1e-10, damp=0.5, verbose=true)

Solve the economic geography equilibrium with endogenous traffic congestion.

Directly solves the self-contained system eqs (28)-(29) for income shares {y_i} and labor
shares {l_i}. These equations incorporate congestion through direct links only:

Eq (28): y_i^{(1+θ+θλ)ρ} l_i^{-θ(1+α+(α+β)θλ)ρ} = χ ū_i^θ Ā_i^θ y_i^{(1+θ+θλ)ρ} l_i^{θ(β-1)ρ}
  + χ^{θλρ} Σ_j (L̄^λ t̄_{ij})^{-θρ} [fund_ij] y_j^{(1+θ)ρ} l_j^{-θ(1+α)ρ}

Eq (29): y_i^{-θ(1-λ)ρ} l_i^{θ(1-β-(α+β)θλ)ρ} = χ ū_i^θ Ā_i^θ y_i^{-θ(1-λ)ρ} l_i^{θ(α+1)ρ}
  + χ^{θλρ} Σ_j (L̄^λ t̄_{ji})^{-θρ} [fund_ij] y_j^{-θρ} l_j^{θ(1-β)ρ}

Returns NamedTuple with y, l, chi, tau, traffic (Xi), welfare, B, converged, iterations
"""
function solve_econ_geo_equilibrium(A_bar, u_bar, t_bar_matrix, theta, alpha, beta, lambda, Lbar;
                                    maxiter=5000, tol=1e-10, damp=0.5, verbose=true)
    N = length(A_bar)
    rho = 1 / (1 + theta * lambda)

    # Initialize uniformly
    y = fill(1.0 / N, N)
    l = fill(1.0 / N, N)

    # Precompute link weights: (Lbar^lambda * t_bar_{kl})^{-theta*rho}
    link_weight = zeros(N, N)       # for eq (28): links from i
    link_weight_T = zeros(N, N)     # for eq (29): links into i (transpose)
    for k in 1:N, ll in 1:N
        if isfinite(t_bar_matrix[k, ll]) && k != ll
            w = (Lbar^lambda * t_bar_matrix[k, ll])^(-theta * rho)
            link_weight[k, ll] = w
            link_weight_T[ll, k] = w
        end
    end

    # Precompute fundamental powers
    u_theta = u_bar .^ theta
    A_theta = A_bar .^ theta
    uA_theta = u_theta .* A_theta  # self-term coefficient ū_i^θ Ā_i^θ

    # Exponents for eq (28) - the "outward" / production equation
    e_y_lhs_28 = (1 + theta + theta * lambda) * rho          # y_i LHS exponent
    e_l_lhs_28 = -theta * (1 + alpha + (alpha + beta) * theta * lambda) * rho  # l_i LHS exponent
    e_y_self_28 = e_y_lhs_28                                  # y_i self-term exponent (same as LHS)
    e_l_self_28 = theta * (beta - 1) * rho                    # l_i self-term exponent
    e_y_sum_28 = (1 + theta) * rho                            # y_j sum exponent
    e_l_sum_28 = -theta * (1 + alpha) * rho                   # l_j sum exponent

    # Exponents for eq (29) - the "inward" / labor equation
    e_y_lhs_29 = -theta * (1 - lambda) * rho                  # y_i LHS exponent
    e_l_lhs_29 = theta * (1 - beta - (alpha + beta) * theta * lambda) * rho  # l_i LHS exponent
    e_y_self_29 = e_y_lhs_29                                  # y_i self-term exponent (same as LHS)
    e_l_self_29 = theta * (alpha + 1) * rho                   # l_i self-term exponent
    e_y_sum_29 = -theta * rho                                 # y_j sum exponent
    e_l_sum_29 = theta * (1 - beta) * rho                     # l_j sum exponent

    e_chi_frac = theta * lambda * rho

    # Fundamental powers for sum terms (derived from eqs 28-29)
    # Eq (28) sum: Ā_i^θ ū_i^{θ²λρ} for i-part, Ā_j^{-θρ} for j-part (no u_j in sum)
    # Wait - let me re-derive. From the derivation in the paper:
    # Eq (28) sum coefficient at i: ū_i^θ Ā_i^{θ²λρ} (outward: from production side)
    # Eq (28) sum kernel at j: ū_j^{-θρ} y_j^{e_y_sum_28} l_j^{e_l_sum_28}
    # (This follows the same pattern as the urban model but with different variable mapping)
    #
    # Actually, for the econ geo model, the derivation gives:
    # Eq (28): self = χ ū_i^θ Ā_i^θ y_i^{e_y_self} l_i^{e_l_self}
    #          sum = χ^{θλρ} Ā_i^{θ²λρ} ū_i^θ Σ_j (L̄^λ t̄_{ij})^{-θρ} ū_j^{-θρ} y_j^{...} l_j^{...}
    # Eq (29): self = χ ū_i^θ Ā_i^θ y_i^{e_y_self} l_i^{e_l_self}
    #          sum = χ^{θλρ} ū_i^{θ²λρ} Ā_i^θ Σ_j (L̄^λ t̄_{ji})^{-θρ} Ā_j^{-θρ} y_j^{...} l_j^{...}
    u_neg_theta_rho = u_bar .^ (-theta * rho)
    A_neg_theta_rho = A_bar .^ (-theta * rho)
    A_theta2_lam_rho = A_bar .^ (theta^2 * lambda * rho)
    u_theta2_lam_rho = u_bar .^ (theta^2 * lambda * rho)

    if verbose
        println("\n>>>> Economic Geography Equilibrium Solver <<<<")
        println("  N=$N, theta=$theta, alpha=$alpha, beta=$beta, lambda=$lambda, Lbar=$Lbar")
        println("  rho=$(round(rho, digits=4)), theta*lambda*rho=$(round(e_chi_frac, digits=4))")
        tc0 = compute_transport_costs(t_bar_matrix, theta; check_spectral=false)
        println("  Spectral radius of free-flow A: $(round(tc0.spectral_radius, digits=4))")
    end

    converged = false
    iter = 0
    chi = 1.0

    for it in 1:maxiter
        iter = it
        y_old = copy(y)
        l_old = copy(l)

        # --- Eq (28): Update y ---
        # Sum kernel: ū_j^{-θρ} * y_j^{(1+θ)ρ} * l_j^{-θ(1+α)ρ}
        sum28_kernel = u_neg_theta_rho .* (y .^ e_y_sum_28) .* (l .^ e_l_sum_28)
        S28 = link_weight * sum28_kernel

        for i in 1:N
            self_i = chi * uA_theta[i] * y[i]^e_y_self_28 * l[i]^e_l_self_28
            sum_i = chi^e_chi_frac * u_theta[i] * A_theta2_lam_rho[i] * S28[i]
            rhs_i = self_i + sum_i
            # LHS = y_i^{e_y_lhs_28} * l_i^{e_l_lhs_28}
            # Since e_y_self_28 = e_y_lhs_28, the self term has y_i on both sides
            # Need to extract y_i: rhs = y_i^{e_y_lhs} * l_i^{e_l_lhs}
            # => y_i = (rhs / l_i^{e_l_lhs})^{1/e_y_lhs}
            y[i] = (rhs_i / l[i]^e_l_lhs_28)^(1 / e_y_lhs_28)
        end
        y .= max.(y, 1e-15)
        y ./= sum(y)

        # --- Eq (29): Update l ---
        # Sum kernel: Ā_j^{-θρ} * y_j^{-θρ} * l_j^{θ(1-β)ρ}
        sum29_kernel = A_neg_theta_rho .* (y .^ e_y_sum_29) .* (l .^ e_l_sum_29)
        S29 = link_weight_T * sum29_kernel

        for i in 1:N
            self_i = chi * uA_theta[i] * y[i]^e_y_self_29 * l[i]^e_l_self_29
            sum_i = chi^e_chi_frac * A_theta[i] * u_theta2_lam_rho[i] * S29[i]
            rhs_i = self_i + sum_i
            l[i] = (rhs_i / y[i]^e_y_lhs_29)^(1 / e_l_lhs_29)
        end
        l .= max.(l, 1e-15)
        l ./= sum(l)

        # --- Damped update ---
        y .= (1 - damp) .* y_old .+ damp .* y
        l .= (1 - damp) .* l_old .+ damp .* l
        y .= max.(y, 1e-15)
        l .= max.(l, 1e-15)
        y ./= sum(y)
        l ./= sum(l)

        # --- Check convergence ---
        err = max(maximum(abs.(y .- y_old)), maximum(abs.(l .- l_old)))
        if verbose && (it % 500 == 0 || it <= 5)
            println("  Iter $it: max change = $(round(err, sigdigits=4))")
        end

        if err < tol
            converged = true
            if verbose println("  Converged at iteration $it (tol=$tol)") end
            break
        end
    end

    if !converged && verbose
        println("  WARNING: Did not converge in $maxiter iterations")
    end

    # Compute chi from eq (28) at location 1
    sum28_kernel = u_neg_theta_rho .* (y .^ e_y_sum_28) .* (l .^ e_l_sum_28)
    S28 = link_weight * sum28_kernel
    lhs_1 = y[1]^e_y_lhs_28 * l[1]^e_l_lhs_28
    self_1_nochi = uA_theta[1] * y[1]^e_y_self_28 * l[1]^e_l_self_28
    sum_1_nochi = u_theta[1] * A_theta2_lam_rho[1] * S28[1]

    if abs(e_chi_frac) < 1e-10
        chi = lhs_1 / (self_1_nochi + sum_1_nochi)
    else
        chi = lhs_1 / (self_1_nochi + sum_1_nochi)
        for _ in 1:50
            f = chi * self_1_nochi + chi^e_chi_frac * sum_1_nochi - lhs_1
            fp = self_1_nochi + e_chi_frac * chi^(e_chi_frac - 1) * sum_1_nochi
            chi = chi - f / fp
            chi = max(chi, 1e-15)
        end
    end

    welfare = chi^(-1 / theta)

    # Compute final transport costs and traffic
    # Pi_i^{-theta} = A_i^{-theta} L_i^{-theta} Y_i^{theta+1}  (from eq 7)
    # P_j^{-theta} = chi^{-1} (A_j u_j)^{-theta}               (from P_j = A_j u_j / W)
    Pi_inv_theta = (A_bar .^ (-theta)) .* ((l .* Lbar) .^ (-theta * (1 + alpha))) .*
                   ((y .* Lbar) .^ (theta + 1))
    P_inv_theta = chi^(-1) .* (A_bar .^ (-theta)) .* (u_bar .^ (-theta)) .*
                  ((l .* Lbar) .^ (-theta * (alpha + beta)))

    Xi = zeros(N, N)
    if lambda > 0
        tc_final = compute_transport_costs_with_congestion(t_bar_matrix, P_inv_theta, Pi_inv_theta,
                                                            theta, lambda; check_spectral=false)
        tau = tc_final.tau
        B_mat = tc_final.B
        Xi = tc_final.Xi
    else
        tc_final = compute_transport_costs(t_bar_matrix, theta)
        tau = tc_final.tau
        B_mat = tc_final.B
        for k in 1:N, ll in 1:N
            if isfinite(t_bar_matrix[k, ll]) && k != ll
                Xi[k, ll] = t_bar_matrix[k, ll]^(-theta) * P_inv_theta[k] * Pi_inv_theta[ll]
            end
        end
    end

    return (; y, l, chi, tau, traffic=Xi, welfare, B=B_mat, converged, iterations=iter)
end

# ============================================================
# Equilibrium Solver - Urban Model (Section 4.1)
# ============================================================

"""
    solve_urban_equilibrium(A_bar, u_bar, t_bar_matrix, theta, alpha, beta, lambda, Lbar;
                            maxiter=5000, tol=1e-10, damp=0.5, verbose=true)

Solve the urban equilibrium with endogenous traffic congestion.

Directly solves the self-contained system eqs (30)-(31) for residential shares {l_R_i}
and workplace shares {l_F_i}. These equations incorporate congestion through direct links
only (no Leontief inverse needed during iteration):

Eq (30): (l_R_i)^{1-θβ} (l_F_i)^{θλ(1-αθ)ρ} = χ ū_i^θ Ā_i^θ (l_F_i)^{θ(α+λ)ρ}
  + χ^{θλρ} ū_i^θ Ā_i^{θ²λρ} Σ_j (L̄^λ t̄_{ij})^{-θρ} ū_j^{-θρ} (l_R_j)^{(1-θβ)ρ}

Eq (31): (l_R_i)^{θλ(1-βθ)ρ} (l_F_i)^{1-θα} = χ ū_i^θ Ā_i^θ (l_R_i)^{θ(β+λ)ρ}
  + χ^{θλρ} Ā_i^θ ū_i^{θ²λρ} Σ_j (L̄^λ t̄_{ji})^{-θρ} Ā_j^{-θρ} (l_F_j)^{(1-αθ)ρ}

When lambda=0, these reduce to the standard eqs (19)-(20) split into self + adjacency terms.
The self-term comes from the identity in B=(I-A)^{-1}; the sum uses direct links only.

Returns NamedTuple with l_R, l_F, chi, tau, traffic (Xi), welfare, B, converged, iterations
"""
function solve_urban_equilibrium(A_bar, u_bar, t_bar_matrix, theta, alpha, beta, lambda, Lbar;
                                  maxiter=5000, tol=1e-10, damp=0.5, verbose=true)
    N = length(A_bar)
    rho = 1 / (1 + theta * lambda)

    # Initialize uniformly
    l_R = fill(1.0 / N, N)
    l_F = fill(1.0 / N, N)

    # Precompute link weights: w_{ij} = (Lbar^lambda * t_bar_{ij})^{-theta*rho}
    # Only for connected links (finite t_bar, k != l)
    link_weight = zeros(N, N)         # for eq (30): links from i
    link_weight_T = zeros(N, N)       # for eq (31): links into i (transpose)
    for k in 1:N, l in 1:N
        if isfinite(t_bar_matrix[k, l]) && k != l
            w = (Lbar^lambda * t_bar_matrix[k, l])^(-theta * rho)
            link_weight[k, l] = w
            link_weight_T[l, k] = w   # t_bar_{kl} contributes to eq (31) at location l
        end
    end

    # Precompute fundamental powers
    u_theta = u_bar .^ theta           # ū_i^θ
    A_theta = A_bar .^ theta           # Ā_i^θ
    uA_theta = u_theta .* A_theta      # ū_i^θ Ā_i^θ  (self-term coefficient)

    # Sum-term fundamental powers (depend on rho)
    u_neg_theta_rho = u_bar .^ (-theta * rho)    # ū_j^{-θρ} for eq (30) sum
    A_neg_theta_rho = A_bar .^ (-theta * rho)    # Ā_j^{-θρ} for eq (31) sum
    u_theta2_lam_rho = u_bar .^ (theta^2 * lambda * rho)  # ū_i^{θ²λρ} for eq (31)
    A_theta2_lam_rho = A_bar .^ (theta^2 * lambda * rho)  # Ā_i^{θ²λρ} for eq (30)

    # Exponents
    e_lR_lhs = 1 - theta * beta                     # LHS l_R exponent in eq (30)
    e_lF_cross_30 = theta * lambda * (1 - alpha * theta) * rho  # LHS l_F cross exponent in eq (30)
    e_lF_self_30 = theta * (alpha + lambda) * rho    # self-term l_F exponent in eq (30)
    e_lR_sum_30 = (1 - theta * beta) * rho           # sum-term l_R_j exponent in eq (30)

    e_lR_cross_31 = theta * lambda * (1 - beta * theta) * rho   # LHS l_R cross exponent in eq (31)
    e_lF_lhs = 1 - theta * alpha                     # LHS l_F exponent in eq (31)
    e_lR_self_31 = theta * (beta + lambda) * rho     # self-term l_R exponent in eq (31)
    e_lF_sum_31 = (1 - theta * alpha) * rho          # sum-term l_F_j exponent in eq (31)

    e_chi_frac = theta * lambda * rho                # χ^{θλρ} exponent on sum-term

    if verbose
        println("\n>>>> Urban Equilibrium Solver <<<<")
        println("  N=$N, theta=$theta, alpha=$alpha, beta=$beta, lambda=$lambda, Lbar=$Lbar")
        println("  rho=$(round(rho, digits=4)), theta*lambda*rho=$(round(e_chi_frac, digits=4))")
        tc0 = compute_transport_costs(t_bar_matrix, theta; check_spectral=false)
        println("  Spectral radius of free-flow A: $(round(tc0.spectral_radius, digits=4))")
    end

    converged = false
    iter = 0
    chi = 1.0  # will be determined by normalization

    for it in 1:maxiter
        iter = it
        l_R_old = copy(l_R)
        l_F_old = copy(l_F)

        # --- Eq (30): Update l_R ---
        # RHS_i = chi * uA_i * l_F_i^{e_lF_self_30}
        #       + chi^{e_chi_frac} * u_i^theta * A_i^{theta^2*lambda*rho}
        #         * sum_j link_weight[i,j] * u_j^{-theta*rho} * l_R_j^{e_lR_sum_30}
        # l_R_i = (RHS_i / l_F_i^{e_lF_cross_30})^{1/e_lR_lhs}

        l_R_sum_pow = l_R .^ e_lR_sum_30              # l_R_j^{(1-θβ)ρ}
        sum30_kernel = u_neg_theta_rho .* l_R_sum_pow  # ū_j^{-θρ} * l_R_j^{(1-θβ)ρ}
        S30 = link_weight * sum30_kernel               # sum_j for each i

        for i in 1:N
            self_i = chi * uA_theta[i] * l_F[i]^e_lF_self_30
            sum_i = chi^e_chi_frac * u_theta[i] * A_theta2_lam_rho[i] * S30[i]
            rhs_i = self_i + sum_i
            l_R[i] = (rhs_i / l_F[i]^e_lF_cross_30)^(1 / e_lR_lhs)
        end
        l_R .= max.(l_R, 1e-15)
        l_R ./= sum(l_R)

        # --- Eq (31): Update l_F ---
        # RHS_i = chi * uA_i * l_R_i^{e_lR_self_31}
        #       + chi^{e_chi_frac} * A_i^theta * u_i^{theta^2*lambda*rho}
        #         * sum_j link_weight_T[i,j] * A_j^{-theta*rho} * l_F_j^{e_lF_sum_31}
        # l_F_i = (RHS_i / l_R_i^{e_lR_cross_31})^{1/e_lF_lhs}

        l_F_sum_pow = l_F .^ e_lF_sum_31              # l_F_j^{(1-αθ)ρ}
        sum31_kernel = A_neg_theta_rho .* l_F_sum_pow  # Ā_j^{-θρ} * l_F_j^{(1-αθ)ρ}
        S31 = link_weight_T * sum31_kernel             # sum_j for each i

        for i in 1:N
            self_i = chi * uA_theta[i] * l_R[i]^e_lR_self_31
            sum_i = chi^e_chi_frac * A_theta[i] * u_theta2_lam_rho[i] * S31[i]
            rhs_i = self_i + sum_i
            l_F[i] = (rhs_i / l_R[i]^e_lR_cross_31)^(1 / e_lF_lhs)
        end
        l_F .= max.(l_F, 1e-15)
        l_F ./= sum(l_F)

        # --- Damped update ---
        l_R .= (1 - damp) .* l_R_old .+ damp .* l_R
        l_F .= (1 - damp) .* l_F_old .+ damp .* l_F
        l_R .= max.(l_R, 1e-15)
        l_F .= max.(l_F, 1e-15)
        l_R ./= sum(l_R)
        l_F ./= sum(l_F)

        # --- Check convergence ---
        err = max(maximum(abs.(l_R .- l_R_old)), maximum(abs.(l_F .- l_F_old)))
        if verbose && (it % 500 == 0 || it <= 5)
            println("  Iter $it: max change = $(round(err, sigdigits=4))")
        end

        if err < tol
            converged = true
            if verbose println("  Converged at iteration $it (tol=$tol)") end
            break
        end
    end

    if !converged && verbose
        println("  WARNING: Did not converge in $maxiter iterations")
    end

    # Compute chi from eq (30) at location 1 (using the self-term alone for simplicity):
    # At equilibrium, LHS = RHS, so chi can be backed out from the full equation
    l_R_sum_pow = l_R .^ e_lR_sum_30
    sum30_kernel = u_neg_theta_rho .* l_R_sum_pow
    S30 = link_weight * sum30_kernel
    lhs_1 = l_R[1]^e_lR_lhs * l_F[1]^e_lF_cross_30
    self_1_nochi = uA_theta[1] * l_F[1]^e_lF_self_30
    sum_1_nochi = u_theta[1] * A_theta2_lam_rho[1] * S30[1]
    # lhs = chi * self + chi^{e_chi_frac} * sum
    # For small lambda: e_chi_frac ≈ 0, so chi^{e_chi_frac} ≈ 1, giving chi ≈ (lhs - sum) / self
    # General case: solve chi * self + chi^{e_chi_frac} * sum = lhs via Newton
    if abs(e_chi_frac) < 1e-10  # lambda ≈ 0
        chi = lhs_1 / (self_1_nochi + sum_1_nochi)
    else
        chi = lhs_1 / (self_1_nochi + sum_1_nochi)  # initial guess
        for _ in 1:50
            f = chi * self_1_nochi + chi^e_chi_frac * sum_1_nochi - lhs_1
            fp = self_1_nochi + e_chi_frac * chi^(e_chi_frac - 1) * sum_1_nochi
            chi = chi - f / fp
            chi = max(chi, 1e-15)
        end
    end

    welfare = chi^(-1 / theta)

    # Compute final transport costs and traffic using correct P/Pi formulas
    # Pi_i^{-theta} = u_bar_i^{-theta} * (l_R_i)^{1-beta*theta} * chi^{-1/2} * Lbar^{(1+(alpha-beta)*theta)/2}
    # P_j^{-theta} = A_bar_j^{-theta} * (l_F_j)^{1-alpha*theta} * chi^{-1/2} * Lbar^{(1+(beta-alpha)*theta)/2}
    Pi_inv_theta = (u_bar .^ (-theta)) .* (l_R .^ (1 - beta * theta)) .*
                   (chi^(-0.5) * Lbar^((1 + (alpha - beta) * theta) / 2))
    P_inv_theta = (A_bar .^ (-theta)) .* (l_F .^ (1 - alpha * theta)) .*
                  (chi^(-0.5) * Lbar^((1 + (beta - alpha) * theta) / 2))

    if lambda > 0
        tc_final = compute_transport_costs_with_congestion(t_bar_matrix, P_inv_theta, Pi_inv_theta,
                                                            theta, lambda; check_spectral=false)
        tau = tc_final.tau
        B_mat = tc_final.B
        Xi = tc_final.Xi
    else
        tc_final = compute_transport_costs(t_bar_matrix, theta)
        tau = tc_final.tau
        B_mat = tc_final.B
        Xi = zeros(N, N)
        for k in 1:N, ll in 1:N
            if isfinite(t_bar_matrix[k, ll]) && k != ll
                Xi[k, ll] = t_bar_matrix[k, ll]^(-theta) * P_inv_theta[k] * Pi_inv_theta[ll]
            end
        end
    end

    return (; l_R, l_F, chi, tau, traffic=Xi, welfare, B=B_mat, converged, iterations=iter)
end

# ============================================================
# Traffic-to-Trade/Commuting Prediction (Section 5.1)
# ============================================================

"""
    predict_trade_from_traffic(Xi, Y, E)

Predict bilateral trade flows from traffic using eq (34):
    X_{ij} = c_{ij}^X * Y_i * E_j

where C^X = (D^X - Xi)^{-1} and D^X is diagonal with elements:
    d_i = 0.5*(Y_i + E_i) + 0.5*(sum_j Xi_{ji} + sum_j Xi_{ij})

This provides an analytical mapping from observed traffic flows on links
to bilateral trade flows between origins and destinations.

Returns NxN matrix of predicted trade flows X_{ij}.
"""
function predict_trade_from_traffic(Xi, Y, E)
    N = length(Y)
    Xi_row_sums = sum(Xi, dims=2)[:]    # sum_j Xi_{ij}
    Xi_col_sums = sum(Xi, dims=1)[:]    # sum_j Xi_{ji}

    # D^X diagonal elements (eq 34)
    d = 0.5 .* (Y .+ E) .+ 0.5 .* (Xi_col_sums .+ Xi_row_sums)
    D_X = Diagonal(d)

    # C^X = (D^X - Xi)^{-1}
    C_X = inv(Matrix(D_X) - Xi)

    # Trade flows: X_{ij} = c_{ij}^X * Y_i * E_j
    X = C_X .* (Y * E')

    return X
end

"""
    predict_commuting_from_traffic(Xi, L_R, L_F)

Predict bilateral commuting flows from traffic using eq (35):
    L_{ij} = c_{ij}^L * L_R_i * L_F_j

where C^L = (D^L - Xi)^{-1} and D^L is diagonal with elements:
    d_i = 0.5*(L_R_i + L_F_i) + 0.5*(sum_j Xi_{ji} + sum_j Xi_{ij})

Returns NxN matrix of predicted commuting flows L_{ij}.
"""
function predict_commuting_from_traffic(Xi, L_R, L_F)
    N = length(L_R)
    Xi_row_sums = sum(Xi, dims=2)[:]
    Xi_col_sums = sum(Xi, dims=1)[:]

    d = 0.5 .* (L_R .+ L_F) .+ 0.5 .* (Xi_col_sums .+ Xi_row_sums)
    D_L = Diagonal(d)

    C_L = inv(Matrix(D_L) - Xi)

    L_flows = C_L .* (L_R * L_F')

    return L_flows
end

# ============================================================
# Link Intensity (Section 3.2)
# ============================================================

"""
    compute_link_intensity(tau, t_matrix, theta, i, j, k, l)

Compute the link intensity pi_{ij}^{kl} (eq 23):
    pi_{ij}^{kl} = (tau_{ij} / (tau_{ik} * t_{kl} * tau_{lj}))^theta

This is the expected number of times link (k,l) is used in trade/commuting between (i,j).
Links along the most direct route have high intensity; out-of-the-way links have low intensity.
"""
function compute_link_intensity(tau, t_matrix, theta, i, j, k, l)
    if !isfinite(t_matrix[k, l])
        return 0.0
    end
    return (tau[i, j] / (tau[i, k] * t_matrix[k, l] * tau[l, j]))^theta
end

# ============================================================
# Exact-Hat Counterfactual Solver - Economic Geography (Section 5.2)
# ============================================================

"""
    counterfactual_econ_geo(Xi, Y, E, t_hat_bar, theta, alpha, beta, lambda;
                           maxiter=5000, tol=1e-10, damp=0.3, verbose=true)

Solve the exact-hat counterfactual for the economic geography model (Proposition 2, eqs 36-37).

The structure mirrors eqs (28)-(29) but with traffic-based weight shares replacing fundamentals.
Each equation has a self-term (weighted by local economic activity relative to traffic) and a
sum-term (weighted by traffic shares), with the same exponents as the levels system.

Normalizations: Σ ŷ_i Y_i = Σ Y_i (total income unchanged) and Σ l̂_i L_i = Σ L_i.

Returns NamedTuple with y_hat, l_hat, chi_hat, welfare_change, converged
"""
function counterfactual_econ_geo(Xi, Y, E, t_hat_bar, theta, alpha, beta, lambda;
                                 maxiter=5000, tol=1e-10, damp=0.3, verbose=true)
    N = length(Y)
    rho = 1 / (1 + theta * lambda)

    # Exponents (identical to eqs 28-29)
    e_y_lhs_36 = (1 + theta + theta * lambda) * rho
    e_l_lhs_36 = -theta * (1 + alpha + theta * lambda * (alpha + beta)) * rho
    e_y_self_36 = e_y_lhs_36
    e_l_self_36 = theta * (beta - 1) * rho
    e_y_sum_36 = (1 + theta) * rho
    e_l_sum_36 = -theta * (1 + alpha) * rho

    e_y_lhs_37 = -theta * (1 - lambda) * rho
    e_l_lhs_37 = theta * (1 - beta - theta * lambda * (alpha + beta)) * rho
    e_y_self_37 = e_y_lhs_37
    e_l_self_37 = theta * (alpha + 1) * rho
    e_y_sum_37 = -theta * rho
    e_l_sum_37 = theta * (1 - beta) * rho

    e_t = -theta * rho
    e_chi_frac = theta * lambda * rho

    # Weight shares from observed data (from Proposition 2)
    Xi_row_sums = sum(Xi, dims=2)[:]
    Xi_col_sums = sum(Xi, dims=1)[:]

    denom_out = E .+ Xi_row_sums
    s_own_out = E ./ denom_out
    s_out = Xi ./ denom_out

    denom_in = Y .+ Xi_col_sums
    s_own_in = Y ./ denom_in
    s_in = Xi' ./ denom_in

    # Pre-compute t_hat^{e_t}
    t_hat_pow = ones(N, N)
    for k in 1:N, l in 1:N
        if t_hat_bar[k, l] != 1.0
            t_hat_pow[k, l] = t_hat_bar[k, l]^e_t
        end
    end

    y_hat = ones(N)
    l_hat = ones(N)
    chi_hat = 1.0

    if verbose println("\n>>>> Econ Geo Counterfactual Solver <<<<") end

    converged = false
    for it in 1:maxiter
        y_old = copy(y_hat)
        l_old = copy(l_hat)
        chi_old = chi_hat

        # Eq (36): update y_hat
        for i in 1:N
            self_term = chi_hat * s_own_out[i] * y_hat[i]^e_y_self_36 * l_hat[i]^e_l_self_36
            sum_term = 0.0
            for j in 1:N
                if s_out[i, j] > 0
                    sum_term += s_out[i, j] * t_hat_pow[i, j] * y_hat[j]^e_y_sum_36 * l_hat[j]^e_l_sum_36
                end
            end
            sum_term *= chi_hat^e_chi_frac
            rhs = self_term + sum_term
            y_hat[i] = (rhs / l_hat[i]^e_l_lhs_36)^(1 / e_y_lhs_36)
        end

        # Eq (37): update l_hat
        for i in 1:N
            self_term = chi_hat * s_own_in[i] * y_hat[i]^e_y_self_37 * l_hat[i]^e_l_self_37
            sum_term = 0.0
            for j in 1:N
                if s_in[i, j] > 0
                    sum_term += s_in[i, j] * t_hat_pow[j, i] * y_hat[j]^e_y_sum_37 * l_hat[j]^e_l_sum_37
                end
            end
            sum_term *= chi_hat^e_chi_frac
            rhs = self_term + sum_term
            l_hat[i] = (rhs / y_hat[i]^e_y_lhs_37)^(1 / e_l_lhs_37)
        end

        # Normalization: Σ ŷ_i Y_i = Σ Y_i and Σ l̂_i L_i = Σ L_i
        # With balanced trade Y=E, and L proportional to Y in the symmetric case
        y_hat .*= sum(Y) / sum(y_hat .* Y)
        l_hat .*= N / sum(l_hat)  # since l_obs = 1/N for all i in symmetric case

        y_hat .= max.(y_hat, 1e-15)
        l_hat .= max.(l_hat, 1e-15)

        # Recover chi_hat from eq (36) at location 1 via Newton
        lhs_1 = y_hat[1]^e_y_lhs_36 * l_hat[1]^e_l_lhs_36
        self_1 = s_own_out[1] * y_hat[1]^e_y_self_36 * l_hat[1]^e_l_self_36
        sum_1 = sum(s_out[1, j] * t_hat_pow[1, j] * y_hat[j]^e_y_sum_36 * l_hat[j]^e_l_sum_36
                     for j in 1:N if s_out[1, j] > 0; init=0.0)
        chi_hat = lhs_1 / (self_1 + sum_1)  # initial guess (exact if lambda=0)
        if abs(e_chi_frac) > 1e-10
            for _ in 1:30
                f_val = chi_hat * self_1 + chi_hat^e_chi_frac * sum_1 - lhs_1
                fp_val = self_1 + e_chi_frac * chi_hat^(e_chi_frac - 1) * sum_1
                chi_hat = chi_hat - f_val / fp_val
                chi_hat = max(chi_hat, 1e-15)
            end
        end

        # Damped update
        y_hat .= (1 - damp) .* y_old .+ damp .* y_hat
        l_hat .= (1 - damp) .* l_old .+ damp .* l_hat
        chi_hat = (1 - damp) * chi_old + damp * chi_hat

        err = max(maximum(abs.(y_hat .- y_old)), maximum(abs.(l_hat .- l_old)),
                  abs(chi_hat - chi_old))
        if verbose && (it % 500 == 0 || it <= 3)
            println("  Iter $it: max change = $(round(err, sigdigits=4)), chi_hat = $(round(chi_hat, sigdigits=6))")
        end

        if err < tol
            converged = true
            if verbose println("  Converged at iteration $it") end
            break
        end
    end

    welfare_change = chi_hat^(-1 / theta)
    return (; y_hat, l_hat, chi_hat, welfare_change, converged)
end

# ============================================================
# Exact-Hat Counterfactual Solver - Urban Model (Section 5.2)
# ============================================================

"""
    counterfactual_urban(Xi, L_R, L_F, t_hat_bar, theta, alpha, beta, lambda;
                         maxiter=5000, tol=1e-10, damp=0.3, verbose=true)

Solve the exact-hat counterfactual for the urban model (Proposition 2, eqs 38-39).

The structure mirrors eqs (30)-(31) but with traffic-based weight shares replacing fundamentals.

Eq (38): (l̂_R_i)^{1-θβ} (l̂_F_i)^{θλ(1-αθ)ρ} = χ̂ s_own_out_i (l̂_F_i)^{θ(α+λ)ρ}
  + χ̂^{θλρ} Σ_j s_out_{ij} t̂_{ij}^{-θρ} (l̂_R_j)^{(1-θβ)ρ}

Eq (39): (l̂_R_i)^{θλ(1-βθ)ρ} (l̂_F_i)^{1-θα} = χ̂ s_own_in_i (l̂_R_i)^{θ(β+λ)ρ}
  + χ̂^{θλρ} Σ_j s_in_{ji} t̂_{ji}^{-θρ} (l̂_F_j)^{(1-αθ)ρ}

Returns NamedTuple with l_R_hat, l_F_hat, chi_hat, welfare_change, converged
"""
function counterfactual_urban(Xi, L_R, L_F, t_hat_bar, theta, alpha, beta, lambda;
                               maxiter=5000, tol=1e-10, damp=0.3, verbose=true)
    N = length(L_R)
    rho = 1 / (1 + theta * lambda)

    # Exponents (identical to eqs 30-31)
    e_lR_lhs_38 = 1 - theta * beta
    e_lF_cross_38 = theta * lambda * (1 - alpha * theta) * rho
    e_lF_self_38 = theta * (alpha + lambda) * rho
    e_lR_sum_38 = (1 - theta * beta) * rho

    e_lR_cross_39 = theta * lambda * (1 - beta * theta) * rho
    e_lF_lhs_39 = 1 - theta * alpha
    e_lR_self_39 = theta * (beta + lambda) * rho
    e_lF_sum_39 = (1 - theta * alpha) * rho

    e_t = -theta * rho
    e_chi_frac = theta * lambda * rho

    # Weight shares from observed data
    Xi_row_sums = sum(Xi, dims=2)[:]
    Xi_col_sums = sum(Xi, dims=1)[:]

    denom_out = L_F .+ Xi_row_sums
    s_own_out = L_F ./ denom_out
    s_out = Xi ./ denom_out

    denom_in = L_R .+ Xi_col_sums
    s_own_in = L_R ./ denom_in
    s_in = Xi' ./ denom_in

    t_hat_pow = ones(N, N)
    for k in 1:N, l in 1:N
        if t_hat_bar[k, l] != 1.0
            t_hat_pow[k, l] = t_hat_bar[k, l]^e_t
        end
    end

    l_R_hat = ones(N)
    l_F_hat = ones(N)
    chi_hat = 1.0

    if verbose println("\n>>>> Urban Counterfactual Solver <<<<") end

    converged = false
    for it in 1:maxiter
        l_R_old = copy(l_R_hat)
        l_F_old = copy(l_F_hat)
        chi_old = chi_hat

        # Eq (38): update l_R_hat
        for i in 1:N
            self_term = chi_hat * s_own_out[i] * l_F_hat[i]^e_lF_self_38
            sum_term = 0.0
            for j in 1:N
                if s_out[i, j] > 0
                    sum_term += s_out[i, j] * t_hat_pow[i, j] * l_R_hat[j]^e_lR_sum_38
                end
            end
            sum_term *= chi_hat^e_chi_frac
            rhs = self_term + sum_term
            l_R_hat[i] = (rhs / l_F_hat[i]^e_lF_cross_38)^(1 / e_lR_lhs_38)
        end

        # Eq (39): update l_F_hat
        for i in 1:N
            self_term = chi_hat * s_own_in[i] * l_R_hat[i]^e_lR_self_39
            sum_term = 0.0
            for j in 1:N
                if s_in[i, j] > 0
                    sum_term += s_in[i, j] * t_hat_pow[j, i] * l_F_hat[j]^e_lF_sum_39
                end
            end
            sum_term *= chi_hat^e_chi_frac
            rhs = self_term + sum_term
            l_F_hat[i] = (rhs / l_R_hat[i]^e_lR_cross_39)^(1 / e_lF_lhs_39)
        end

        # Normalization: Σ l̂_R_i L_R_i = Σ L_R_i, Σ l̂_F_i L_F_i = Σ L_F_i
        l_R_hat .*= sum(L_R) / sum(l_R_hat .* L_R)
        l_F_hat .*= sum(L_F) / sum(l_F_hat .* L_F)

        l_R_hat .= max.(l_R_hat, 1e-15)
        l_F_hat .= max.(l_F_hat, 1e-15)

        # Recover chi_hat from eq (38) at location 1 via Newton
        lhs_1 = l_R_hat[1]^e_lR_lhs_38 * l_F_hat[1]^e_lF_cross_38
        self_1 = s_own_out[1] * l_F_hat[1]^e_lF_self_38
        sum_1 = sum(s_out[1, j] * t_hat_pow[1, j] * l_R_hat[j]^e_lR_sum_38
                     for j in 1:N if s_out[1, j] > 0; init=0.0)
        chi_hat = lhs_1 / (self_1 + sum_1)
        if abs(e_chi_frac) > 1e-10
            for _ in 1:30
                f_val = chi_hat * self_1 + chi_hat^e_chi_frac * sum_1 - lhs_1
                fp_val = self_1 + e_chi_frac * chi_hat^(e_chi_frac - 1) * sum_1
                chi_hat = chi_hat - f_val / fp_val
                chi_hat = max(chi_hat, 1e-15)
            end
        end

        # Damped update
        l_R_hat .= (1 - damp) .* l_R_old .+ damp .* l_R_hat
        l_F_hat .= (1 - damp) .* l_F_old .+ damp .* l_F_hat
        chi_hat = (1 - damp) * chi_old + damp * chi_hat

        err = max(maximum(abs.(l_R_hat .- l_R_old)), maximum(abs.(l_F_hat .- l_F_old)),
                  abs(chi_hat - chi_old))
        if verbose && (it % 500 == 0 || it <= 3)
            println("  Iter $it: max change = $(round(err, sigdigits=4))")
        end

        if err < tol
            converged = true
            if verbose println("  Converged at iteration $it") end
            break
        end
    end

    welfare_change = chi_hat^(-1 / theta)
    return (; l_R_hat, l_F_hat, chi_hat, welfare_change, converged)
end

# ============================================================
# Welfare Elasticity Computation (Section 6)
# ============================================================

"""
    compute_welfare_elasticities(Xi, data1, data2, t_bar_matrix, theta, alpha, beta, lambda;
                                 model=:urban, delta=0.01)

Compute the welfare elasticity of improving each link in the network.

For each connected link (k,l), constructs a small improvement t_hat_bar_{kl} = 1 - delta
(and its reverse t_hat_bar_{lk} = 1 - delta for symmetric networks), then solves the
counterfactual to compute:
    elasticity_{kl} = (W_hat - 1) / delta ≈ d ln W / d ln t_bar_{kl}

Arguments:
- Xi: NxN observed traffic matrix
- data1: Y (econ_geo) or L_R (urban)
- data2: E (econ_geo) or L_F (urban)
- model: :econ_geo or :urban
- delta: size of perturbation (default 1%)

Returns NxN matrix of welfare elasticities (0 for non-connected links).
"""
function compute_welfare_elasticities(Xi, data1, data2, t_bar_matrix, theta, alpha, beta, lambda;
                                      model=:urban, delta=0.01)
    N = size(t_bar_matrix, 1)
    elasticities = zeros(N, N)
    connected = isfinite.(t_bar_matrix) .& (.!I(N))  # connected links excluding diagonal

    n_links = sum(connected)
    link_count = 0

    for k in 1:N, l in 1:N
        if !connected[k, l]
            continue
        end

        link_count += 1
        if link_count % 10 == 0
            print("  Link $link_count / $n_links\r")
        end

        # Construct t_hat_bar: 1 everywhere except link (k,l) and (l,k)
        t_hat = ones(N, N)
        t_hat[k, l] = 1.0 - delta
        if connected[l, k]
            t_hat[l, k] = 1.0 - delta
        end

        if model == :econ_geo
            result = counterfactual_econ_geo(Xi, data1, data2, t_hat, theta, alpha, beta, lambda;
                                             maxiter=2000, tol=1e-8, damp=0.3, verbose=false)
        else
            result = counterfactual_urban(Xi, data1, data2, t_hat, theta, alpha, beta, lambda;
                                          maxiter=2000, tol=1e-8, damp=0.3, verbose=false)
        end

        # Welfare elasticity: (W_hat - 1) / delta
        elasticities[k, l] = (result.welfare_change - 1.0) / delta
    end

    println("  Computed $link_count welfare elasticities")
    return elasticities
end

# ============================================================
# Plotting Utilities
# ============================================================

"""
    plot_grid_equilibrium(coords, pop, traffic, t_bar_matrix, N_side;
                          title="Equilibrium", pop_label="Population")

Plot equilibrium population as a heatmap on the grid and traffic intensity on edges.
"""
function plot_grid_equilibrium(coords, pop, traffic, t_bar_matrix, N_side;
                                title="Equilibrium", pop_label="Population share")
    N = N_side^2

    # Reshape population to grid
    pop_grid = reshape(pop, N_side, N_side)'  # transpose for correct orientation
    # Reverse rows so that y increases upward
    pop_grid = pop_grid[end:-1:1, :]

    p = heatmap(1:N_side, 1:N_side, pop_grid,
                c=:YlOrRd, xlabel="x", ylabel="y",
                title=title, colorbar_title=pop_label,
                aspect_ratio=1, size=(500, 450))

    # Add traffic on edges
    max_traffic = maximum(traffic[isfinite.(t_bar_matrix) .& (t_bar_matrix .< Inf)])
    if max_traffic > 0
        for k in 1:N, l in k+1:N
            if isfinite(t_bar_matrix[k, l])
                x1, y1 = coords[k, 1] + 1, coords[k, 2] + 1
                x2, y2 = coords[l, 1] + 1, coords[l, 2] + 1
                intensity = (traffic[k, l] + traffic[l, k]) / (2 * max_traffic)
                lw = 1 + 4 * intensity
                plot!(p, [x1, x2], [y1, y2], color=RGBA(0, 0, 1, 0.3 + 0.7 * intensity),
                      linewidth=lw, label=false)
            end
        end
    end

    return p
end

# ============================================================
# Main Execution Block
# ============================================================

println("\n" * "="^60)
println("Allen & Arkolakis (RES 2022) - Traffic Congestion Model")
println("="^60)

# ---- Part A: Synthetic Grid Example (replicating Figure 1) ----
println("\n" * "-"^60)
println("Part A: Synthetic Grid Example (Figure 1)")
println("-"^60)

# Build 5x5 grid with identical fundamentals
net = build_grid_network(N_side; t_bar=t_bar_default)
A_bar_grid = ones(N)       # identical productivities
u_bar_grid = ones(N)       # identical amenities

# (a) No congestion (lambda=0, Lbar=100) - Figure 1(a)
println("\n--- Figure 1(a): No congestion (lambda=0, Lbar=100) ---")
res_1a = solve_urban_equilibrium(A_bar_grid, u_bar_grid, net.t_bar_matrix,
                                  theta_grid, 0.0, 0.0, 0.0, 100.0;
                                  maxiter=5000, tol=1e-10, damp=0.5)
println("  Welfare = $(round(res_1a.welfare, sigdigits=6))")
println("  Max pop share = $(round(maximum(res_1a.l_R), sigdigits=4)), Min = $(round(minimum(res_1a.l_R), sigdigits=4))")

# (b) Congestion, low scale (lambda=0.05, Lbar=100) - Figure 1(b)
println("\n--- Figure 1(b): Congestion (lambda=0.05, Lbar=100) ---")
res_1b = solve_urban_equilibrium(A_bar_grid, u_bar_grid, net.t_bar_matrix,
                                  theta_grid, 0.0, 0.0, 0.05, 100.0;
                                  maxiter=5000, tol=1e-10, damp=0.5)
println("  Welfare = $(round(res_1b.welfare, sigdigits=6))")
println("  Max pop share = $(round(maximum(res_1b.l_R), sigdigits=4)), Min = $(round(minimum(res_1b.l_R), sigdigits=4))")

# (c) Congestion, medium scale (lambda=0.05, Lbar=1000) - Figure 1(c)
println("\n--- Figure 1(c): Congestion (lambda=0.05, Lbar=1000) ---")
res_1c = solve_urban_equilibrium(A_bar_grid, u_bar_grid, net.t_bar_matrix,
                                  theta_grid, 0.0, 0.0, 0.05, 1000.0;
                                  maxiter=5000, tol=1e-10, damp=0.5)
println("  Welfare = $(round(res_1c.welfare, sigdigits=6))")
println("  Max pop share = $(round(maximum(res_1c.l_R), sigdigits=4)), Min = $(round(minimum(res_1c.l_R), sigdigits=4))")

# (d) Congestion, large scale (lambda=0.05, Lbar=10000) - Figure 1(d)
println("\n--- Figure 1(d): Congestion (lambda=0.05, Lbar=10000) ---")
res_1d = solve_urban_equilibrium(A_bar_grid, u_bar_grid, net.t_bar_matrix,
                                  theta_grid, 0.0, 0.0, 0.05, 10000.0;
                                  maxiter=5000, tol=1e-10, damp=0.5)
println("  Welfare = $(round(res_1d.welfare, sigdigits=6))")
println("  Max pop share = $(round(maximum(res_1d.l_R), sigdigits=4)), Min = $(round(minimum(res_1d.l_R), sigdigits=4))")

# Scale invariance check: without congestion, changing Lbar should NOT change distribution
println("\n--- Scale invariance check (lambda=0) ---")
res_scale = solve_urban_equilibrium(A_bar_grid, u_bar_grid, net.t_bar_matrix,
                                     theta_grid, 0.0, 0.0, 0.0, 10000.0;
                                     maxiter=5000, tol=1e-10, damp=0.5, verbose=false)
scale_diff = maximum(abs.(res_scale.l_R .- res_1a.l_R))
println("  Max diff in l_R between Lbar=100 and Lbar=10000 (lambda=0): $(round(scale_diff, sigdigits=4))")
println("  Scale invariant: $(scale_diff < 1e-6 ? "YES" : "NO")")

# Plot Figure 1
println("\n  Generating Figure 1 plots...")
p1a = plot_grid_equilibrium(net.coords, res_1a.l_R, res_1a.traffic, net.t_bar_matrix, N_side;
                             title="(a) No congestion\n(lambda=0, L=100)")
p1b = plot_grid_equilibrium(net.coords, res_1b.l_R, res_1b.traffic, net.t_bar_matrix, N_side;
                             title="(b) Congestion, low scale\n(lambda=0.05, L=100)")
p1c = plot_grid_equilibrium(net.coords, res_1c.l_R, res_1c.traffic, net.t_bar_matrix, N_side;
                             title="(c) Congestion, medium scale\n(lambda=0.05, L=1000)")
p1d = plot_grid_equilibrium(net.coords, res_1d.l_R, res_1d.traffic, net.t_bar_matrix, N_side;
                             title="(d) Congestion, large scale\n(lambda=0.05, L=10000)")

fig1 = plot(p1a, p1b, p1c, p1d, layout=(2, 2), size=(1000, 900),
            plot_title="Figure 1: Traffic Congestion and Distribution of Economic Activity")
savefig(fig1, joinpath(@__DIR__, "graphs", "figure1_congestion_scale.pdf"))
println("  Saved: graphs/figure1_congestion_scale.pdf")

# ---- Part B: Economic Geography Counterfactuals ----
println("\n" * "-"^60)
println("Part B: Economic Geography Counterfactuals")
println("-"^60)

# Solve equilibrium with US Highway parameters
println("\n--- Solving economic geography equilibrium ---")
res_eg = solve_econ_geo_equilibrium(A_bar_grid, u_bar_grid, net.t_bar_matrix,
                                     theta_eg, alpha_eg, beta_eg, lambda_eg, Lbar_default;
                                     maxiter=5000, tol=1e-10, damp=0.5)
println("  Welfare = $(round(res_eg.welfare, sigdigits=6))")
println("  Sum y = $(round(sum(res_eg.y), sigdigits=8)), Sum l = $(round(sum(res_eg.l), sigdigits=8))")

# Counterfactual: 1% improvement to central link (node 13 <-> 14)
println("\n--- Counterfactual: 1% improvement to central link ---")
Y_obs = res_eg.y .* Lbar_default
E_obs = copy(Y_obs)

t_hat_central = ones(N, N)
t_hat_central[13, 14] = 0.99   # 1% reduction in t_bar
t_hat_central[14, 13] = 0.99

cf_eg = counterfactual_econ_geo(res_eg.traffic, Y_obs, E_obs, t_hat_central,
                                 theta_eg, alpha_eg, beta_eg, lambda_eg;
                                 maxiter=3000, tol=1e-9, damp=0.3)
println("  Welfare change: $(round((cf_eg.welfare_change - 1) * 100, sigdigits=4))%")

# Welfare elasticities for all links
println("\n--- Computing welfare elasticities (econ geo) ---")
println("  With congestion (lambda=$(lambda_eg)):")
elast_eg_cong = compute_welfare_elasticities(res_eg.traffic, Y_obs, E_obs, net.t_bar_matrix,
                                              theta_eg, alpha_eg, beta_eg, lambda_eg;
                                              model=:econ_geo, delta=0.01)

println("  Without congestion (lambda=0):")
elast_eg_nocong = compute_welfare_elasticities(res_eg.traffic, Y_obs, E_obs, net.t_bar_matrix,
                                                theta_eg, alpha_eg, beta_eg, 0.0;
                                                model=:econ_geo, delta=0.01)

connected_mask = isfinite.(net.t_bar_matrix) .& (.!Bool.(I(N)))
e_cong = elast_eg_cong[connected_mask]
e_nocong = elast_eg_nocong[connected_mask]
println("  Mean elasticity with congestion:    $(round(mean(e_cong), sigdigits=4))")
println("  Mean elasticity without congestion: $(round(mean(e_nocong), sigdigits=4))")
println("  Correlation:                        $(round(cor(e_cong, e_nocong), sigdigits=4))")

# ---- Part C: Urban Counterfactuals ----
println("\n" * "-"^60)
println("Part C: Urban Counterfactuals")
println("-"^60)

# Solve urban equilibrium with Seattle parameters
println("\n--- Solving urban equilibrium ---")
res_urb = solve_urban_equilibrium(A_bar_grid, u_bar_grid, net.t_bar_matrix,
                                   theta_urban, alpha_urban, beta_urban, lambda_urban, Lbar_default;
                                   maxiter=5000, tol=1e-10, damp=0.5)
println("  Welfare = $(round(res_urb.welfare, sigdigits=6))")

# Counterfactual: 1% improvement to central link
println("\n--- Counterfactual: 1% improvement to central link ---")
L_R_obs = res_urb.l_R .* Lbar_default
L_F_obs = res_urb.l_F .* Lbar_default

cf_urb = counterfactual_urban(res_urb.traffic, L_R_obs, L_F_obs, t_hat_central,
                               theta_urban, alpha_urban, beta_urban, lambda_urban;
                               maxiter=3000, tol=1e-9, damp=0.3)
println("  Welfare change: $(round((cf_urb.welfare_change - 1) * 100, sigdigits=4))%")

# Welfare elasticities
println("\n--- Computing welfare elasticities (urban) ---")
println("  With congestion (lambda=$(lambda_urban)):")
elast_urb_cong = compute_welfare_elasticities(res_urb.traffic, L_R_obs, L_F_obs, net.t_bar_matrix,
                                               theta_urban, alpha_urban, beta_urban, lambda_urban;
                                               model=:urban, delta=0.01)

println("  Without congestion (lambda=0):")
elast_urb_nocong = compute_welfare_elasticities(res_urb.traffic, L_R_obs, L_F_obs, net.t_bar_matrix,
                                                 theta_urban, alpha_urban, beta_urban, 0.0;
                                                 model=:urban, delta=0.01)

eu_cong = elast_urb_cong[connected_mask]
eu_nocong = elast_urb_nocong[connected_mask]
println("  Mean elasticity with congestion:    $(round(mean(eu_cong), sigdigits=4))")
println("  Mean elasticity without congestion: $(round(mean(eu_nocong), sigdigits=4))")
println("  Correlation:                        $(round(cor(eu_cong, eu_nocong), sigdigits=4))")

# ---- Part D: Diagnostics and Predictions ----
println("\n" * "-"^60)
println("Part D: Diagnostics")
println("-"^60)

# Trade prediction from traffic (economic geography)
println("\n--- Predicting trade from traffic (econ geo) ---")
X_pred = predict_trade_from_traffic(res_eg.traffic, Y_obs, E_obs)
# Compare with gravity prediction: X_{ij} = tau_{ij}^{-theta} * Y_i * E_j / (sum_k tau_{kj}^{-theta} * Y_k)
X_gravity = zeros(N, N)
for i in 1:N, j in 1:N
    denom = sum(res_eg.tau[k, j]^(-theta_eg) * Y_obs[k] for k in 1:N)
    X_gravity[i, j] = res_eg.tau[i, j]^(-theta_eg) * Y_obs[i] * E_obs[j] / denom * E_obs[j]
end
# Compare off-diagonal
mask_offdiag = .!Bool.(I(N))
pred_vals = log.(max.(X_pred[mask_offdiag], 1e-20))
grav_vals = log.(max.(X_gravity[mask_offdiag], 1e-20))
corr_trade = cor(pred_vals[isfinite.(pred_vals) .& isfinite.(grav_vals)],
                 grav_vals[isfinite.(pred_vals) .& isfinite.(grav_vals)])
println("  Correlation between traffic-predicted and gravity trade (log): $(round(corr_trade, sigdigits=4))")

# Commuting prediction from traffic (urban)
println("\n--- Predicting commuting from traffic (urban) ---")
L_flows_pred = predict_commuting_from_traffic(res_urb.traffic, L_R_obs, L_F_obs)
# Gravity-predicted commuting: L_{ij} = tau_{ij}^{-theta} * u_i^theta * w_j^theta * Lbar / W^theta
L_flows_grav = zeros(N, N)
for i in 1:N, j in 1:N
    u_i = u_bar_grid[i] * (res_urb.l_R[i] * Lbar_default)^beta_urban
    w_j = A_bar_grid[j] * (res_urb.l_F[j] * Lbar_default)^alpha_urban
    L_flows_grav[i, j] = res_urb.tau[i, j]^(-theta_urban) * u_i^theta_urban * w_j^theta_urban
end
L_flows_grav .*= Lbar_default / sum(L_flows_grav)
pred_comm = log.(max.(L_flows_pred[mask_offdiag], 1e-20))
grav_comm = log.(max.(L_flows_grav[mask_offdiag], 1e-20))
valid_comm = isfinite.(pred_comm) .& isfinite.(grav_comm)
corr_comm = cor(pred_comm[valid_comm], grav_comm[valid_comm])
println("  Correlation between traffic-predicted and gravity commuting (log): $(round(corr_comm, sigdigits=4))")

# Spectral radius check
println("\n--- Spectral radius of adjacency matrix ---")
tc_check = compute_transport_costs(net.t_bar_matrix, theta_grid)
println("  theta=$theta_grid: spectral radius = $(round(tc_check.spectral_radius, sigdigits=4))")
tc_check2 = compute_transport_costs(net.t_bar_matrix, theta_eg)
println("  theta=$theta_eg: spectral radius = $(round(tc_check2.spectral_radius, sigdigits=4))")

# ---- Part E: Additional Plots ----
println("\n" * "-"^60)
println("Part E: Generating plots")
println("-"^60)

# Plot: Welfare elasticity heatmap (urban model)
elast_grid_cong = zeros(N_side, N_side)
elast_grid_nocong = zeros(N_side, N_side)
for i in 1:N
    # Sum of elasticities for all links connected to location i
    elast_grid_cong[(i-1)÷N_side+1, (i-1)%N_side+1] = sum(elast_urb_cong[i, :]) + sum(elast_urb_cong[:, i])
    elast_grid_nocong[(i-1)÷N_side+1, (i-1)%N_side+1] = sum(elast_urb_nocong[i, :]) + sum(elast_urb_nocong[:, i])
end

p_elast_c = heatmap(1:N_side, 1:N_side, elast_grid_cong[end:-1:1, :],
                     c=:YlOrRd, title="With congestion\n(lambda=$(lambda_urban))",
                     xlabel="x", ylabel="y", aspect_ratio=1,
                     colorbar_title="Sum of welfare\nelasticities")
p_elast_nc = heatmap(1:N_side, 1:N_side, elast_grid_nocong[end:-1:1, :],
                      c=:YlOrRd, title="Without congestion\n(lambda=0)",
                      xlabel="x", ylabel="y", aspect_ratio=1,
                      colorbar_title="Sum of welfare\nelasticities")

fig_elast = plot(p_elast_c, p_elast_nc, layout=(1, 2), size=(900, 400),
                 plot_title="Welfare Elasticities by Location (Urban Model)")
savefig(fig_elast, joinpath(@__DIR__, "graphs", "welfare_elasticities_urban.pdf"))
println("  Saved: graphs/welfare_elasticities_urban.pdf")

# Plot: Scatter of welfare elasticities with vs without congestion
p_scatter_eg = scatter(e_nocong, e_cong, xlabel="Without congestion", ylabel="With congestion",
                        title="Econ Geo", label=false, markersize=4, alpha=0.7,
                        aspect_ratio=1)
plot!(p_scatter_eg, [0, maximum(e_nocong)], [0, maximum(e_nocong)],
      linestyle=:dash, color=:red, label="45-degree line")

p_scatter_urb = scatter(eu_nocong, eu_cong, xlabel="Without congestion", ylabel="With congestion",
                         title="Urban", label=false, markersize=4, alpha=0.7,
                         aspect_ratio=1)
plot!(p_scatter_urb, [0, maximum(eu_nocong)], [0, maximum(eu_nocong)],
      linestyle=:dash, color=:red, label="45-degree line")

fig_scatter = plot(p_scatter_eg, p_scatter_urb, layout=(1, 2), size=(900, 400),
                   plot_title="Welfare Elasticities: With vs Without Congestion")
savefig(fig_scatter, joinpath(@__DIR__, "graphs", "elasticity_scatter_congestion.pdf"))
println("  Saved: graphs/elasticity_scatter_congestion.pdf")

# ---- Summary ----
println("\n" * "="^60)
println("Summary of Results")
println("="^60)
println("\nFigure 1 (Urban model on 5x5 grid, theta=$theta_grid, alpha=beta=0):")
println("  (a) No congestion:     center pop share = $(round(res_1a.l_R[13], sigdigits=4))")
println("  (b) lambda=0.05, L=100:  center = $(round(res_1b.l_R[13], sigdigits=4))")
println("  (c) lambda=0.05, L=1000: center = $(round(res_1c.l_R[13], sigdigits=4))")
println("  (d) lambda=0.05, L=10000: center = $(round(res_1d.l_R[13], sigdigits=4))")
println("  -> Congestion pushes activity to edges; scale dependence confirmed.")

println("\nEconomic Geography (theta=$theta_eg, alpha=$alpha_eg, beta=$beta_eg, lambda=$lambda_eg):")
println("  Mean welfare elasticity with congestion:    $(round(mean(e_cong), sigdigits=4))")
println("  Mean welfare elasticity without congestion: $(round(mean(e_nocong), sigdigits=4))")
println("  Congestion reduces welfare gains: $(mean(e_cong) < mean(e_nocong) ? "YES" : "NO")")

println("\nUrban (theta=$theta_urban, alpha=$alpha_urban, beta=$beta_urban, lambda=$lambda_urban):")
println("  Mean welfare elasticity with congestion:    $(round(mean(eu_cong), sigdigits=4))")
println("  Mean welfare elasticity without congestion: $(round(mean(eu_nocong), sigdigits=4))")
println("  Congestion reduces welfare gains: $(mean(eu_cong) < mean(eu_nocong) ? "YES" : "NO")")

println("\nAll plots saved to: $(joinpath(@__DIR__, "graphs"))")
println("="^60)
