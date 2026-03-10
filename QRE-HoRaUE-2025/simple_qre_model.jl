# Simple Quantitative Regional Model
# Based on Allen and Arkolakis (2014), Section 4 of Allen & Arkolakis (2025)
# "Quantitative Regional Economics", Handbook of Regional and Urban Economics

# Implements: Armington trade, CES preferences, agglomeration/congestion forces,
# perfect labor mobility across locations.

using Random
using Statistics
using Plots
using LinearAlgebra
using StatsBase: geomean

# ****************************
# **** Parameters ****
# ****************************

N = 10                  # Number of locations
sigma = 4.0             # Elasticity of substitution between varieties (> 1)
alpha = 0.05            # Agglomeration elasticity on productivity (eq 22)
beta = -0.1             # Congestion elasticity on amenities (eq 23, negative = congestion)
Lbar = 1000.0           # Total labor supply

# Derived
sigma_tilde = (sigma - 1) / (2 * sigma - 1)  # Used in eq 47

# Uniqueness check (condition after eq 47)
rho_unique = (1 + sigma * alpha + beta * (sigma - 1)) / (1 - alpha * (sigma - 1) - beta * sigma)
println("Uniqueness parameter rho = $(round(rho_unique, digits=4)), |rho| = $(round(abs(rho_unique), digits=4))")
abs(rho_unique) <= 1 ? println("Unique equilibrium guaranteed.") : println("WARNING: |rho| > 1, uniqueness not guaranteed!")

# ****************************
# **** Geography ****
# ****************************

Random.seed!(1)

# Innate productivities (eq 22: A_i = Abar_i * L_i^alpha)
Abar = exp.(randn(N))
Abar ./= geomean(Abar)

# Innate amenities (eq 23: u_i = ubar_i * L_i^beta)
ubar = exp.(randn(N))
ubar ./= geomean(ubar)

@show round.(Abar, digits=3)
@show round.(ubar, digits=3)

# Symmetric trade cost matrix (iceberg: tau_ij >= 1, tau_ii = 1)
coords = range(0, stop=4, length=N)
dist_mat = [abs(coords[i] - coords[j]) for i in 1:N, j in 1:N]
tau = exp.(0.3 .* dist_mat)
for i in 1:N; tau[i,i] = 1.0; end

# Geography matrix T_ij = (Abar_i * ubar_j / tau_ij)^(sigma-1)  (used in eqs 39-40)
T_mat = [(Abar[i] * ubar[j] / tau[i,j])^(sigma - 1) for i in 1:N, j in 1:N]

# ============================================================
# CASE 1: Eigenvector Solution (alpha = beta = 0)
# Section 4.4, equations 39-40
# ============================================================

"""
    solve_eigenvector(T_mat, sigma, Lbar)

Solve the model without agglomeration/congestion (alpha=beta=0) using the
eigenvector characterization (eqs 39-40).

Returns NamedTuple with w, L, W (welfare), lambda, x (right eigvec), y (left eigvec).
"""
function solve_eigenvector(T_mat, sigma, Lbar)
    N = size(T_mat, 1)

    # Eq 39: lambda * x = T * x  (right eigenvector)
    # Eq 40: lambda * y = T' * y (left eigenvector)
    eig_r = eigen(T_mat)
    eig_l = eigen(T_mat')

    # Largest real eigenvalue (Perron-Frobenius: positive matrix => real, positive)
    idx_r = argmax(real.(eig_r.values))
    idx_l = argmax(real.(eig_l.values))

    lambda = real(eig_r.values[idx_r])
    x = real.(eig_r.vectors[:, idx_r])
    y = real.(eig_l.vectors[:, idx_l])

    # Ensure positive (flip sign if needed)
    if minimum(x) < 0; x = -x; end
    if minimum(y) < 0; y = -y; end

    # Recover wages and populations
    # y_i = w_i^(1-sigma) => w_i = y_i^(1/(1-sigma)) = y_i^(-1/(sigma-1))
    w = y .^ (-1 / (sigma - 1))
    w ./= geomean(w)  # numeraire normalization

    # x_i = w_i^sigma * L_i => L_i = x_i / w_i^sigma
    L = x ./ (w .^ sigma)
    L .*= (Lbar / sum(L))  # labor market clearing (eq 27)

    # Welfare: lambda = W^(sigma-1) => W = lambda^(1/(sigma-1))
    W = lambda ^ (1 / (sigma - 1))

    return (; w, L, W, lambda, x, y)
end

# ============================================================
# CASE 2: General Solver with Agglomeration
# Section 4.2, equations 28-32
# ============================================================

"""
    solve_general(Abar, ubar, tau, sigma, alpha, beta, Lbar; maxiter, damping, tol)

Solve the model with agglomeration and congestion forces using nested fixed-point
iteration on wages (inner loop) and population (outer loop).

Follows the solver pattern from Redding (2016) solveLw.jl.

Returns NamedTuple with w, L, tradesh, converged_w, converged_L, elapsed.
"""
function solve_general(Abar, ubar, tau, sigma, alpha, beta, Lbar;
                       maxiter=2000, damping=0.25, tol=6)
    tic = time()
    N = length(Abar)
    converged_w = false
    converged_L = false

    # Initialize
    L = fill(Lbar / N, N)
    w = ones(N)

    # Precompute tau^(1-sigma)
    tau_1ms = tau .^ (1 - sigma)

    println(">>>> Start Wage and Population Convergence <<<<")

    tradesh = zeros(N, N)
    nummat = zeros(N, N)

    for outer in 1:maxiter
        # Total productivity: A_i = Abar_i * L_i^alpha (eq 22)
        A = Abar .* L .^ alpha

        # ---- Inner loop: wage convergence ----
        converged_w = false
        for inner in 1:maxiter
            # Trade share numerator: tau_ij^(1-sigma) * (w_i / A_i)^(1-sigma)
            unit_cost_1ms = (w ./ A) .^ (1 - sigma)  # N-vector
            nummat .= tau_1ms .* (unit_cost_1ms * ones(1, N))  # N x N

            # Trade shares: pi_ij = num_ij / sum_k num_kj
            denom = sum(nummat, dims=1)  # 1 x N
            tradesh = nummat ./ denom     # pi_ij = share of j's expenditure on i's goods

            # Income and expenditure
            income = w .* L
            expend = tradesh * income  # expenditure on i's goods = sum_j pi_ij * w_j * L_j

            # Convergence check (round to tol decimal places)
            income_r = round.(income .* 10.0^tol)
            expend_r = round.(expend .* 10.0^tol)

            if income_r == expend_r
                converged_w = true
                break
            else
                # Wage update (following solveLw.jl pattern)
                w_e = w .* (expend ./ income) .^ (1 / sigma)
                w = damping .* w_e + (1 - damping) .* w
                w ./= geomean(w)
            end
        end

        # ---- Outer loop: population convergence ----
        # Compute welfare in each location: W_i = w_i * u_i / P_i (eq 18)
        # Price index: P_j^(1-sigma) = sum_i p_ij^(1-sigma) = sum_i nummat_ij (the denominator above)
        P_1ms = vec(sum(nummat, dims=1))  # P_j^(1-sigma)
        P = P_1ms .^ (1 / (1 - sigma))

        # Total amenities: u_i = ubar_i * L_i^beta (eq 23)
        u = ubar .* L .^ beta

        # Welfare: W_i = w_i * u_i / P_i (eq 18)
        W_i = w .* u ./ P

        # Implied population from welfare equalization
        # From eq 30: W = w_i * ubar_i * L_i^beta / P_i
        # Higher W_i => location i attracts more workers
        # Update L proportional to (W_i / W_mean)^(exponent) * L
        W_mean = geomean(W_i)
        # The relevant exponent comes from the power of L in the welfare expression
        # W_i propto L_i^(beta + alpha*(sigma-1)/sigma - ...), use a safe general update:
        L_e = L .* (W_i ./ W_mean) .^ (1.0 / max(abs(beta) + abs(alpha), 0.1))
        L_e .*= (Lbar / sum(L_e))

        # Convergence check
        L_r = round.(L .* 10.0^tol)
        L_e_r = round.(L_e .* 10.0^tol)

        if L_r == L_e_r
            println(">>>> Population Convergence Achieved at outer iteration $outer <<<<")
            converged_L = true
            L = L_e
            break
        else
            # Damped update
            L = damping .* L_e + (1 - damping) .* L
            L .*= (Lbar / sum(L))
        end
    end

    elapsed = time() - tic
    if converged_w
        println(">>>> Wage Convergence Achieved <<<<")
    else
        println("WARNING: Wage convergence NOT achieved")
    end
    if !converged_L
        println("WARNING: Population convergence NOT achieved")
    end
    println("Elapsed time: $(round(elapsed, digits=2))s")

    return (; w, L, tradesh, converged_w, converged_L, elapsed)
end

# ============================================================
# Derived Quantities
# ============================================================

"""
    compute_derived(w, L, Abar, ubar, tau, sigma, alpha, beta, Lbar)

Compute all derived equilibrium quantities from wages and populations.

Returns NamedTuple with A, u, P, X, tradesh, Pi_out, MA_in, W, W_i, Y_W, dtradesh.
"""
function compute_derived(w, L, Abar, ubar, tau, sigma, alpha, beta, Lbar)
    N = length(w)

    # Total productivity and amenity (eqs 22-23)
    A = Abar .* L .^ alpha
    u = ubar .* L .^ beta

    # Trade share matrix
    tau_1ms = tau .^ (1 - sigma)
    unit_cost_1ms = (w ./ A) .^ (1 - sigma)
    nummat = tau_1ms .* (unit_cost_1ms * ones(1, N))
    denom = sum(nummat, dims=1)
    tradesh = nummat ./ denom

    # Price indices: P_j^(1-sigma) = sum_i (tau_ij * w_i/A_i)^(1-sigma)
    P_1ms = vec(denom)
    P = P_1ms .^ (1 / (1 - sigma))

    # Trade flows: X_ij = pi_ij * w_j * L_j (eq 21 via trade shares)
    income = w .* L
    X = tradesh .* (ones(N) * income')

    # Domestic trade share (diagonal)
    dtradesh = diag(tradesh)

    # Outward market access: Pi_i^(1-sigma) = sum_j tau_ij^(1-sigma) * P_j^(sigma-1) * w_j * L_j (eq 34)
    P_sm1 = P_1ms .^ (-1)  # P_j^(sigma-1) = (P_j^(1-sigma))^(-1)
    Pi_out = [sum(tau_1ms[i,j] * P_sm1[j] * income[j] for j in 1:N) for i in 1:N]

    # Inward market access = P_i^(1-sigma) (already computed)
    MA_in = P_1ms

    # Welfare: W_i = w_i * u_i / P_i (eq 18)
    W_i = w .* u ./ P
    W = geomean(W_i)  # should be equalized across locations

    # World income
    Y_W = sum(income)

    # Labor demand curve (eq 35): ln w_i = -(1/sigma - alpha*(sigma-1)/sigma)*ln L_i + (1/sigma)*ln Pi_out_i + (sigma-1)/sigma * ln Abar_i
    ld_slope = -(1/sigma - alpha*(sigma-1)/sigma)
    ld_intercept = (1/sigma) .* log.(Pi_out) .+ ((sigma-1)/sigma) .* log.(Abar)

    # Labor supply curve (eq 36): ln w_i = -beta*ln L_i + 1/(1-sigma)*ln P_1ms_i + ln W - ln ubar_i
    ls_slope = -beta
    ls_intercept = (1/(1-sigma)) .* log.(P_1ms) .+ log(W) .- log.(ubar)

    return (; A, u, P, X, tradesh, Pi_out, MA_in, W, W_i, Y_W, dtradesh,
              ld_slope, ld_intercept, ls_slope, ls_intercept)
end

# ============================================================
# Counterfactual Analysis
# ============================================================

"""
    counterfactual(tau_new, Abar, ubar, sigma, alpha, beta, Lbar, baseline;
                   maxiter, damping, tol)

Solve the model under new trade costs and compare with baseline.

Returns NamedTuple with counterfactual results and welfare gains.
"""
function counterfactual(tau_new, Abar, ubar, sigma, alpha, beta, Lbar, baseline;
                        maxiter=2000, damping=0.25, tol=6)
    println("\n>>>> Counterfactual: Solving with new trade costs <<<<")
    cf = solve_general(Abar, ubar, tau_new, sigma, alpha, beta, Lbar;
                       maxiter=maxiter, damping=damping, tol=tol)
    cf_derived = compute_derived(cf.w, cf.L, Abar, ubar, tau_new, sigma, alpha, beta, Lbar)

    # Changes (hat variables)
    w_hat = cf.w ./ baseline.w
    L_hat = cf.L ./ baseline.L
    P_hat = cf_derived.P ./ baseline.P
    W_hat = cf_derived.W / baseline.W

    # Welfare gains
    println("\nWelfare change: W_hat = $(round(W_hat, digits=6))")

    return (; cf..., cf_derived..., w_hat, L_hat, P_hat, W_hat)
end

# ============================================================
# MAIN EXECUTION
# ============================================================

println("\n", "="^60)
println("CASE 1: No agglomeration (alpha=beta=0), eigenvector solution")
println("="^60)

result_eigen = solve_eigenvector(T_mat, sigma, Lbar)
derived_eigen = compute_derived(result_eigen.w, result_eigen.L, Abar, ubar, tau, sigma, 0.0, 0.0, Lbar)

@show round.(result_eigen.w, digits=4)
@show round.(result_eigen.L, digits=2)
@show round(result_eigen.W, digits=6)

# Verify welfare equalization
@show round.(derived_eigen.W_i, digits=6)
println("Max welfare deviation: $(round(maximum(abs.(derived_eigen.W_i .- derived_eigen.W)), digits=10))")

# Verify trade balance
income_eigen = result_eigen.w .* result_eigen.L
exports_eigen = vec(sum(derived_eigen.X, dims=2))
println("Max trade balance deviation: $(round(maximum(abs.(income_eigen .- exports_eigen)), digits=10))")

# Welfare elasticity (eq 42): -d ln W / d ln tau_ij = X_ij / Y_W
welfare_elasticity = derived_eigen.X ./ derived_eigen.Y_W
@show round.(welfare_elasticity, digits=4)

println("\n", "="^60)
println("CASE 2: With agglomeration (alpha=$alpha, beta=$beta)")
println("="^60)

result_gen = solve_general(Abar, ubar, tau, sigma, alpha, beta, Lbar)
derived_gen = compute_derived(result_gen.w, result_gen.L, Abar, ubar, tau, sigma, alpha, beta, Lbar)

@show round.(result_gen.w, digits=4)
@show round.(result_gen.L, digits=2)
@show round(derived_gen.W, digits=6)

# Verify welfare equalization
println("Max welfare deviation: $(round(maximum(derived_gen.W_i) - minimum(derived_gen.W_i), digits=8))")

# Verify trade balance
income_gen = result_gen.w .* result_gen.L
exports_gen = vec(sum(derived_gen.X, dims=2))
println("Max trade balance deviation: $(round(maximum(abs.(income_gen .- exports_gen)), digits=8))")

println("\n", "="^60)
println("COUNTERFACTUAL: 10% reduction in trade costs")
println("="^60)

# Reduce excess trade costs by 10%
tau_cf = 1.0 .+ 0.9 .* (tau .- 1.0)
cf_result = counterfactual(tau_cf, Abar, ubar, sigma, alpha, beta, Lbar,
                           (; w=result_gen.w, L=result_gen.L, P=derived_gen.P, W=derived_gen.W))

@show round.(cf_result.w_hat, digits=4)
@show round.(cf_result.L_hat, digits=4)
@show round(cf_result.W_hat, digits=6)

# ============================================================
# Plots
# ============================================================

# Population distribution: baseline vs counterfactual
p1 = bar(1:N, result_gen.L, label="Baseline", alpha=0.7,
    title="Population by Location", xlabel="Location", ylabel="Population",
    xticks=1:N)
bar!(1:N, cf_result.L, label="Counterfactual", alpha=0.5)

# Wages: baseline vs counterfactual
p2 = bar(1:N, result_gen.w, label="Baseline", alpha=0.7,
    title="Wages by Location", xlabel="Location", ylabel="Wage",
    xticks=1:N)
bar!(1:N, cf_result.w, label="Counterfactual", alpha=0.5)

# Welfare changes
p3 = bar(log.(cf_result.L_hat), label="Log population change",
    title="Log Changes from Trade Cost Reduction", xlabel="Location", ylabel="Log change",
    xticks=1:N)

# Labor supply and demand curves (in the general case)
ln_L = log.(result_gen.L)
ln_w = log.(result_gen.w)
p4 = scatter(ln_L, ln_w, label="Equilibrium", markersize=6,
    title="Labor Supply and Demand", xlabel="ln L", ylabel="ln w")
# Demand curve: ln w = ld_slope * ln L + ld_intercept
L_range = range(minimum(ln_L) - 0.5, stop=maximum(ln_L) + 0.5, length=100)
for i in 1:N
    plot!(p4, L_range, derived_gen.ld_slope .* L_range .+ derived_gen.ld_intercept[i],
        label=(i == 1 ? "Demand" : ""), color=:blue, alpha=0.3, linewidth=0.5)
    plot!(p4, L_range, derived_gen.ls_slope .* L_range .+ derived_gen.ls_intercept[i],
        label=(i == 1 ? "Supply" : ""), color=:red, alpha=0.3, linewidth=0.5)
end

plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 700))
savefig("QRE-HoRaUE-2025/simple_qre_model.pdf")
println("\nPlots saved to QRE-HoRaUE-2025/simple_qre_model.pdf")
