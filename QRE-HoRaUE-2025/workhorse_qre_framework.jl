# Workhorse Quantitative Economic Geography Framework
# Based on Section 6 of Allen & Arkolakis (2025)
# "Quantitative Regional Economics", Handbook of Regional and Urban Economics

# Implements: General gravity-based framework with four model elasticities,
# nesting Allen-Arkolakis, Krugman, Helpman, Eaton-Kortum, Redding, and other models.
# Includes equilibrium solver, exact hat algebra counterfactuals, and linear comparative statics.

using Random
using Statistics
using Plots
using LinearAlgebra
using StatsBase: geomean

# ============================================================
# Model Elasticity Mappings (Table 1)
# ============================================================

"""
    model_elasticities(model; sigma, alpha, beta, theta, mu, epsilon_pref)

Return the four workhorse elasticities (eps_S_local, eps_S_global, eps_D_local, eps_D_global)
for a named model from Table 1 of Allen & Arkolakis (2025).
"""
function model_elasticities(model::String; sigma=4.0, alpha=0.05, beta=-0.1,
                            theta=4.0, mu=0.75, epsilon_pref=3.0)
    if model == "Allen-Arkolakis"
        return (; eps_S_local = -beta, eps_S_global = 1/(sigma-1),
                  eps_D_local = (1-alpha*(sigma-1))/sigma, eps_D_global = 1/sigma)
    elseif model == "Rosen-Roback"
        return (; eps_S_local = -beta, eps_S_global = 0.0,
                  eps_D_local = alpha, eps_D_global = 0.0)
    elseif model == "Krugman"
        return (; eps_S_local = 0.0, eps_S_global = 1/(sigma-1),
                  eps_D_local = 0.0, eps_D_global = 1/sigma)
    elseif model == "Helpman"
        return (; eps_S_local = (1-mu)/mu, eps_S_global = 1/(sigma-1),
                  eps_D_local = 0.0, eps_D_global = 1/sigma)
    elseif model == "Eaton-Kortum"
        return (; eps_S_local = 0.0, eps_S_global = 1/theta,
                  eps_D_local = 1/(1+theta), eps_D_global = 1/(1+theta))
    elseif model == "Donaldson-Hornbeck"
        gamma = 1 - alpha - mu  # labor share (residual)
        return (; eps_S_local = 0.0, eps_S_global = 1/theta,
                  eps_D_local = (1+alpha*theta)/(1+(alpha+gamma)*theta),
                  eps_D_global = 1/(1+(alpha+gamma)*theta))
    elseif model == "Redding-CRS"
        return (; eps_S_local = (1+(1-mu)*epsilon_pref)/(mu*epsilon_pref),
                  eps_S_global = 1/theta,
                  eps_D_local = 1/(1+theta), eps_D_global = 1/(1+theta))
    elseif model == "Redding-IRS"
        return (; eps_S_local = (1+(1-mu)*epsilon_pref)/(mu*epsilon_pref),
                  eps_S_global = 1/(sigma-1),
                  eps_D_local = 0.0, eps_D_global = 1/sigma)
    else
        error("Unknown model: $model. Choose from: Allen-Arkolakis, Rosen-Roback, Krugman, Helpman, Eaton-Kortum, Donaldson-Hornbeck, Redding-CRS, Redding-IRS")
    end
end

# ============================================================
# A and B Matrices (Section 6.3-6.5)
# ============================================================

"""
    compute_AB_matrices(eps_S_local, eps_S_global, eps_D_local, eps_D_global)

Compute the B matrix (change of variables, eqs 86-87), A matrix (model elasticities,
eqs 71-72), and check uniqueness condition rho(|A|) <= 1.

Returns NamedTuple with A, B, B_inv, rho_A, unique.
"""
function compute_AB_matrices(eps_S_local, eps_S_global, eps_D_local, eps_D_global)
    # B matrix: maps (ln y, ln l) -> (ln x, ln z) from eqs 86-87
    # ln x = (1/eps_D_global)*ln y - ((1-eps_D_local)/eps_D_global)*ln l     (eq 86)
    # ln z = (-1/eps_S_global)*ln y + ((eps_S_local+1)/eps_S_global)*ln l    (eq 87)
    B = [1/eps_D_global       (eps_D_local-1)/eps_D_global;
        -1/eps_S_global       (eps_S_local+1)/eps_S_global]

    # Check if B is singular (some models like Krugman with eps_D_local=0 have det(B)=0)
    d = det(B)
    if abs(d) < 1e-12
        # B is singular; cannot invert. Return with NaN A matrix.
        A = fill(NaN, 2, 2)
        return (; A, B, B_inv=fill(NaN, 2, 2), rho_A=NaN, unique=false, singular_B=true)
    end

    B_inv = inv(B)

    # RHS exponent matrix from eqs 68-69 (exponents on y_j and l_j on the right-hand side)
    RHS = [(eps_S_global+1)/eps_S_global    -(1+eps_S_local)/eps_S_global;
           (eps_D_global-1)/eps_D_global     (1-eps_D_local)/eps_D_global]

    # A matrix: A = B^(-1) * RHS (eqs 71-72)
    A = B_inv * RHS

    # Uniqueness: spectral radius of |A| (Section 6.5)
    A_abs = abs.(A)
    rho_A = maximum(abs.(eigvals(A_abs)))

    is_unique = rho_A <= 1.0 + 1e-10  # small tolerance for floating point

    return (; A, B, B_inv, rho_A, unique=is_unique, singular_B=false)
end

# ============================================================
# Equilibrium Solver (eqs 68-69 / 71-72)
# ============================================================

"""
    solve_workhorse(K, eps_S_local, eps_S_global, eps_D_local, eps_D_global, Lbar;
                    maxiter, damping, tol)

Solve the workhorse equilibrium system (eqs 71-72) using fixed-point iteration.

K: N x N composite geography matrix, K_ij = T_ij * C_D_i^(1/eps_D_global) * C_S_j^(1/eps_S_global)

Returns NamedTuple with w, L, y, l, x, z, lambda, MA_out, MA_in, X, Y, E, A_mat, B_mat, converged.
"""
function solve_workhorse(K, eps_S_local, eps_S_global, eps_D_local, eps_D_global, Lbar;
                         maxiter=5000, damping=0.25, tol=6)
    tic = time()
    N = size(K, 1)

    # Compute A and B matrices
    ab = compute_AB_matrices(eps_S_local, eps_S_global, eps_D_local, eps_D_global)
    A_mat = ab.A
    B_mat = ab.B
    B_inv = ab.B_inv
    a11, a12 = A_mat[1,1], A_mat[1,2]
    a21, a22 = A_mat[2,1], A_mat[2,2]

    println("A matrix of model elasticities:")
    display(round.(A_mat, digits=4))
    println("Spectral radius rho(|A|) = $(round(ab.rho_A, digits=4))")
    ab.unique ? println("Unique equilibrium guaranteed.") : println("WARNING: rho(|A|) > 1, uniqueness not guaranteed!")

    # Initialize x, z (positive vectors)
    x = ones(N)
    z = ones(N)
    converged = false

    println(">>>> Start Equilibrium Convergence <<<<")

    for iter in 1:maxiter
        # Eq 71: lambda * x_i = sum_j K_ij * x_j^a11 * z_j^a12
        rhs_x = K * (x .^ a11 .* z .^ a12)

        # Eq 72: lambda * z_i = sum_j K_ji * x_j^a21 * z_j^a22
        rhs_z = K' * (x .^ a21 .* z .^ a22)

        # Compute lambda as the geometric mean of implied lambdas
        lambda_x = sum(rhs_x) / sum(x)
        lambda_z = sum(rhs_z) / sum(z)
        lambda = sqrt(lambda_x * lambda_z)

        # Update
        x_new = rhs_x ./ lambda
        z_new = rhs_z ./ lambda

        # Normalize to prevent drift
        x_new ./= geomean(x_new)
        z_new ./= geomean(z_new)

        # Damped update
        x_new = damping .* x_new + (1 - damping) .* x
        z_new = damping .* z_new + (1 - damping) .* z

        # Convergence check
        x_r = round.(x_new .* 10.0^tol)
        z_r = round.(z_new .* 10.0^tol)
        x_old_r = round.(x .* 10.0^tol)
        z_old_r = round.(z .* 10.0^tol)

        if x_r == x_old_r && z_r == z_old_r
            println(">>>> Equilibrium Convergence Achieved at iteration $iter <<<<")
            converged = true
            x = x_new
            z = z_new
            break
        end

        x = x_new
        z = z_new
    end

    if !converged
        println("WARNING: Equilibrium convergence NOT achieved after $maxiter iterations")
    end

    # Recover (y, l) from (x, z) using B^(-1)
    # [ln y; ln l] = B^(-1) * [ln x; ln z]
    ln_x = log.(x)
    ln_z = log.(z)
    ln_y = B_inv[1,1] .* ln_x + B_inv[1,2] .* ln_z
    ln_l = B_inv[2,1] .* ln_x + B_inv[2,2] .* ln_z

    y = exp.(ln_y)
    l = exp.(ln_l)

    # Normalize shares to sum to 1
    y ./= sum(y)
    l ./= sum(l)

    # Recover physical quantities
    L = l .* Lbar
    # w_i * L_i / Y_W = y_i, and Y_W is just a scale factor
    # w_i = y_i * Y_W / L_i = y_i / l_i * (Y_W / Lbar)
    # For now set Y_W = Lbar (numeraire), so w_i = y_i / l_i
    w = y ./ l
    w ./= geomean(w)

    # Recalculate y consistently
    Y = w .* L
    Y_W = sum(Y)
    y = Y ./ Y_W

    # Compute lambda from the equilibrium
    # Recompute final lambda
    rhs_x_final = K * (x .^ a11 .* z .^ a12)
    lambda = sum(rhs_x_final) / sum(x)

    # Compute market access (eqs 64-65)
    E = Y  # balanced trade: E_i = Y_i
    MA_out, MA_in = solve_market_access(K, Y, E, eps_S_global, eps_D_global;
                                         maxiter=maxiter, damping=damping, tol=tol)

    # Trade flows (eq 61): X_ij = T_ij * Y_i / MA_out_i * E_j / MA_in_j
    # Note: K_ij includes C_D and C_S, but trade flows use T_ij.
    # Since we only have K, we compute X_ij = K_ij * Y_i / MA_out_i * E_j / MA_in_j
    # (proportional to true trade flows, fine for shares)
    X = zeros(N, N)
    for i in 1:N, j in 1:N
        X[i,j] = K[i,j] * (Y[i] / MA_out[i]) * (E[j] / MA_in[j])
    end
    # Normalize rows to match Y_i
    for i in 1:N
        X[i,:] .*= Y[i] / sum(X[i,:])
    end

    elapsed = time() - tic
    println("Elapsed time: $(round(elapsed, digits=2))s")

    return (; w, L, y, l, x, z, lambda, MA_out, MA_in, X, Y, E, Y_W,
              A_mat, B_mat, B_inv, converged, elapsed)
end

"""
    solve_market_access(K, Y, E, eps_S_global, eps_D_global; maxiter, damping, tol)

Solve the market access equations (64-65) by fixed-point iteration.
MA_out_i = sum_j K_ij * E_j / MA_in_j
MA_in_i = sum_j K_ji * Y_j / MA_out_j
"""
function solve_market_access(K, Y, E, eps_S_global, eps_D_global;
                              maxiter=2000, damping=0.25, tol=6)
    N = length(Y)
    MA_out = ones(N)
    MA_in = ones(N)

    for iter in 1:maxiter
        MA_out_new = [sum(K[i,j] * E[j] / MA_in[j] for j in 1:N) for i in 1:N]
        MA_in_new = [sum(K[j,i] * Y[j] / MA_out[j] for j in 1:N) for i in 1:N]

        # Normalize
        MA_out_new ./= geomean(MA_out_new)
        MA_in_new ./= geomean(MA_in_new)

        # Convergence check
        out_r = round.(MA_out_new .* 10.0^tol)
        in_r = round.(MA_in_new .* 10.0^tol)
        out_old_r = round.(MA_out .* 10.0^tol)
        in_old_r = round.(MA_in .* 10.0^tol)

        if out_r == out_old_r && in_r == in_old_r
            MA_out = MA_out_new
            MA_in = MA_in_new
            break
        end

        MA_out = damping .* MA_out_new + (1 - damping) .* MA_out
        MA_in = damping .* MA_in_new + (1 - damping) .* MA_in
    end

    return MA_out, MA_in
end

# ============================================================
# Counterfactual: Exact Hat Algebra (eqs 80-81)
# ============================================================

"""
    counterfactual_hat(K_hat, result; maxiter, damping, tol)

Perform counterfactual analysis using exact hat algebra (eqs 80-81).

K_hat: N x N matrix of changes in composite geography (K_ij^B / K_ij^A).
result: NamedTuple from solve_workhorse (baseline equilibrium).

Returns NamedTuple with y_hat, l_hat, w_hat, L_hat, lambda_hat.
"""
function counterfactual_hat(K_hat, result;
                            maxiter=5000, damping=0.25, tol=6)
    N = length(result.w)
    eps_S_local = result.B_mat[2,2] * result.B_mat[2,2] # recover from stored matrices
    # Actually, we need the elasticities. Extract from A and B matrices.
    # Better: pass them explicitly or recover from the B matrix
    # B = [1/eD_g, (eD_l-1)/eD_g; -1/eS_g, (eS_l+1)/eS_g]
    eD_g = 1 / result.B_mat[1,1]
    eD_l = result.B_mat[1,2] * eD_g + 1
    eS_g = -1 / result.B_mat[2,1]
    eS_l = result.B_mat[2,2] * eS_g - 1

    # Export shares: T_ij = X_ij / Y_i
    T_share = result.X ./ result.Y
    # Import shares: S_ij = X_ij / E_j  => S_ji = X_ji / E_i
    # For eq 81: sum_j (X_ji^A / E_i^A) * K_hat_ji * ...
    # D_ji = X_ji / E_i
    E = result.E

    # Exponents from eqs 80-81
    # LHS of eq 80: lambda_hat * y_hat_i^(1/eD_g) * l_hat_i^((eD_l-1)/eD_g)
    # RHS of eq 80: sum_j (X_ij/Y_i) * K_hat_ij * y_hat_j^((eS_g+1)/eS_g) * l_hat_j^(-(1+eS_l)/eS_g)
    # LHS of eq 81: lambda_hat * y_hat_i^(-1/eS_g) * l_hat_i^((eS_l+1)/eS_g)
    # RHS of eq 81: sum_j (X_ji/E_i) * K_hat_ji * y_hat_j^((eD_g-1)/eD_g) * l_hat_j^((1-eD_l)/eD_g)

    # Use the (x_hat, z_hat) formulation via A matrix
    A_mat = result.A_mat
    B_inv = result.B_inv
    a11, a12 = A_mat[1,1], A_mat[1,2]
    a21, a22 = A_mat[2,1], A_mat[2,2]

    # Weighted K_hat matrices using export and import shares
    # For eq 80 (export side): C_ij = (X_ij/Y_i) * K_hat_ij
    C = T_share .* K_hat
    # For eq 81 (import side): D_ji = (X_ji/E_i) * K_hat_ji
    # D[j,i] = X[j,i]/E[i] * K_hat[j,i]
    S_by_E = result.X ./ (ones(N) * E')  # S[j,i] = X[j,i] / E[i]  -- this is N x N
    D = S_by_E .* K_hat

    # Initialize
    x_hat = ones(N)
    z_hat = ones(N)
    converged = false

    y_A = result.y
    l_A = result.l

    for iter in 1:maxiter
        # RHS of eq 80 in (x,z) form: sum_j C_ij * x_hat_j^a11 * z_hat_j^a12
        rhs_x = C * (x_hat .^ a11 .* z_hat .^ a12)

        # RHS of eq 81 in (x,z) form: sum_j D_ji * x_hat_j^a21 * z_hat_j^a22
        # D is stored as D[j,i], so D' gives us D[i,j] indexed by i
        rhs_z = D' * (x_hat .^ a21 .* z_hat .^ a22)

        # Compute lambda_hat from adding-up constraints:
        # sum(y_A .* y_hat) = 1 and sum(l_A .* l_hat) = 1
        # y_hat and l_hat are recovered from x_hat, z_hat via B_inv
        # For now, use normalization on x_hat and z_hat
        lambda_x = sum(rhs_x .* x_hat) / sum(x_hat .^ 2)
        lambda_z = sum(rhs_z .* z_hat) / sum(z_hat .^ 2)
        lambda_hat = sqrt(abs(lambda_x * lambda_z))
        if lambda_hat == 0; lambda_hat = 1.0; end

        x_hat_new = rhs_x ./ lambda_hat
        z_hat_new = rhs_z ./ lambda_hat

        # Recover y_hat, l_hat to enforce adding-up constraints
        ln_xh = log.(max.(x_hat_new, 1e-15))
        ln_zh = log.(max.(z_hat_new, 1e-15))
        ln_yh = B_inv[1,1] .* ln_xh + B_inv[1,2] .* ln_zh
        ln_lh = B_inv[2,1] .* ln_xh + B_inv[2,2] .* ln_zh
        y_hat_raw = exp.(ln_yh)
        l_hat_raw = exp.(ln_lh)

        # Enforce: sum(y_A .* y_hat) = 1
        y_hat_raw ./= sum(y_A .* y_hat_raw)
        l_hat_raw ./= sum(l_A .* l_hat_raw)

        # Back to (x_hat, z_hat)
        ln_yh = log.(y_hat_raw)
        ln_lh = log.(l_hat_raw)
        x_hat_new = exp.(result.B_mat[1,1] .* ln_yh + result.B_mat[1,2] .* ln_lh)
        z_hat_new = exp.(result.B_mat[2,1] .* ln_yh + result.B_mat[2,2] .* ln_lh)

        # Damped update
        x_hat_new = damping .* x_hat_new + (1 - damping) .* x_hat
        z_hat_new = damping .* z_hat_new + (1 - damping) .* z_hat

        # Convergence
        if round.(x_hat_new .* 10.0^tol) == round.(x_hat .* 10.0^tol) &&
           round.(z_hat_new .* 10.0^tol) == round.(z_hat .* 10.0^tol)
            println(">>>> Counterfactual Convergence at iteration $iter <<<<")
            converged = true
            x_hat = x_hat_new
            z_hat = z_hat_new
            break
        end

        x_hat = x_hat_new
        z_hat = z_hat_new
    end

    if !converged
        println("WARNING: Counterfactual convergence NOT achieved")
    end

    # Final recovery
    ln_xh = log.(x_hat)
    ln_zh = log.(z_hat)
    ln_yh = B_inv[1,1] .* ln_xh + B_inv[1,2] .* ln_zh
    ln_lh = B_inv[2,1] .* ln_xh + B_inv[2,2] .* ln_zh
    y_hat = exp.(ln_yh)
    l_hat = exp.(ln_lh)

    # Normalize
    y_hat ./= sum(y_A .* y_hat)
    l_hat ./= sum(l_A .* l_hat)

    # Welfare change
    lambda_hat = sum(C * (x_hat .^ a11 .* z_hat .^ a12) ./ x_hat) / N

    # Recover wage and population changes
    w_hat = (y_hat ./ l_hat) ./ (y_A ./ l_A) .* (y_A ./ l_A)
    # More precisely: w_new*L_new/Y_W_new = y_A*y_hat and L_new/Lbar = l_A*l_hat
    # => w_hat = y_hat / l_hat (up to aggregate scale)
    w_hat = y_hat ./ l_hat
    L_hat = l_hat

    return (; y_hat, l_hat, w_hat, L_hat, lambda_hat, converged, x_hat, z_hat)
end

# ============================================================
# Comparative Statics (eq 93)
# ============================================================

"""
    comparative_statics(result, dln_K; eps_S_local, eps_S_global, eps_D_local, eps_D_global)

Compute first-order comparative statics for small changes in geography (eq 93).

result: NamedTuple from solve_workhorse.
dln_K: N x N matrix of log changes in K_ij.

Returns NamedTuple with dln_x, dln_z, dln_y, dln_l, dln_lambda, dln_w.
"""
function comparative_statics(result, dln_K)
    N = length(result.w)
    A_mat = result.A_mat
    B_mat = result.B_mat
    B_inv = result.B_inv
    a11, a12 = A_mat[1,1], A_mat[1,2]
    a21, a22 = A_mat[2,1], A_mat[2,2]
    b11, b12 = B_inv[1,1], B_inv[1,2]  # Note: using B_inv for recovery
    b21, b22 = B_inv[2,1], B_inv[2,2]

    # From B_mat for Q construction
    bm11, bm12 = B_mat[1,1], B_mat[1,2]
    bm21, bm22 = B_mat[2,1], B_mat[2,2]

    # Export share matrix: T[i,j] = X_ij / Y_i
    T_share = result.X ./ result.Y
    # Import share matrix: S[i,j] = X_ij / E_j
    S_share = result.X ./ (ones(N) * result.E')

    y = result.y
    l = result.l

    # Q matrix (2N x 2N) - ensures adding-up constraints (eq after 92)
    # Q = [b11/(b11+b12) * ones(N)*y'   b12/(b11+b12) * ones(N)*y' ;
    #      b21/(b21+b22) * ones(N)*l'   b22/(b21+b22) * ones(N)*l' ]
    # But here b's refer to B_mat (not B_inv). Let me re-read the article...
    # From eq 88-90: the B used for Q is actually the B_mat (the original B, not inverse)
    # bm = B_mat entries
    Q_tl = (bm11 / (bm11 + bm12)) .* ones(N) * y'
    Q_tr = (bm12 / (bm11 + bm12)) .* ones(N) * y'
    Q_bl = (bm21 / (bm21 + bm22)) .* ones(N) * l'
    Q_br = (bm22 / (bm21 + bm22)) .* ones(N) * l'
    Q = [Q_tl Q_tr; Q_bl Q_br]

    # Block coefficient matrix (2N x 2N)
    M = [a11 .* T_share   a12 .* T_share;
         a21 .* S_share'   a22 .* S_share']

    # Shock vectors
    # diag(T * dln_K') = sum over columns: T_share .* dln_K' row sums
    shock_x = vec(sum(T_share .* dln_K, dims=2))   # diag(T * dln_K^T) when dln_K is symmetric
    shock_z = vec(sum(S_share' .* dln_K', dims=2))  # diag(S^T * dln_K)
    shock = [shock_x; shock_z]

    # Solve eq 93
    I_2N = Matrix{Float64}(I, 2N, 2N)
    coeff = I_2N - (I_2N - Q) * M
    rhs = (I_2N - Q) * shock
    dln_xz = coeff \ rhs

    dln_x = dln_xz[1:N]
    dln_z = dln_xz[N+1:2N]

    # Recover dln_y and dln_l using B^(-1)
    dln_y = B_inv[1,1] .* dln_x + B_inv[1,2] .* dln_z
    dln_l = B_inv[2,1] .* dln_x + B_inv[2,2] .* dln_z

    # Wage change: dln_w = dln_y - dln_l (since y_i = w_i*L_i/Y_W and l_i = L_i/Lbar)
    dln_w = dln_y - dln_l

    # Welfare change (eq 91)
    bm_sum1 = bm11 + bm12
    dln_lambda = (y' * (bm11 .* (shock_x + T_share * (a11 .* dln_x + a12 .* dln_z)) +
                         bm12 .* (shock_z + S_share' * (a21 .* dln_x + a22 .* dln_z)))) / bm_sum1

    return (; dln_x, dln_z, dln_y, dln_l, dln_w, dln_lambda=dln_lambda[1])
end

# ============================================================
# Neumann Series Decomposition (eq 94)
# ============================================================

"""
    neumann_decomposition(result, dln_K; K_max=10)

Compute the Neumann series expansion of the comparative statics (eq 94),
decomposing the total effect into direct (k=0) and indirect (k=1,...) effects.

Returns vector of 2N-vectors, one per degree k.
"""
function neumann_decomposition(result, dln_K; K_max=10)
    N = length(result.w)
    A_mat = result.A_mat
    B_mat = result.B_mat
    a11, a12 = A_mat[1,1], A_mat[1,2]
    a21, a22 = A_mat[2,1], A_mat[2,2]
    bm11, bm12 = B_mat[1,1], B_mat[1,2]
    bm21, bm22 = B_mat[2,1], B_mat[2,2]

    T_share = result.X ./ result.Y
    S_share = result.X ./ (ones(N) * result.E')
    y = result.y
    l = result.l

    Q_tl = (bm11 / (bm11 + bm12)) .* ones(N) * y'
    Q_tr = (bm12 / (bm11 + bm12)) .* ones(N) * y'
    Q_bl = (bm21 / (bm21 + bm22)) .* ones(N) * l'
    Q_br = (bm22 / (bm21 + bm22)) .* ones(N) * l'
    Q = [Q_tl Q_tr; Q_bl Q_br]

    M = [a11 .* T_share   a12 .* T_share;
         a21 .* S_share'   a22 .* S_share']

    shock_x = vec(sum(T_share .* dln_K, dims=2))
    shock_z = vec(sum(S_share' .* dln_K', dims=2))
    shock = [shock_x; shock_z]

    I_2N = Matrix{Float64}(I, 2N, 2N)
    IQ = I_2N - Q
    IQ_M = IQ * M

    contributions = Vector{Vector{Float64}}()
    base = IQ * shock
    push!(contributions, copy(base))

    current = base
    for k in 1:K_max
        current = IQ_M * current
        push!(contributions, copy(current))
    end

    return contributions
end

# ============================================================
# MAIN EXECUTION
# ============================================================

# ****************************
# **** Parameters ****
# ****************************

N = 10
Lbar = 1000.0

# Default: Allen-Arkolakis (2014) parameterization
sigma = 4.0
alpha_param = 0.05
beta_param = -0.1

elast = model_elasticities("Allen-Arkolakis"; sigma=sigma, alpha=alpha_param, beta=beta_param)
eps_S_local = elast.eps_S_local
eps_S_global = elast.eps_S_global
eps_D_local = elast.eps_D_local
eps_D_global = elast.eps_D_global

println("Model elasticities (Allen-Arkolakis):")
println("  eps_S_local  = $(round(eps_S_local, digits=4))")
println("  eps_S_global = $(round(eps_S_global, digits=4))")
println("  eps_D_local  = $(round(eps_D_local, digits=4))")
println("  eps_D_global = $(round(eps_D_global, digits=4))")

# ****************************
# **** Geography ****
# ****************************

Random.seed!(1)

# Location-specific supply and demand shifters
C_S = exp.(randn(N)); C_S ./= geomean(C_S)
C_D = exp.(randn(N)); C_D ./= geomean(C_D)

# Trade friction matrix T_ij <= 1 (distance-based)
coords = range(0, stop=4, length=N)
dist_mat = [abs(coords[i] - coords[j]) for i in 1:N, j in 1:N]
T_ij = exp.(-0.3 .* dist_mat)
for i in 1:N; T_ij[i,i] = 1.0; end

# Composite geography: K_ij = T_ij * C_D_i^(1/eps_D_global) * C_S_j^(1/eps_S_global)
K = [T_ij[i,j] * C_D[i]^(1/eps_D_global) * C_S[j]^(1/eps_S_global) for i in 1:N, j in 1:N]

# ============================================================
# Solve Equilibrium
# ============================================================

println("\n", "="^60)
println("Solving Workhorse Framework (Allen-Arkolakis parameterization)")
println("="^60)

result = solve_workhorse(K, eps_S_local, eps_S_global, eps_D_local, eps_D_global, Lbar)

println("\nEquilibrium results:")
@show round.(result.w, digits=4)
@show round.(result.L, digits=2)
@show round.(result.y, digits=4)
@show round.(result.l, digits=4)
@show round(result.lambda, digits=6)

# Verify trade balance
exports = vec(sum(result.X, dims=2))
imports = vec(sum(result.X, dims=1))
println("Max trade imbalance (exports): $(round(maximum(abs.(exports .- result.Y)), digits=8))")
println("Max trade imbalance (imports): $(round(maximum(abs.(imports .- result.E)), digits=8))")

# ============================================================
# Counterfactual: 10% reduction in trade frictions
# ============================================================

println("\n", "="^60)
println("Counterfactual: 10% increase in off-diagonal K (trade liberalization)")
println("="^60)

K_hat = ones(N, N)
K_hat .+= 0.1 .* (1.0 .- Matrix{Float64}(I, N, N))  # 10% increase in off-diagonal

cf = counterfactual_hat(K_hat, result)

println("\nCounterfactual changes:")
@show round.(cf.y_hat, digits=4)
@show round.(cf.l_hat, digits=4)
@show round.(cf.w_hat, digits=4)
@show round(cf.lambda_hat, digits=6)

# ============================================================
# Comparative Statics (small shock)
# ============================================================

println("\n", "="^60)
println("Comparative Statics: 1% increase in off-diagonal K")
println("="^60)

dln_K = 0.01 .* (1.0 .- Matrix{Float64}(I, N, N))
cs = comparative_statics(result, dln_K)

println("\nFirst-order changes:")
@show round.(cs.dln_y, digits=6)
@show round.(cs.dln_l, digits=6)
@show round.(cs.dln_w, digits=6)
@show round(cs.dln_lambda, digits=6)

# Compare with small exact counterfactual
K_hat_small = ones(N, N) .+ dln_K
cf_small = counterfactual_hat(K_hat_small, result)

println("\nComparison (linear approx vs exact for 1% shock):")
println("  dln_y (linear):  $(round.(cs.dln_y, digits=6))")
println("  dln_y (exact):   $(round.(log.(cf_small.y_hat), digits=6))")
println("  Max |diff|:      $(round(maximum(abs.(cs.dln_y .- log.(cf_small.y_hat))), digits=8))")

# ============================================================
# Neumann Series Decomposition
# ============================================================

println("\n", "="^60)
println("Neumann Series Decomposition (eq 94)")
println("="^60)

contributions = neumann_decomposition(result, dln_K; K_max=10)
cumulative = zeros(2N)
println("Degree | Cumulative norm (first N components)")
for (k, c) in enumerate(contributions)
    cumulative .+= c
    println("  k=$(k-1)   | $(round(norm(cumulative[1:N]), digits=8))")
end
println("  Exact  | $(round(norm(cs.dln_x), digits=8))")

# ============================================================
# Model Comparison
# ============================================================

println("\n", "="^60)
println("Model Comparison: Uniqueness conditions across models")
println("="^60)

for model_name in ["Allen-Arkolakis", "Krugman", "Helpman", "Eaton-Kortum", "Redding-CRS", "Redding-IRS"]
    e = model_elasticities(model_name; sigma=sigma, alpha=alpha_param, beta=beta_param)
    if e.eps_S_global == 0.0 || e.eps_D_global == 0.0
        println("  $model_name: No spatial linkages (trivial equilibrium)")
        continue
    end
    ab = compute_AB_matrices(e.eps_S_local, e.eps_S_global, e.eps_D_local, e.eps_D_global)
    if ab.singular_B
        println("  $model_name: B matrix is singular (degenerate change of variables)")
    else
        println("  $model_name: rho(|A|) = $(round(ab.rho_A, digits=4)) $(ab.unique ? "(unique)" : "(WARNING: not unique)")")
    end
end

# ============================================================
# Plots
# ============================================================

# Population and income shares
p1 = bar(1:N, result.y, label="Income share y", alpha=0.7,
    title="Equilibrium Shares", xlabel="Location", ylabel="Share",
    xticks=1:N)
bar!(1:N, result.l, label="Population share l", alpha=0.5)

# Counterfactual changes
p2 = bar(1:N, log.(cf.y_hat), label="d ln y", alpha=0.7,
    title="Counterfactual (10% trade liberalization)", xlabel="Location", ylabel="Log change",
    xticks=1:N)
bar!(1:N, log.(cf.l_hat), label="d ln l", alpha=0.5)

# Linear vs exact comparison
p3 = scatter(cs.dln_l, log.(cf_small.l_hat),
    xlabel="Linear approximation", ylabel="Exact (log)",
    title="Comparative statics vs exact (1% shock)",
    label="d ln l", markersize=5)
plot!(p3, [-0.01, 0.01], [-0.01, 0.01], label="45-degree line", linestyle=:dash, color=:gray)

# Neumann convergence
norms = [norm(sum(contributions[1:k])[1:N]) for k in 1:length(contributions)]
p4 = plot(0:length(norms)-1, norms,
    xlabel="Degree k", ylabel="Cumulative norm",
    title="Neumann Series Convergence",
    marker=:circle, label="Cumulative effect")
hline!([norm(cs.dln_x)], label="Exact", linestyle=:dash, color=:red)

plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 700))
savefig("QRE-HoRaUE-2025/workhorse_qre_framework.pdf")
println("\nPlots saved to QRE-HoRaUE-2025/workhorse_qre_framework.pdf")
