# =====================================================================================
# The Economics of Density: a Julia port of the ARSW (2015) quantitative urban model
# =====================================================================================
#
# Ahlfeldt, G. M., Redding, S. J., Sturm, D. M., Wolf, N. (2015),
# "The Economics of Density: Evidence from the Berlin Wall",
# Econometrica 83(6), 2127-2189.   https://doi.org/10.3982/ECTA10876
#
# Ported from Gabriel M. Ahlfeldt's didactic MATLAB toolkit (ARSW2015-toolkit), whose 16
# numbered algorithms map ~1:1 onto MATLAB files. Every Julia function below names the
# MATLAB file and the supplement equation it implements.
#
# -------------------------------------------------------------------------------------
# RUNS ON SYNTHETIC DATA.
# -------------------------------------------------------------------------------------
# The toolkit needs three .mat files (prepdata_big_TD.mat, prepdata_big_TD86.mat,
# ttpublic_2006_ren.mat) that are external downloads and are NOT present. This file
# therefore generates a synthetic Berlin-like block geography, complete with a "wall"
# that severs commuting between the two halves of the city, plants known fundamentals,
# solves the model FORWARD to manufacture "observed" data, and then runs the toolkit's
# INVERSION on that data to check the planted fundamentals come back.
#
# That is a sharper test than matching Berlin's numbers. Supplement Prop. S.3 (sec.
# S.3.1.6) proves a one-to-one map {alpha,beta,mu,eps,kappa} + {Q,H_M,H_R,K,tau} -> {A,B,phi};
# Prop. S.4 (sec. S.3.2) extends it with {lambda,delta,eta,rho} to {a,b}. If the port is
# correct, recovery must hold to machine precision. It does (see the validation block).
#
# None of the paper's EMPIRICAL results are reproduced or claimed here: lambda = 0.0710,
# eps = 6.83/6.694, nu = 0.07/0.0987 are used as parameter INPUTS to the synthetic
# exercise, never as targets. To run on the real Berlin data, implement the single
# clearly marked function `load_berlin_data` (see its docstring for the schema).
#
# =====================================================================================
# THE MODEL
# =====================================================================================
# N blocks; each block is simultaneously a candidate workplace j and residence i.
# Workers draw i.i.d. Frechet(eps) shocks over (residence, workplace) pairs.
#
#   Commuting / residence choice
#     (4)     phi_ij = exp(-eps*kappa*tau_ij) * B_i^eps * Q_i^(-(1-beta)*eps) * w_j^eps
#             pi_ij  = phi_ij / Phi,   Phi = sum_ij phi_ij
#     (5)     H_Ri = (sum_j phi_ij / Phi) * H ;   H_Mj = (sum_i phi_ij / Phi) * H
#     (6)     pi_{j|i} = phi_ij / sum_j' phi_ij'
#     (S.20)  E[w|i] = sum_j pi_{j|i} w_j ;   total worker income vv_i = E[w|i] * H_Ri
#     (9)     Ubar = gamma * Phi^(1/eps),  gamma = Gamma((eps-1)/eps)
#
#   Production and the floor-space market
#     (10)    Y_j = A_j * H_Mj^alpha * (theta_j L_j)^(1-alpha)
#     (12)    q_j = (1-alpha) (alpha/w_j)^(alpha/(1-alpha)) A_j^(1/(1-alpha))
#     (13)    no arbitrage: theta_j = 1 if q>xi*Q, in [0,1] if q=xi*Q, 0 if q<xi*Q
#     (18)/(S.30)  theta_j L_j     = ((1-alpha) A_j / q_j)^(1/alpha) * H_Mj
#     (19)/(S.29)  (1-theta_i) L_i = (1-beta) E[w|i] H_Ri / Q_i
#     (15)/(S.31)  L_i = phi_i * K_i^(1-mu)
#
#   Inversion (this is what the toolkit's Algorithms 4-9 actually solve)
#     (S.44)  H_Mj = sum_i [ omega_j e^{-nu tau_ij} / sum_s omega_s e^{-nu tau_is} ] H_Ri,
#             omega_j = w_j^eps,  nu = kappa*eps     <- omega is identified WITHOUT eps
#     (S.46)  W_i = sum_s (w_s / e^{kappa tau_is})^eps      (commuting market access, CMA)
#     (S.47)  B_i/Bbar = (H_Ri/Hbar_R)^(1/eps) (Q_i/Qbar)^(1-beta) (W_i/Wbar)^(-1/eps)
#     (S.48)  A_j = (q_j/(1-alpha))^(1-alpha) (w_j/alpha)^alpha
#
#   Endogenous agglomeration (Section 7)
#     (20)    A_j = a_j * Upsilon_j^lambda, Upsilon_j = sum_s e^{-delta tau_js} (H_Ms/K_s)
#     (21)    B_i = b_i * Omega_i^eta,      Omega_i    = sum_s e^{-rho tau_is}   (H_Rs/K_s)
#
# =====================================================================================
# !!! NOTATION COLLISIONS WITH THE REST OF THIS REPOSITORY -- READ BEFORE REUSING !!!
# =====================================================================================
#  lambda  ARSW: AGGLOMERATION elasticity of productivity, 0.0710. HIGHER lambda RAISES
#          productivity via spillovers.
#          AllenArkolakis-RES-2022: traffic CONGESTION elasticity, 0.07-0.09. HIGHER
#          lambda RAISES transport cost. Same letter, near-identical magnitude, OPPOSITE
#          economic content. A value transplanted between the two looks plausible and is
#          wrong. NEVER copy a calibrated lambda across these two models.
#  epsilon ARSW: Frechet shape for bilateral COMMUTING shocks, ~6.83.
#          This repo elsewhere: migration elasticity, 3-5.
#  alpha   ARSW: LABOUR share in production, 0.80 (floor space in firm costs = 1-alpha
#          = 0.20). NOT the floor-space share -- see the note at ARSWParams.
#          AllenArkolakis-RES-2022: productivity externality, 0.1.
#  beta    ARSW: expenditure share on the tradable numeraire, 0.75 (residential floor
#          space = 1-beta = 0.25).
#          AllenArkolakis-RES-2022: amenity externality, -0.3.
#  theta   ARSW: COMMERCIAL FLOOR-SPACE SHARE of a block, theta_i in [0,1].
#          Redding-JIE-2016 / AA2022: trade / Frechet elasticity, 6-10.
#  pi_ij   ARSW: COMMUTING probability (live i, work j).
#          Redding-JIE-2016: bilateral TRADE share. The code looks near-identical; the
#          economics is not.
#
# =====================================================================================
# RELATION TO THE OTHER MODELS IN THIS REPOSITORY
# =====================================================================================
#  * QSE-ARE-2017 (Helpman_E.jl) is the same family (CES consumption, Cobb-Douglas floor
#    space, agglomeration/amenity externalities) but WITHOUT bilateral commuting choice:
#    one location = one workplace and residence, no H_M != H_R, no pi_ij, no eps.
#  * Redding-JIE-2016 is the closest CODE analogue: thetaepsopt.jl <-> optimepsilon_TD86.m,
#    solveab.jl/solveHab.jl <-> cprod.m/cres.m. But its pi_ij is a goods-TRADE share.
#  * AllenArkolakis-RES-2022 endogenises transport COST (route choice + congestion);
#    here tau_ij is exogenous data. Conversely AA2022 has no floor-space market and no
#    separate workplace/residence margin. The two are complements, not nests.
#  * This model supplies the two margins the rest of the repository lacks: bilateral
#    commuting (residence x workplace Frechet) and an explicit floor-space market with
#    an endogenous commercial/residential land-use split.
#
# =====================================================================================
# UPSTREAM DEFECTS HANDLED HERE (each flagged again at the site)
# =====================================================================================
#  D1  cftualprep_end_TD.m:94-96 -- the inline comments for delta and eta are SWAPPED.
#      The code is right; only the comments lie. Not propagated.
#  D2  META.m hard-codes the repo root behind a `user` switch. This file is self-locating
#      via @__DIR__.
#  D3  mu is never a named variable upstream: it is the literal exponent 0.75 in four
#      files. Here it is a named parameter and the exponent is written (1 - mu).
#  D4  comegaoptC.m:206-214 ("New code by GA") fills cprob by indexing with the LOGICAL
#      VALUES of the mask; a 0 subscript is a MATLAB error. The original, commented-out,
#      correct line is used instead.
#  D5  smodexog.m:152 uses the FIXED observed QT inside the commuting probability where
#      smodendog.m:167 correctly uses the guess Q_i. Fixed here; `legacy_QT=true`
#      restores the upstream behaviour so the difference can be measured.
#  D6  smodendog.m:110 `toc(xtic./60)` divides the timer handle. Cosmetic; not ported.
#
# Dependencies (all already installed in this repo's global environment):
#   LinearAlgebra, Statistics, Printf, Random, SpecialFunctions, Optim, Plots
#
# Run with:  julia ahlfeldt_redding_sturm_wolf_density_model.jl
# =====================================================================================

using LinearAlgebra, Statistics, Printf, Random, SpecialFunctions, Optim

const MAKE_PLOTS = get(ENV, "ARSW_PLOTS", "1") != "0"
if MAKE_PLOTS
    ENV["GKSwstype"] = "100"          # headless GR; no window is opened
    using Plots
end

# =====================================================================================
# 1. PARAMETERS AND DATA STRUCTURES
# =====================================================================================

"""
    ARSWParams

Structural parameters. Greek names in comments, ASCII in code (repo convention).

`alpha` is the **labour** share in production: `Y = A * H_M^alpha * (theta*L)^(1-alpha)`,
so the floor-space share of firm costs is `1 - alpha = 0.20` (paper p. 2167:
"the share of firm expenditure on commercial floor space (1-alpha) is 0.20"). The
upstream toolkit's CLAUDE.md calls alpha the "floor-space input share", which is the
complement of the truth; `cftualprep_end_TD.m` itself comments
`alpha=0.80; % Set input share of labour in production`, which is correct, and
`comegaoptC.m:222`'s `A = (q/(1-alpha))^(1-alpha) (w/alpha)^alpha` only inverts Eq. (12)
under the labour-share reading.

`mu` (defect D3): upstream hard-codes the exponent on land as the literal `0.75`, which
reads as `1-mu = 0.75`, i.e. mu = 0.25. The paper (p. 2167) sets "the share of land in
construction costs (1-mu) equal to 0.25", i.e. mu = 0.75 and the exponent is 0.25. Here
mu is a parameter and the exponent is written `(1 - mu)`; `mu = 0.25` reproduces
upstream's phi exactly. The choice is numerically inert everywhere else (check C8).
"""
Base.@kwdef struct ARSWParams
    alpha    :: Float64 = 0.80    # alpha: labour share in production, Eq. (10)
    beta     :: Float64 = 0.75    # beta: expenditure share on the tradable good, Eq. (1)
    mu       :: Float64 = 0.75    # mu: non-land input share in floor-space supply, Eq. (15)
    epsilon  :: Float64 = 6.83    # epsilon: Frechet shape, commuting
    kappaeps :: Float64 = 0.07    # nu = kappa*epsilon: commuting decay per minute
    lambda   :: Float64 = 0.0     # lambda: density elasticity of PRODUCTIVITY, Eq. (20)
    delta    :: Float64 = 0.0     # delta: PRODUCTIVITY-spillover decay, Eq. (20)   [D1]
    eta      :: Float64 = 0.0     # eta: density elasticity of AMENITY, Eq. (21)    [D1]
    rho      :: Float64 = 0.0     # rho: AMENITY-spillover decay, Eq. (21)
end

"""
Damping weight used wherever the solver has to cross a LARGE shock (the wall going up or
down, a cold start from uniform guesses).

Upstream's `smodexog.m:286` fixes the weight on the predicted value at 0.5, and
`smodendog.m:293-304` uses 0.5 unless the max log gap is already below 0.1. On this
synthetic geography that fixed 0.5 does not converge once the wall is raised: the
iteration settles into a period-2 limit cycle with a max log gap stuck at ~0.63 for
20,000 iterations (the West/East population split is a very elastic margin at
eps = 6.83, so the map's Jacobian picks up an eigenvalue below -1 and the 50/50 blend
cannot damp it). A weight of 0.25 converges to 1e-13 in ~375 iterations. Upstream's own
comment anticipates this -- "it can be beneficial to choose a smaller weight if the
algorithm is bouncing too much" -- it just never lowers the weight while the gap is still
large. Measured explicitly in `run_counterfactuals`.
"""
const DAMP_ROBUST = 0.25

kappa_of(p::ARSWParams)  = p.kappaeps / p.epsilon          # kappa = nu / epsilon
gammaf_of(p::ARSWParams) = gamma((p.epsilon - 1) / p.epsilon)   # gamma in Eq. (9)

"GMM parameter vector of the paper's Section 7 (Table V), loaded upstream from roptimis_all_big.mat."
const PARAMS_SECTION7 = ARSWParams(alpha = 0.80, beta = 0.75, mu = 0.75,
                                   epsilon  = 6.694,   # EThetaA(2)
                                   kappaeps = 0.0987,  # EThetaA(1)  nu
                                   lambda   = 0.0710,  # EThetaA(3)  density elast. of PRODUCTIVITY
                                   delta    = 0.3617,  # EThetaA(4)  PRODUCTIVITY decay      [D1]
                                   eta      = 0.1553,  # EThetaA(5)  density elast. of AMENITY [D1]
                                   rho      = 0.7595)  # EThetaA(6)  AMENITY decay

"Section 6 parameters: no agglomeration spillovers, reduced-form nu, one-step epsilon."
const PARAMS_SECTION6 = ARSWParams()

"""
    ARSWData

The observed data the model is inverted from. Field names map onto the MATLAB variables
built in `prepdata_TD.m` and assembled into `obsvar06` by every driver script:

| field    | MATLAB      | obsvar06 col | units                                |
|----------|-------------|--------------|--------------------------------------|
| `Q`      | `floor06`   | 1            | floor-space price                     |
| `HM`     | `empwpl06`  | 2            | workplace employment                  |
| `HR`     | `emprsd06`  | 3            | residence employment                  |
| `K`      | `area06`    | 4            | geographic area                       |
| `tau`    | `tt06`      | --           | N x N bilateral travel time (minutes) |
| `west`   | `dummywestr`| --           | former-West-Berlin flag               |
| `coords` | `Xr`,`Yr`   | --           | block centroids, plotting only        |
"""
struct ARSWData
    N      :: Int
    Q      :: Vector{Float64}
    HM     :: Vector{Float64}
    HR     :: Vector{Float64}
    K      :: Vector{Float64}
    tau    :: Matrix{Float64}
    west   :: Vector{Bool}
    coords :: Matrix{Float64}
end

"""
    Fund

Bundle handed to the equilibrium solvers. Mirrors the `fund` matrix that
`smodexog.m` / `smodendog.m` / `ussmodendog.m` take, but with names instead of column
numbers (upstream's columns differ between the two solvers, which is a standing trap).

* For `solve_equilibrium_exog` (Alg. 11), `A`/`B` are TOTAL adjusted productivity and
  amenity `{A~, B~}` -- there are no spillovers.
* For `solve_equilibrium_endog` (Alg. 15/16), `A`/`B` are the FUNDAMENTAL components
  `{a~, b~}` and the total `{A,B}` is rebuilt each iteration from Eqs. (20)-(21).

`LD` is total floor space SUPPLY (exogenous, `= V .* K.^(1-mu)`); `theta`, `wage`, `vv`,
`Q`, `HM`, `HR` are starting guesses for the target variables.
"""
Base.@kwdef struct Fund
    A     :: Vector{Float64}
    B     :: Vector{Float64}
    V     :: Vector{Float64}       # phi, density of development
    K     :: Vector{Float64}
    Q     :: Vector{Float64}
    HM    :: Vector{Float64}
    HR    :: Vector{Float64}
    LD    :: Vector{Float64}
    theta :: Vector{Float64}
    wage  :: Vector{Float64}
    vv    :: Vector{Float64}
end

"Solved equilibrium. Column order of upstream's `endog` matrix is kept in `as_matrix`."
Base.@kwdef struct Equilibrium
    wage  :: Vector{Float64}
    vv    :: Vector{Float64}
    theta :: Vector{Float64}
    Y     :: Vector{Float64}
    Q     :: Vector{Float64}
    q     :: Vector{Float64}
    HM    :: Vector{Float64}
    HR    :: Vector{Float64}
    rent  :: Vector{Float64}       # Crent: q where commercial, Q where residential
    A     :: Vector{Float64}       # total productivity (= a * Upsilon^lambda)
    B     :: Vector{Float64}       # total amenity      (= b * Omega^eta)
    ucprob:: Matrix{Float64}       # pi_ij, rows = residences, cols = workplaces
    Phi   :: Float64
    Ubar  :: Float64
    HH    :: Float64
    converged :: Bool
    iters :: Int
    path  :: Matrix{Float64}       # convergence path: [maxLDwage maxLDq maxLDQ maxLDtheta maxLD iter]
end

# --- small helpers -------------------------------------------------------------------

"Geometric mean over strictly positive entries. MATLAB's `geomean` is always called on
already-masked subvectors in this toolkit, so this matches its behaviour."
function geomean_pos(x::AbstractVector)
    v = filter(>(0.0), x)
    isempty(v) && return 1.0
    return exp(mean(log.(v)))
end

"max_i |log(a_i / b_i)| over entries where both are strictly positive. Scale-free."
function maxlogdiff(a::AbstractVector, b::AbstractVector)
    m = 0.0
    @inbounds for i in eachindex(a)
        if a[i] > 0 && b[i] > 0
            d = abs(log(a[i] / b[i]))
            d > m && (m = d)
        end
    end
    return m
end
maxlogdiff(a::AbstractVector, b::AbstractVector, mask::AbstractVector{Bool}) =
    maxlogdiff(a[mask], b[mask])

"Upstream's rounding stopping rule: round(x*scale) elementwise and require exact equality."
rounded_equal(a, b; scale = 100.0) = round.(a .* scale) == round.(b .* scale)

# =====================================================================================
# 2. SHARED COMMUTING KERNEL
# =====================================================================================
# Every algorithm below builds the same object. Upstream re-derives it in eight separate
# files (comegaoptC, camen, expincome, cmodexog, smodexog, smodendog, ussmodendog, ubar),
# which is where most of the copy-paste risk in the MATLAB lives.
#
# CONVENTION, fixed once and never varied: phi_ij has RESIDENCES in rows, WORKPLACES in
# columns. Upstream flips this between files (comegaopt* has workplaces in rows); every
# transpose below is deliberate.

"""
    commuting_phi(p, tau_sub, B_from, Q_from, w_to)

Numerator of Eq. (4) on a residence-subset x workplace-subset rectangle:

    phi_ij = exp(-eps*kappa*tau_ij) * B_i^eps * Q_i^(-(1-beta)*eps) * w_j^eps

`tau_sub` is `tau[Ifrom, Ito]` (nfrom x nto). Returns nfrom x nto.
"""
function commuting_phi(p::ARSWParams, tau_sub::AbstractMatrix,
                       B_from::AbstractVector, Q_from::AbstractVector, w_to::AbstractVector)
    eps = p.epsilon
    d_ij_eps = exp.(-p.kappaeps .* tau_sub)                       # nu = eps*kappa
    return d_ij_eps .* (B_from .^ eps) .* (Q_from .^ (-(1 - p.beta) * eps)) .*
           transpose(w_to .^ eps)
end

"""
    commuting_block(p, phi_ij, HH)

Everything Eqs. (5), (6), (9) and (S.20) need, from the phi_ij rectangle.
Returns a NamedTuple with

* `Phi`      scalar sum over all bilaterals (denominator of Eq. 4)
* `phi_i`, `phi_j` row / column sums (numerators of Eq. 5)
* `pp_ij`    unconditional commuting probabilities pi_ij
* `pp_iji`   pi_{j|i}: probability of working at j GIVEN residence i (rows sum to 1)
* `pp_ijj`   pi_{i|j}: probability of living at i GIVEN workplace j (nto x nfrom, rows sum to 1)
* `HR`,`HM`  implied residence / workplace employment
* `Ubar`     reservation utility, Eq. (9)
"""
function commuting_block(p::ARSWParams, phi_ij::AbstractMatrix, HH::Real)
    Phi   = sum(phi_ij)                                # MATLAB `sum(sum(.))`; Julia `sum` is already whole-array
    phi_i = vec(sum(phi_ij, dims = 2))                 # over workplaces  -> per residence
    phi_j = vec(sum(phi_ij, dims = 1))                 # over residences  -> per workplace
    pp_ij  = phi_ij ./ Phi
    pp_i   = phi_i ./ Phi
    pp_j   = phi_j ./ Phi
    pp_iji = phi_ij ./ phi_i                           # broadcast column vector down the rows
    pp_ijj = transpose(phi_ij) ./ phi_j
    return (Phi = Phi, phi_i = phi_i, phi_j = phi_j,
            pp_ij = pp_ij, pp_i = pp_i, pp_j = pp_j,
            pp_iji = pp_iji, pp_ijj = pp_ijj,
            HR = pp_i .* HH, HM = pp_j .* HH,
            Ubar = gammaf_of(p) * Phi^(1 / p.epsilon))
end

"""
    commuting_market_access(p, tau_sub, w_to)

Residential commuting market access, Eq. (S.46) / paper Eq. (29):
`W_i = sum_s (w~_s / e^{kappa tau_is})^eps = sum_s e^{-nu tau_is} w~_s^eps`.
`camen.m:68` calls this `ECMA`. `tau_sub` is `tau[Ifrom, Ito]`.
"""
commuting_market_access(p::ARSWParams, tau_sub::AbstractMatrix, w_to::AbstractVector) =
    exp.(-p.kappaeps .* tau_sub) * (w_to .^ p.epsilon)

# =====================================================================================
# 3. ALGORITHMS 1 AND 4 -- TRANSFORMED WAGES, ADJUSTED WAGES, PRODUCTIVITY
#    comegaoptO.m (Alg. 1)  and  comegaoptC.m (Alg. 4)
# =====================================================================================

"""
    solve_omega_fixed_point(p, HM, HR, tau; ...)

Core of BOTH Algorithm 1 (`comegaoptO.m`) and Algorithm 4 (`comegaoptC.m`): the
commuting-market-clearing fixed point in transformed wages, Eq. (S.44)

    H_Mj = sum_i [ omega_j e^{-nu tau_ij} / sum_s omega_s e^{-nu tau_is} ] H_Ri

The conditional commuting probability pi_{j|i} does not contain B_i or Q_i (they cancel
out of the conditional), so omega is identified from {H_M, H_R, tau} ALONE -- and from nu
alone, not epsilon. That is the observation the paper's one-step estimation of epsilon
rests on (p. 2167): "transformed wages are determined independently of epsilon".

Lemmas S.6-S.7 establish that this system satisfies gross substitution, so the fixed
point is unique up to scale; we normalise `geomean(omega) = 1`.

Orientation here follows upstream: `cprob` has WORKPLACES in rows and RESIDENCES in
columns, and each COLUMN sums to one.

`rule = :tol` (default) stops on `max_j |HMC_j/HMT_j - 1| < tol`; `rule = :round`
reproduces upstream's `round(x*round_scale)` equality test exactly (`1e5` in
`comegaoptO.m`, `1e4` in `comegaoptC.m`).
"""
function solve_omega_fixed_point(p::ARSWParams, HM::Vector{Float64}, HR::Vector{Float64},
                                 tau::Matrix{Float64};
                                 omega0::Union{Nothing,Vector{Float64}} = nothing,
                                 maxiter::Int = 20_000, tol::Float64 = 1e-14,
                                 rule::Symbol = :tol, round_scale::Float64 = 1e4,
                                 damp::Float64 = 0.5, verbose::Bool = false,
                                 rng::AbstractRNG = Random.default_rng())
    N     = length(HM)
    Iwpl  = HM .!= 0                  # blocks with positive workplace employment
    Irsd  = HR .!= 0                  # blocks with positive residence employment
    nto   = count(Iwpl)
    nfrom = count(Irsd)
    EHMT  = HM[Iwpl]
    EHRT  = HR[Irsd]

    # nu*tau on the (workplace x residence) rectangle. Upstream builds Ecc = exp(+nu*tau)
    # and divides; multiplying by exp(-nu*tau) is the same and avoids an overflow when a
    # wall penalty makes nu*tau large.
    Ecc = exp.(-p.kappaeps .* tau[Iwpl, Irsd])                      # nto x nfrom

    Eomega = omega0 === nothing ? ones(nto) : omega0[Iwpl]
    Eomega = Eomega ./ geomean_pos(Eomega)

    # Pre-declared: assigned inside the loop, needed afterwards (Julia loop scoping trap).
    Ecprob    = zeros(nto, nfrom)
    EHMC      = zeros(nto)
    gap       = Inf
    converged = false
    iters     = 0

    for x in 1:maxiter
        iters = x
        Ecnum   = Eomega .* Ecc                                      # numerator of Eq. S.44
        Ecdenom = sum(Ecnum, dims = 1)                               # MATLAB `sum(M)` = COLUMN sums
        Ecprob  = Ecnum ./ Ecdenom                                   # pi_{j|i}; columns sum to 1
        EHMC    = Ecprob * EHRT                                      # sum over i in Eq. S.44

        if rule === :round
            converged = rounded_equal(EHMC, EHMT; scale = round_scale)
            gap = maximum(abs.(round.(EHMC .* round_scale) .- round.(EHMT .* round_scale)))
        else
            gap = maximum(abs.(EHMC ./ EHMT .- 1))
            converged = gap < tol
        end
        if converged
            verbose && println(">>>> Wage System Converged <<<< (iter $x, gap $gap)")
            break
        end

        Eomega_e = (EHMT ./ EHMC) .* Eomega       # inflate where we under-predict H_M
        if isnan(gap)                             # upstream's random-restart safety net
            Eomega_e = 0.95 .+ 0.10 .* rand(rng, nto)
            Eomega   = 0.95 .+ 0.10 .* rand(rng, nto)
        end
        Eomega = damp .* Eomega_e .+ (1 - damp) .* Eomega
        Eomega = Eomega ./ geomean_pos(Eomega)    # omega identified up to scale
    end

    omega = zeros(N);  omega[Iwpl] = Eomega       # zero omega rationalises zero employment
    HMC   = zeros(N);  HMC[Iwpl]   = EHMC
    # [D4] Upstream's "New code by GA" sparse fill is broken (indexes with the logical
    # VALUES of the mask). The original, correct, commented-out line is used instead.
    cprob = zeros(N, N)
    cprob[Iwpl, Irsd] = Ecprob

    return (omega = omega, cprob = cprob, converged = converged, HMC = HMC,
            gap = gap, iters = iters)
end

"""
    solve_transformed_wages(p, data; kwargs...)  --  ALGORITHM 1, `comegaoptO.m`

Thin wrapper: returns transformed wages `omega_j = w~_j^eps` given `{H_M, H_R, tau}`.
Used by the epsilon estimation, where the point is precisely that omega does not depend
on epsilon.
"""
solve_transformed_wages(p::ARSWParams, data::ARSWData; kwargs...) =
    solve_omega_fixed_point(p, data.HM, data.HR, data.tau; round_scale = 1e5, kwargs...)

"""
    solve_wages_and_productivity(p, data; kwargs...)  --  ALGORITHM 4, `comegaoptC.m`

Solves Eq. (S.44) for transformed wages, maps back to adjusted wages
`w~_j = omega_j^(1/eps)` normalised to `geomean(w~) = 1`, then recovers adjusted
productivity from the zero-profit condition Eq. (S.48) inverted:

    A~_j = (q_j / (1-alpha))^(1-alpha) * (w~_j / alpha)^alpha

with `q_j` taken from the observed floor price. Blocks with `H_M = 0` get
theory-consistent zeros (Lemmas S.1-S.5).
"""
function solve_wages_and_productivity(p::ARSWParams, data::ARSWData; kwargs...)
    N    = data.N
    Iwpl = data.HM .!= 0
    res  = solve_omega_fixed_point(p, data.HM, data.HR, data.tau; round_scale = 1e4, kwargs...)

    wage = res.omega .^ (1 / p.epsilon)                               # omega -> w~
    pos  = wage .> 0
    wage[pos] ./= geomean_pos(wage[pos])                              # normalisation of Sec. S.3.1.5

    A = zeros(N)
    A[Iwpl] = ((data.Q[Iwpl] ./ (1 - p.alpha)) .^ (1 - p.alpha)) .*
              ((wage[Iwpl] ./ p.alpha) .^ p.alpha)                    # Eq. (S.48) inverted
    return (wage = wage, A = A, cprob = res.cprob, converged = res.converged,
            HMC = res.HMC, gap = res.gap, iters = res.iters)
end

# =====================================================================================
# 4. ALGORITHM 5 -- ADJUSTED AMENITIES FROM COMMUTING MARKET ACCESS
#    camen.m
# =====================================================================================

"""
    solve_amenities(p, data, wage)  --  ALGORITHM 5, `camen.m`

Inverts Eq. (S.47),

    B~_i / Bbar = (H_Ri/Hbar_R)^(1/eps) * (Q_i/Qbar)^(1-beta) * (W_i/Wbar)^(-1/eps)

where `W_i` is residential commuting market access, Eq. (S.46). Intuitively (supplement
p. 30): high residence employment and high floor prices must be explained either by good
commuting market access or by attractive residential amenities.

Every factor is divided by its geometric mean first, so the returned `B` automatically
satisfies `geomean(B) = 1`; Algorithm 6 then puts it on the level that makes `Phi = H`.
Returns `B`, the UN-normalised `CMA`, and residence-employment shares `HRS`.
"""
function solve_amenities(p::ARSWParams, data::ARSWData, wage::Vector{Float64})
    N    = data.N
    Iwpl = data.HM .!= 0
    Irsd = data.HR .!= 0

    tau_sub = data.tau[Irsd, Iwpl]                                   # nfrom x nto
    EHRT    = data.HR[Irsd]

    EHRS = EHRT ./ sum(EHRT)                                         # residence shares
    EHRSn = EHRS ./ geomean_pos(EHRS)
    EQTn  = data.Q[Irsd] ./ geomean_pos(data.Q[Irsd])
    ECMA  = commuting_market_access(p, tau_sub, wage[Iwpl])          # Eq. (S.46)
    ECMAn = ECMA ./ geomean_pos(ECMA)

    EB = (EHRSn .^ (1 / p.epsilon)) .* (EQTn .^ (1 - p.beta)) ./ (ECMAn .^ (1 / p.epsilon))

    B   = zeros(N);  B[Irsd]   = EB                                   # zero H_R -> zero B~
    CMA = zeros(N);  CMA[Irsd] = ECMA
    HRS = zeros(N);  HRS[Irsd] = EHRS
    return (B = B, CMA = CMA, HRS = HRS)
end

# =====================================================================================
# 5. ALGORITHM 6 -- RESCALING {A~, B~} ONTO THE SOLVER'S NORMALISATION
#    calcal_adj_TD.m  (NOT part of the original replication directory)
# =====================================================================================

"""
    rescale_fundamentals(p, data, A, B)  --  ALGORITHM 6, `calcal_adj_TD.m`

The sequential procedure (Algs. 4-5) returns `{A~, B~}` on normalisations that are fine
for cross-location comparisons but not for seeding an equilibrium solver. This puts them
on the same footing as the simultaneous procedure (Alg. 9):

1. `A~ <- A~ / geomean(A~)`, since Alg. 4 normalised *wages* to geometric mean one, which
   does not imply that productivities have geometric mean one;
2. recompute wages consistent with the rescaled `A~` from Eq. (12) solved for `w`,
   `w = alpha ((1-alpha)/q)^((1-alpha)/alpha) A^(1/alpha)`;
3. `B~ <- (H/Phi)^(1/eps) * B~`, so that model city population `Phi` matches observed
   total employment `H`. Because `B` enters Eq. (4) only as `B_i^eps`, a common factor
   pulls straight out of the sum: `H = c^eps * Phi` gives `c = (H/Phi)^(1/eps)`.

After this step the sequential and simultaneous procedures agree numerically -- which is
exactly what test R4 checks.
"""
function rescale_fundamentals(p::ARSWParams, data::ARSWData,
                              A::Vector{Float64}, B::Vector{Float64})
    N    = data.N
    Iwpl = data.HM .!= 0
    Irsd = data.HR .!= 0
    HH   = sum(data.HM[Iwpl])

    EA = A[Iwpl] ./ geomean_pos(A[Iwpl])
    # Eq. (12) solved for w: combines the FOC for labour with zero profits.
    Ewage = (((1 - p.alpha) ./ data.Q[Iwpl]) .^ ((1 - p.alpha) / p.alpha)) .*
            p.alpha .* (EA .^ (1 / p.alpha))

    phi_ij = commuting_phi(p, data.tau[Irsd, Iwpl], B[Irsd], data.Q[Irsd], Ewage)
    Phi    = sum(phi_ij)
    EB     = ((HH / Phi) ^ (1 / p.epsilon)) .* B[Irsd]

    Aout = zeros(N);  Aout[Iwpl] = EA
    Bout = zeros(N);  Bout[Irsd] = EB
    wout = zeros(N);  wout[Iwpl] = Ewage
    return (A = Aout, B = Bout, wage = wout)
end

# =====================================================================================
# 6. ALGORITHM 7 -- TOTAL EXPECTED WORKER INCOME
#    expincome.m
# =====================================================================================

"""
    expected_income(p, data, wage, B)  --  ALGORITHM 7, `expincome.m`

Eq. (S.20): `E[w|i] = sum_j pi_{j|i} w~_j`, then total worker income
`vv_i = E[w|i] * H_Ri`. This is the object that enters residential land market clearing,
Eq. (19)/(S.29).

Upstream replaces the vectorised `repmat` construction with an explicit double `for`
loop purely to save memory at 15,000 blocks; the loop and the vectorised version compute
the same thing, so only the vectorised form is kept here.
"""
function expected_income(p::ARSWParams, data::ARSWData, wage::Vector{Float64},
                         B::Vector{Float64})
    N    = data.N
    Iwpl = data.HM .!= 0
    Irsd = data.HR .!= 0

    phi_ij = commuting_phi(p, data.tau[Irsd, Iwpl], B[Irsd], data.Q[Irsd], wage[Iwpl])
    cb     = commuting_block(p, phi_ij, sum(data.HM[Iwpl]))

    EEWI = cb.pp_iji * wage[Iwpl]                    # E[w|i], Eq. (S.20)
    ETWI = EEWI .* data.HR[Irsd]                     # total worker income

    vv = zeros(N);  vv[Irsd] = ETWI
    Ew = zeros(N);  Ew[Irsd] = EEWI
    return (vv = vv, Ew = Ew, Phi = cb.Phi, HMC = cb.HM, HRC = cb.HR)
end

# =====================================================================================
# 7. ALGORITHMS 8 AND 10 -- DENSITY OF DEVELOPMENT, FLOOR SPACE, LAND-USE SHARE
#    cdensity.m (Alg. 8)  and  cdensityE.m (Alg. 10, identical but also returns LM, LR)
# =====================================================================================

"""
    solve_density(p, data, A, vv)  --  ALGORITHMS 8 and 10, `cdensity.m` / `cdensityE.m`

Applies land market clearing in reverse: given the observed floor price and the inverted
`{A~, vv}`, read off the floor space each use demands, add them up, and the total must be
the supply.

    (S.30)/(18)  LM_i = theta_i L_i     = ((1-alpha) A~_i / q_i)^(1/alpha) * H_Mi
    (S.29)/(19)  LR_i = (1-theta_i) L_i = (1-beta) vv_i / Q_i
    (S.31)/(15)  L_i  = LM_i + LR_i = phi_i * K_i^(1-mu)   ->   phi_i = L_i / K_i^(1-mu)
                 theta_i = LM_i / L_i            (COMMERCIAL floor-space share)

Note both lines use the SINGLE observed floor price. That is coherent because a block
with `H_M = 0` contributes no `LM` and a block with `H_R = 0` no `LR`, while an
incompletely specialised block has `q_i = Q_i` exactly (supplement S.2.6). Asserted in
the validation block rather than left implicit.

[D3] The exponent is `(1 - mu)`, not the hard-coded literal `0.75` of `cdensity.m:57`.
"""
function solve_density(p::ARSWParams, data::ARSWData, A::Vector{Float64}, vv::Vector{Float64})
    N    = data.N
    Iwpl = data.HM .!= 0
    Irsd = data.HR .!= 0

    LM = zeros(N)
    LM[Iwpl] = ((((1 - p.alpha) .* A[Iwpl]) ./ data.Q[Iwpl]) .^ (1 / p.alpha)) .* data.HM[Iwpl]

    LR = zeros(N)
    LR[Irsd] = ((1 - p.beta) .* vv[Irsd]) ./ data.Q[Irsd]

    LD    = LM .+ LR
    V     = LD ./ (data.K .^ (1 - p.mu))              # phi, Eq. (S.31)
    theta = LM ./ LD                                  # commercial share
    return (V = V, LD = LD, LM = LM, LR = LR, theta = theta)
end

# =====================================================================================
# 8. ALGORITHM 9 -- SIMULTANEOUS INVERSION OF {A~, B~}
#    cmodexog.m
# =====================================================================================

"""
    invert_simultaneous(p, data; ...)  --  ALGORITHM 9, `cmodexog.m`

One outer fixed point that jointly updates guesses of `{A~, B~}` until the model's
PREDICTED `{H_Mj, H_Ri}` (Eqs. 4-5 given the guesses and the observed floor prices) match
the OBSERVED values. This, not the sequential route, is the quantification the paper
actually uses to seed counterfactuals (`cftualexog_TD.m:45`).

Each pass:
1. wages from Eq. (12) given the productivity guess and the observed floor price;
2. `phi_ij`, `Phi`, and the implied `{H_M, H_R}` from Eqs. (4)-(5);
3. inflate `A~` by `(H_M^obs/H_M^pred)^(1/eps)` and `B~` by `(H_R^obs/H_R^pred)^(1/eps)`;
4. blend 50/50 with the old guess, renormalise `geomean(A~) = 1` and `Phi = H`.

The two renormalisations are what makes the map well posed: scaling all `A~` by a
constant leaves `{H_M, H_R, theta, pi_ij}` untouched, and scaling all `B~` leaves the
entire allocation untouched (both verified numerically in check C7).

Upstream caps this at 200 iterations with a `round(gap*1e4) <= 0` rule; here `rule=:tol`
drives it to machine precision, which is what the recovery test needs.
"""
function invert_simultaneous(p::ARSWParams, data::ARSWData;
                             A0::Union{Nothing,Vector{Float64}} = nothing,
                             B0::Union{Nothing,Vector{Float64}} = nothing,
                             maxiter::Int = 20_000, tol::Float64 = 1e-14,
                             rule::Symbol = :tol, round_scale::Float64 = 1e4,
                             damp::Float64 = 0.5, verbose::Bool = false,
                             rng::AbstractRNG = Random.default_rng())
    N    = data.N
    Iwpl = data.HM .!= 0
    Irsd = data.HR .!= 0
    EHMT = data.HM[Iwpl]
    EHRT = data.HR[Irsd]
    HH   = sum(EHMT)
    EQw  = data.Q[Iwpl]
    EQr  = data.Q[Irsd]
    tau_sub = data.tau[Irsd, Iwpl]

    EA = A0 === nothing ? ones(count(Iwpl)) : A0[Iwpl]
    EB = B0 === nothing ? ones(count(Irsd))  : B0[Irsd]

    # Pre-declared for use after the loop.
    Ewage = zeros(length(EA));  cb = nothing
    mAgap = Inf;  mBgap = Inf;  converged = false;  iters = 0

    for x in 1:maxiter
        iters = x
        # Eq. (12) solved for w: w = alpha ((1-alpha)/q)^((1-alpha)/alpha) A^(1/alpha)
        Ewage = (((1 - p.alpha) ./ EQw) .^ ((1 - p.alpha) / p.alpha)) .*
                p.alpha .* (EA .^ (1 / p.alpha))

        phi_ij = commuting_phi(p, tau_sub, EB, EQr, Ewage)
        cb     = commuting_block(p, phi_ij, HH)

        if rule === :round
            mAgap = maximum(abs.(round.((EHMT .- cb.HM) .* round_scale)))
            mBgap = maximum(abs.(round.((EHRT .- cb.HR) .* round_scale)))
            converged = (mAgap <= 0) && (mBgap <= 0)
        else
            mAgap = maximum(abs.(cb.HM ./ EHMT .- 1))
            mBgap = maximum(abs.(cb.HR ./ EHRT .- 1))
            converged = (mAgap < tol) && (mBgap < tol)
        end
        if converged
            verbose && println(">>>> Calibration Convergence Achieved <<<< (iter $x)")
            break
        end

        EA_e = ((EHMT ./ cb.HM) .^ (1 / p.epsilon)) .* EA
        EB_e = ((EHRT ./ cb.HR) .^ (1 / p.epsilon)) .* EB
        if isnan(mBgap) || isnan(mAgap)                     # upstream's safety net
            EA_e = 0.95 .+ 0.10 .* rand(rng, length(EA))
            EB_e = 0.95 .+ 0.10 .* rand(rng, length(EB))
            EA, EB = EA_e, EB_e
        end
        EA = damp .* EA_e .+ (1 - damp) .* EA
        EB = damp .* EB_e .+ (1 - damp) .* EB
        EA = EA ./ geomean_pos(EA)                          # geomean(A~) = 1
        EB = ((HH / cb.Phi) ^ (1 / p.epsilon)) .* EB        # Phi = H
    end

    A    = zeros(N);  A[Iwpl]    = EA
    B    = zeros(N);  B[Irsd]    = EB
    wage = zeros(N);  wage[Iwpl] = Ewage
    HMC  = zeros(N);  HMC[Iwpl]  = cb.HM
    HRC  = zeros(N);  HRC[Irsd]  = cb.HR

    EEWI = cb.pp_iji * Ewage                                # Eq. (S.20)
    vv   = zeros(N);  vv[Irsd] = EEWI .* EHRT

    CMA  = zeros(N);  CMA[Irsd] = commuting_market_access(p, tau_sub, Ewage)
    ucprob = zeros(N, N);  ucprob[Irsd, Iwpl] = cb.pp_ij

    return (A = A, B = B, wage = wage, ucprob = ucprob, vv = vv, HMC = HMC, HRC = HRC,
            CMA = CMA, Phi = cb.Phi, HH = HH, converged = converged,
            mAgap = mAgap, mBgap = mBgap, iters = iters)
end

# =====================================================================================
# 9. ALGORITHM 11 -- EQUILIBRIUM WITH EXOGENOUS FUNDAMENTALS (CLOSED CITY)
#    smodexog.m
# =====================================================================================

"""
    solve_equilibrium_exog(p, f, tau; ...)  --  ALGORITHM 11, `smodexog.m`

Closed-city equilibrium for given `{A~, B~, L}`: total employment `H` is fixed and `Ubar`
adjusts. Target variables are `{w, q, Q, theta}`; the map is

  guess {w,q,Q,theta}
    -> commuting probabilities and {H_M, H_R}, Eqs. (4)-(5)
    -> output Y, Eq. (10), and predicted wage from the labour FOC, w = alpha Y / H_M
    -> total worker income vv, Eq. (S.20)
    -> floor prices from land market clearing, Eqs. (18)-(19)
    -> land-use share theta from the floor-space input demand, Eq. (18)
  -> compare, damp, repeat.

Three land-use regimes are handled exactly as upstream, via the index sets
`IcsA` (A>0, B=0: completely specialised commercial), `IcsB` (B>0, A=0: completely
specialised residential) and `Iis` (both positive: incompletely specialised, where
q_i = Q_i and theta is interior).

[D5] `smodexog.m:152` builds the commuting term from the FIXED observed `QT` rather than
the current guess `Q_i`, although its own comment says "guesses of residence floor space
prices"; `smodendog.m:167` correctly uses `Q_i`. With `QT` the residential price feedback
into location choice is switched off, so counterfactual price changes do not affect
commuting and check C5 (lambda=delta=eta=rho=0 => smodendog == smodexog) cannot hold.
This port uses `Q_i`; pass `legacy_QT = true` to reproduce upstream and measure the gap.
"""
function solve_equilibrium_exog(p::ARSWParams, f::Fund, tau::Matrix{Float64};
                                maxiter::Int = 20_000, tol::Float64 = 1e-13,
                                rule::Symbol = :tol, damp::Float64 = 0.5,
                                legacy_QT::Bool = false, legacy_Qe_init::Bool = false,
                                verbose::Bool = false)
    N = length(f.A)
    Ito   = f.A .!= 0                     # positive productivity  -> possible workplace
    Ifrom = f.B .!= 0                     # positive amenity       -> possible residence
    IcsA  = Ito  .& .!Ifrom
    IcsB  = Ifrom .& .!Ito
    Iis   = Ito  .& Ifrom
    EA, EB0 = f.A[Ito], f.B[Ifrom]
    L  = f.LD
    HH = sum(f.HM)

    tau_sub = tau[Ifrom, Ito]

    Q_i = zeros(N);  Q_i[Ifrom] = f.Q[Ifrom]
    q_i = zeros(N);  q_i[Ito]   = f.Q[Ito]
    # [D9, benign] smodexog.m:116-117 initialises Q_e = q_e = QT (the FULL observed
    # vector) where smodendog.m:107-108 uses Q_e = Q_i, q_e = q_i. Inside the loop only
    # Q_e[IcsB], Q_e[Iis], q_e[IcsA] and q_e[Iis] are written, so upstream leaves
    # Q_e[IcsA] = QT[IcsA] against Q_i[IcsA] = 0, and its stopping rule compares the FULL
    # vectors -- which looks like it must fail forever on any completely specialised block.
    # It does not: the damping blend q_i <- w*q_e + (1-w)*q_i drags the structurally unused
    # entries onto their own initial values too, so the test is eventually satisfied. The
    # cost is only extra iterations (32 vs 1 on the :specialised scenario) plus meaningless
    # nonzero entries in the returned q[IcsB] / Q[IcsA] -- which `Crent` never reads.
    # Measured in run_model_checks (C9), not assumed: the first guess here was that it
    # broke convergence outright, and the measurement said otherwise.
    Q_e = legacy_Qe_init ? copy(f.Q) : copy(Q_i)
    q_e = legacy_Qe_init ? copy(f.Q) : copy(q_i)
    theta_i = copy(f.theta)
    Ewage_i = f.wage[Ito]
    QT_fixed = copy(f.Q)                  # only used when legacy_QT = true  [D5]

    Y = zeros(N); vv = copy(f.vv)
    wage_i = zeros(N); wage_e = zeros(N)
    theta_e = copy(theta_i)
    cb = nothing;  EHM = f.HM[Ito];  EHR = f.HR[Ifrom]
    converged = false;  iters = 0
    path = Matrix{Float64}(undef, 0, 6)
    pathv = Vector{NTuple{6,Float64}}()

    for x in 1:maxiter
        iters = x
        Qrow   = legacy_QT ? QT_fixed[Ifrom] : Q_i[Ifrom]      # [D5]
        phi_ij = commuting_phi(p, tau_sub, EB0, Qrow, Ewage_i)
        cb     = commuting_block(p, phi_ij, HH)
        EHR, EHM = cb.HR, cb.HM

        # Output, Eq. (10), at guessed land use and predicted employment
        Y[Ito] = EA .* (EHM .^ p.alpha) .* ((theta_i[Ito] .* L[Ito]) .^ (1 - p.alpha))
        # Predicted wage from the labour FOC of Eq. (10)
        Ewage_e = (p.alpha .* Y[Ito]) ./ EHM
        # Total worker income, Eq. (S.20)
        vv[Ifrom] = (cb.pp_iji * Ewage_i) .* EHR

        # Floor prices from land market clearing
        q_e[IcsA] = ((1 - p.alpha) .* Y[IcsA]) ./ (theta_i[IcsA] .* L[IcsA])
        Q_e[IcsB] = ((1 - p.beta)  .* vv[IcsB]) ./ ((1 .- theta_i[IcsB]) .* L[IcsB])
        mixed_rent = (((1 - p.alpha) .* Y[Iis]) .+ ((1 - p.beta) .* vv[Iis])) ./ L[Iis]
        q_e[Iis] = mixed_rent
        Q_e[Iis] = mixed_rent

        # Land-use share from commercial floor-space input demand, Eq. (18)
        theta_e = copy(theta_i)
        theta_e[Iis] = ((1 - p.alpha) .* Y[Iis]) ./ (q_e[Iis] .* L[Iis])

        wage_i[Ito] = Ewage_i
        wage_e[Ito] = Ewage_e

        if rule === :round
            converged = rounded_equal(wage_e, wage_i) && rounded_equal(q_e, q_i) &&
                        rounded_equal(Q_e, Q_i)       && rounded_equal(theta_e, theta_i)
            mLDw = maxlogdiff(wage_e, wage_i, Ito);  mLDq = maxlogdiff(q_e, q_i, Ito)
            mLDQ = maxlogdiff(Q_e, Q_i, Ifrom);      mLDt = maxlogdiff(theta_e, theta_i, Iis)
        else
            mLDw = maxlogdiff(wage_e, wage_i, Ito);  mLDq = maxlogdiff(q_e, q_i, Ito)
            mLDQ = maxlogdiff(Q_e, Q_i, Ifrom);      mLDt = maxlogdiff(theta_e, theta_i, Iis)
            converged = max(mLDw, mLDq, mLDQ, mLDt) < tol
        end
        mLD = max(mLDw, mLDq, mLDQ, mLDt)
        push!(pathv, (mLDw, mLDq, mLDQ, mLDt, mLD, Float64(x)))
        if converged
            verbose && println(">>>> Convergence Achieved <<<< (iter $x, maxLD $mLD)")
            break
        end

        Ewage_i = damp .* Ewage_e .+ (1 - damp) .* Ewage_i
        q_i     = damp .* q_e     .+ (1 - damp) .* q_i
        Q_i     = damp .* Q_e     .+ (1 - damp) .* Q_i
        theta_i = damp .* theta_e .+ (1 - damp) .* theta_i
    end

    path = isempty(pathv) ? zeros(0, 6) : reduce(vcat, [collect(t)' for t in pathv])

    HM = zeros(N); HM[Ito]   = EHM
    HR = zeros(N); HR[Ifrom] = EHR
    ucprob = zeros(N, N); ucprob[Ifrom, Ito] = cb.pp_ij
    rent = zeros(N); rent[Ito] = q_i[Ito]; rent[Ifrom] = Q_i[Ifrom]

    return Equilibrium(wage = wage_i, vv = vv, theta = theta_i, Y = Y, Q = Q_i, q = q_i,
                       HM = HM, HR = HR, rent = rent, A = copy(f.A), B = copy(f.B),
                       ucprob = ucprob, Phi = cb.Phi, Ubar = cb.Ubar, HH = HH,
                       converged = converged, iters = iters, path = path)
end

# =====================================================================================
# 10. ALGORITHMS 12, 13, 14 -- DECOMPOSITION AND RESERVATION UTILITY
#     cprod.m, cres.m, ubar.m
# =====================================================================================

"""
    decompose_productivity(p, data, A)  --  ALGORITHM 12, `cprod.m`

Eq. (20) / (S.55): split total adjusted productivity into a density-driven externality
and an exogenous fundamental,

    Upsilon_j = sum_s e^{-delta tau_js} (H_Ms / K_s) ,   a~_j = A~_j * Upsilon_j^(-lambda)

[D1] `delta` is the PRODUCTIVITY-spillover distance decay and `lambda` the density
elasticity of productivity. The inline comments in `cftualprep_end_TD.m:94-96` label
`delta` and `eta` with each other's meaning; the code there is right, the comments are
not, and the error is not propagated here.
"""
function decompose_productivity(p::ARSWParams, data::ARSWData, A::Vector{Float64})
    N    = data.N
    Iwpl = data.HM .!= 0
    dd   = exp.(-p.delta .* data.tau[Iwpl, Iwpl])
    EUps = dd * (data.HM[Iwpl] ./ data.K[Iwpl])
    Ups  = zeros(N);  Ups[Iwpl] = EUps
    a    = zeros(N);  a[Iwpl]   = A[Iwpl] ./ (EUps .^ p.lambda)
    return (a = a, Ups = Ups)
end

"""
    decompose_amenity(p, data, B)  --  ALGORITHM 13, `cres.m`

Eq. (21) / (S.56):
`Omega_i = sum_s e^{-rho tau_is} (H_Rs / K_s)`,  `b~_i = B~_i * Omega_i^(-eta)`.
[D1] `eta` is the density elasticity of AMENITY, `rho` its distance decay.
"""
function decompose_amenity(p::ARSWParams, data::ARSWData, B::Vector{Float64})
    N    = data.N
    Irsd = data.HR .!= 0
    rr   = exp.(-p.rho .* data.tau[Irsd, Irsd])
    EOme = rr * (data.HR[Irsd] ./ data.K[Irsd])
    Ome  = zeros(N);  Ome[Irsd] = EOme
    b    = zeros(N);  b[Irsd]   = B[Irsd] ./ (EOme .^ p.eta)
    return (b = b, Ome = Ome)
end

"""
    reservation_utility(p, B, Q, wage, tau, a, b)  --  ALGORITHM 14, `ubar.m`

Eq. (9): `Ubar = gamma * Phi^(1/eps)`, `gamma = Gamma((eps-1)/eps)`. Used as the target
in the OPEN-city solver (Alg. 16), where `Ubar` is pinned and total employment `H` is the
object that adjusts.
"""
function reservation_utility(p::ARSWParams, B::Vector{Float64}, Q::Vector{Float64},
                             wage::Vector{Float64}, tau::Matrix{Float64},
                             a::Vector{Float64}, b::Vector{Float64})
    Ito   = a .!= 0
    Ifrom = b .!= 0
    phi_ij = commuting_phi(p, tau[Ifrom, Ito], B[Ifrom], Q[Ifrom], wage[Ito])
    return gammaf_of(p) * sum(phi_ij) ^ (1 / p.epsilon)
end

# =====================================================================================
# 11. ALGORITHMS 15 AND 16 -- EQUILIBRIUM WITH ENDOGENOUS AGGLOMERATION
#     smodendog.m (closed city, Alg. 15)  and  ussmodendog.m (open city, Alg. 16)
# =====================================================================================

"""
    solve_equilibrium_endog(p, f, tau; closed = true, Utarget = nothing, ...)

ALGORITHMS 15 and 16, `smodendog.m` and `ussmodendog.m`, merged: the open-city solver is
the closed-city one plus a single extra update line for total employment, so upstream's
two 400-line near-duplicates collapse into one function with a `closed` keyword.

Here `f.A`, `f.B` hold the FUNDAMENTALS `{a~, b~}`; total productivity and amenity are
rebuilt from Eqs. (20)-(21) EVERY outer iteration from the current `{H_M, H_R}`:

    Upsilon_j = sum_s e^{-delta tau_js} (H_Ms/K_s) ,  A_j = a~_j Upsilon_j^lambda
    Omega_i   = sum_s e^{-rho tau_is}   (H_Rs/K_s) ,  B_i = b~_i Omega_i^eta

This recomputation is exactly what makes agglomeration endogenous relative to Alg. 11,
and is Section 7's main extra cost (two dense N x N products per pass).

* `closed = true`  (Alg. 15): `H` fixed, `Ubar` adjusts. With
  `lambda = delta = eta = rho = 0` this reduces to `solve_equilibrium_exog` exactly
  (check C5).
* `closed = false` (Alg. 16): `Ubar` pinned at `Utarget`, `H` adjusts through
  `H_up = (Ubar_model/Ubar_target)^eps * H`, then `H = 0.05*H_up + 0.95*H`. The weight is
  deliberately small: epsilon is large, so `Ubar/Utarget` close to one still implies big
  moves in `H`, and a larger weight makes the solver oscillate (upstream's own comment at
  `ussmodendog.m:338`).

Damping follows `smodendog.m`: 0.25 when the relevant max log gap is already below 0.1
(near convergence, avoid overshoot), 0.5 otherwise. Note that upstream's comment on that
branch reads "If we are far from convergence", which is the wrong way round; the code is
right.
"""
function solve_equilibrium_endog(p::ARSWParams, f::Fund, tau::Matrix{Float64};
                                 closed::Bool = true,
                                 Utarget::Union{Nothing,Float64} = nothing,
                                 maxiter::Int = 20_000, tol::Float64 = 1e-13,
                                 rule::Symbol = :tol, damp_far::Float64 = 0.5,
                                 damp_near::Float64 = 0.25, verbose::Bool = false)
    closed || Utarget !== nothing ||
        error("open-city solver (Algorithm 16) needs a reservation-utility target Utarget")

    N = length(f.A)
    ea_full, eb_full = f.A, f.B
    Ito   = ea_full .!= 0
    Ifrom = eb_full .!= 0
    IcsA  = Ito   .& .!Ifrom
    IcsB  = Ifrom .& .!Ito
    Iis   = Ito   .& Ifrom
    ea, eb = ea_full[Ito], eb_full[Ifrom]

    # [D3] exponent (1 - mu), not the literal 0.75 of smodendog.m:85
    L  = f.V .* (f.K .^ (1 - p.mu))
    HH = sum(f.HM)

    tau_sub = tau[Ifrom, Ito]
    dd_ij   = exp.(-p.delta .* tau[Ito,   Ito])       # spatial weights in Eq. (20)
    cc_ij   = exp.(-p.rho   .* tau[Ifrom, Ifrom])     # spatial weights in Eq. (21)
    EKM, EKR = f.K[Ito], f.K[Ifrom]

    Q_i = zeros(N);  Q_i[Ifrom] = f.Q[Ifrom]
    q_i = zeros(N);  q_i[Ito]   = f.Q[Ito]
    Q_e = copy(Q_i); q_e = copy(q_i)
    theta_i = copy(f.theta)
    Ewage_i = f.wage[Ito]
    EHM, EHR = f.HM[Ito], f.HR[Ifrom]

    EUps = dd_ij * (EHM ./ EKM);  EA = ea .* (EUps .^ p.lambda)      # Eq. (20)
    EOme = cc_ij * (EHR ./ EKR);  EB = eb .* (EOme .^ p.eta)         # Eq. (21)

    Y = zeros(N);  vv = copy(f.vv)
    wage_i = zeros(N); wage_e = zeros(N); theta_e = copy(theta_i)
    cb = nothing;  Ubar = NaN
    converged = false; iters = 0
    pathv = Vector{NTuple{6,Float64}}()

    for x in 1:maxiter
        iters = x
        phi_ij = commuting_phi(p, tau_sub, EB, Q_i[Ifrom], Ewage_i)
        cb     = commuting_block(p, phi_ij, HH)
        EHR, EHM = cb.HR, cb.HM
        Ubar = cb.Ubar

        # Agglomeration recomputed from the CURRENT employment distribution
        EUps = dd_ij * (EHM ./ EKM);  EA = ea .* (EUps .^ p.lambda)
        EOme = cc_ij * (EHR ./ EKR);  EB = eb .* (EOme .^ p.eta)

        Y[Ito]  = EA .* (EHM .^ p.alpha) .* ((theta_i[Ito] .* L[Ito]) .^ (1 - p.alpha))
        Ewage_e = (p.alpha .* Y[Ito]) ./ EHM
        vv[Ifrom] = (cb.pp_iji * Ewage_i) .* EHR

        q_e[IcsA] = ((1 - p.alpha) .* Y[IcsA]) ./ (theta_i[IcsA] .* L[IcsA])
        Q_e[IcsB] = ((1 - p.beta)  .* vv[IcsB]) ./ ((1 .- theta_i[IcsB]) .* L[IcsB])
        mixed_rent = (((1 - p.alpha) .* Y[Iis]) .+ ((1 - p.beta) .* vv[Iis])) ./ L[Iis]
        q_e[Iis] = mixed_rent
        Q_e[Iis] = mixed_rent

        theta_e = copy(theta_i)
        theta_e[Iis] = ((1 - p.alpha) .* Y[Iis]) ./ (Q_e[Iis] .* L[Iis])

        wage_i[Ito] = Ewage_i
        wage_e[Ito] = Ewage_e

        mLDw = maxlogdiff(wage_e, wage_i, Ito);  mLDq = maxlogdiff(q_e, q_i, Ito)
        mLDQ = maxlogdiff(Q_e, Q_i, Ifrom);      mLDt = maxlogdiff(theta_e, theta_i, Iis)
        mLDU = closed ? 0.0 : abs(log(Ubar / Utarget))
        mLD  = max(mLDw, mLDq, mLDQ, mLDt, mLDU)
        push!(pathv, (mLDw, mLDq, mLDQ, mLDt, mLD, Float64(x)))

        if rule === :round
            converged = rounded_equal(wage_e, wage_i) && rounded_equal(q_e, q_i) &&
                        rounded_equal(Q_e, Q_i)       && rounded_equal(theta_e, theta_i) &&
                        (closed || rounded_equal([Ubar], [Utarget]))
        else
            converged = mLD < tol
        end
        if converged
            verbose && println(">>>> Convergence Achieved <<<< (iter $x, maxLD $mLD)")
            break
        end

        # Adaptive damping (smodendog.m:293-304): 0.25 when close, 0.5 otherwise.
        cw = mLDw < 0.1 ? damp_near : damp_far
        cQ = mLDQ < 0.1 ? damp_near : damp_far
        ct = mLDt < 0.1 ? damp_near : damp_far
        Ewage_i = cw .* Ewage_e .+ (1 - cw) .* Ewage_i
        q_i     = cQ .* q_e     .+ (1 - cQ) .* q_i
        Q_i     = cQ .* Q_e     .+ (1 - cQ) .* Q_i
        theta_i = ct .* theta_e .+ (1 - ct) .* theta_i

        if !closed                                    # ALGORITHM 16 only
            HH_up = (Ubar / Utarget) ^ p.epsilon * HH
            HH    = 0.05 * HH_up + 0.95 * HH
        end
    end

    path = isempty(pathv) ? zeros(0, 6) : reduce(vcat, [collect(t)' for t in pathv])

    HM = zeros(N); HM[Ito]   = EHM
    HR = zeros(N); HR[Ifrom] = EHR
    A  = zeros(N); A[Ito]    = EA
    B  = zeros(N); B[Ifrom]  = EB
    ucprob = zeros(N, N); ucprob[Ifrom, Ito] = cb.pp_ij
    rent = zeros(N); rent[Ito] = q_i[Ito]; rent[Ifrom] = Q_i[Ifrom]

    return Equilibrium(wage = wage_i, vv = vv, theta = theta_i, Y = Y, Q = Q_i, q = q_i,
                       HM = HM, HR = HR, rent = rent, A = A, B = B, ucprob = ucprob,
                       Phi = cb.Phi, Ubar = Ubar, HH = HH,
                       converged = converged, iters = iters, path = path)
end

# =====================================================================================
# 12. ALGORITHMS 2 AND 3 -- ONE-STEP ESTIMATION OF EPSILON  (stretch goal)
#     cdensityoptren.m (Alg. 2)  and  optimepsilon_TD86.m (Alg. 3)
# =====================================================================================

"""
    epsilon_objective(eps_try, omega, HM, district, var_log_wage_data)

ALGORITHM 2, `cdensityoptren.m`. Objective for the one-step estimation of epsilon.

The identification (paper p. 2167, supplement Eq. S.64): transformed wages `omega_j` are
determined by commuting market clearing alone, independently of epsilon; epsilon only
scales the map `w~_j = omega_j^(1/eps)`, so
`var(log w~) = (1/eps)^2 var(log omega)`. Matching the model's variance of log DISTRICT
wages to the variance observed in district wage data therefore pins epsilon down.

Given a trial epsilon: convert omega to adjusted wages, normalise to geometric mean one,
form the block wage bill `w~_j H_Mj`, aggregate to districts (district wage = district
wage bill / district employment), take the variance of demeaned log district wages, and
return the squared deviation from the data moment, scaled by 1e6 (as upstream does, to
keep the residual sum of squares away from numerical noise).
"""
function epsilon_objective(eps_try::Float64, omega::Vector{Float64}, HM::Vector{Float64},
                           district::Vector{Int}, var_log_wage_data::Float64)
    wage = omega .^ (1 / eps_try)
    pos  = wage .> 0
    wage[pos] ./= geomean_pos(wage[pos])
    wbill = wage .* HM

    ids = sort(unique(district))
    dw  = Float64[]
    for d in ids
        m = district .== d
        emp = sum(HM[m])
        emp > 0 && push!(dw, sum(wbill[m]) / emp)
    end
    lbw = log.(dw);  lbw .-= mean(lbw)
    var_mod = var(lbw)
    return ((var_mod - var_log_wage_data) ^ 2) * 1e6
end

"""
    estimate_epsilon(p, data, district, var_log_wage_data; ...)

ALGORITHM 3, `optimepsilon_TD86.m`. Upstream uses MATLAB's `patternsearch` (Global
Optimization Toolbox) over `eps in [2, 24]` from a start of 4. The objective is
one-dimensional and smooth, so here `Optim.Brent()` over the same bracket is used, with a
`NelderMead` run from the same start as a cross-check that the two agree.

Note that Algorithm 1 is run ONCE, before the search: omega does not depend on epsilon,
which is the whole point of the "one-step" procedure.

Upstream rounds the estimate to two decimals afterwards (`epsilon=round(epsilon*100)/100`,
giving the paper's 6.83). That rounding is reported but not applied to the returned value,
so the recovery error is visible.
"""
function estimate_epsilon(p::ARSWParams, data::ARSWData, district::Vector{Int},
                          var_log_wage_data::Float64;
                          lower::Float64 = 2.0, upper::Float64 = 24.0,
                          start::Float64 = 4.0, verbose::Bool = true)
    # ALGORITHM 1: transformed wages, independent of epsilon.
    om = solve_transformed_wages(p, data; tol = 1e-14, maxiter = 50_000)
    verbose && @printf("    Algorithm 1 (omega): converged=%s, iters=%d, gap=%.3e\n",
                       om.converged, om.iters, om.gap)

    obj(e) = epsilon_objective(e, copy(om.omega), data.HM, district, var_log_wage_data)
    res_b  = Optim.optimize(obj, lower, upper, Optim.Brent(); abs_tol = 1e-12)
    eps_b  = Optim.minimizer(res_b)
    res_n  = Optim.optimize(x -> obj(clamp(x[1], lower, upper)), [start], Optim.NelderMead(),
                            Optim.Options(g_tol = 1e-14, iterations = 10_000))
    eps_n  = clamp(Optim.minimizer(res_n)[1], lower, upper)

    return (epsilon = eps_b, epsilon_neldermead = eps_n, objective = Optim.minimum(res_b),
            omega = om.omega)
end

# =====================================================================================
# 13. DATA LAYER -- THE ONE SEAM BETWEEN SYNTHETIC AND REAL
# =====================================================================================

"""
    load_berlin_data(path) -> ARSWData

THE SINGLE FUNCTION TO REPLACE to run this file on the real Berlin data. Everything
downstream takes an `ARSWData` and does not care where it came from.

Expected inputs (external downloads, see the toolkit README; NOT fetched by this script):

    <path>/prepdata_big_TD.mat        2006 cross-section
        floor06    nobs06 x 1   floor-space price        -> ARSWData.Q
        empwpl06   nobs06 x 1   workplace employment     -> ARSWData.HM
        emprsd06   nobs06 x 1   residence employment     -> ARSWData.HR
        area06     nobs06 x 1   block area               -> ARSWData.K
        tt06       nobs06^2     bilateral travel time    -> ARSWData.tau   (minutes)
        dummywestr nobs06 x 1   1 = former West Berlin   -> ARSWData.west
        Xr, Yr     nobs06 x 1   block centroids          -> ARSWData.coords
        nobs06     scalar       number of blocks         -> ARSWData.N

    <path>/prepdata_big_TD86.mat      1986 cross-section, same schema with `86rw`
                                      suffixes, West Berlin only; needed only for the
                                      epsilon estimation (Algorithms 1-3), together with
                                      `wageworker1986.csv` (12 Bezirke wages) and the
                                      block-to-Bezirk key `bzk86rw`.
    <path>/ttpublic_2006_ren.mat      public-transport-only travel times `ttpub06`,
                                      needed only for the car-ban counterfactual.

Implementation is four lines with `MAT.jl`:

    using MAT
    v = matread(joinpath(path, "prepdata_big_TD.mat"))
    ARSWData(Int(v["nobs06"]), vec(v["floor06"]), vec(v["empwpl06"]), vec(v["emprsd06"]),
             vec(v["area06"]), Matrix(v["tt06"]), vec(v["dummywestr"]) .== 1,
             hcat(vec(v["Xr"]), vec(v["Yr"])))

MAT.jl is deliberately NOT a dependency of this file, because the repository has no
Project.toml and the .mat files are not present.

Index alignment is load-bearing upstream: the MATLAB arrays and the `Berlin4matlab`
shapefile are sorted identically. That matters for mapping, which is out of scope here.
"""
function load_berlin_data(path::AbstractString)
    error("""
          load_berlin_data is a stub: the real Berlin .mat files are not in this
          repository and are deliberately not downloaded.

          Expected at $path :
            prepdata_big_TD.mat, prepdata_big_TD86.mat, ttpublic_2006_ren.mat

          See this function's docstring for the full variable -> field map, then
          implement it with MAT.jl (four lines). Nothing else in this file needs to
          change: every algorithm takes an ARSWData.
          """)
end

# =====================================================================================
# 14. SYNTHETIC BERLIN-LIKE GEOGRAPHY, WITH A WALL
# =====================================================================================

"""
    build_grid_geography(; nx, ny, block_km, speed, tau_own, wall_col, wall_penalty,
                           divided, seed)

A Berlin-like block geography: an `nx x ny` grid of blocks on a `block_km` lattice, a CBD
at the centre, block land area falling towards the CBD (finer subdivision downtown), and
bilateral travel times in MINUTES -- the unit in which `nu = 0.07` is calibrated.

THE WALL. A vertical cut between columns `wall_col` and `wall_col+1` splits the city into
West (`col <= wall_col`) and East. With the defaults (`nx = 15`, `wall_col = 7`) the
CBD sits in column 8, i.e. just inside the EAST -- deliberately, because that is
Berlin's situation: the historic centre (Mitte) fell on the eastern side of the wall,
which is why division hit West Berlin's access to the old core hardest. Under `divided = true` every origin-destination pair
that straddles the cut gets `tau_ij += wall_penalty`. With `nu = 0.07` and the default
penalty of 45 minutes that multiplies `phi_ij` by `exp(-3.15) ~ 0.043`: cross-wall
commuting falls by ~96%, a severe barrier, while `tau` stays finite so no `Inf`/`0` ever
enters a log or an exponential kernel.

WHY NOT A LITERAL SEVERING. Raising the penalty further is numerically, not
economically, costly. Commuting market clearing, Eq. (S.44), pins transformed wages only
up to ONE global scale, and the proof that the fixed point is unique (Lemmas S.6-S.7,
gross substitution) needs the commuting matrix to be irreducible. A hard wall makes it
nearly reducible: the two halves' wage LEVELS are then tied to each other only through
commuting flows of order `exp(-nu*penalty)`, so the slow mode of the iteration has
modulus `1 - O(exp(-nu*penalty))` and the iteration count scales like
`exp(+nu*penalty)`. Measured on this geography (`run_wall_conditioning`), the inversion
needs 58 iterations with no wall, 218 at 30 minutes, 984 at 45, 2,668 at 60 and more
than 100,000 at 120 -- all converging to machine precision, so this is conditioning, not
correctness. It is also the reason ARSW estimate epsilon on the WEST-ONLY 1986 sample
rather than on the divided city as a whole. This is the synthetic analogue of the paper's identifying variation -- the same
fundamentals under two different travel-time matrices.

Returns a NamedTuple `(N, coords, K, tau, west, dist_cbd, col, row)`.
"""
function build_grid_geography(; nx::Int = 15, ny::Int = 15, block_km::Float64 = 1.0,
                              speed::Float64 = 0.5, tau_own::Float64 = 1.5,
                              wall_col::Int = 7, wall_penalty::Float64 = 45.0,
                              divided::Bool = false, seed::Int = 20150615,
                              area_scale::Float64 = 0.25, area_gradient::Float64 = 1.0,
                              area_sd::Float64 = 0.15)
    rng = MersenneTwister(seed)
    N = nx * ny
    coords = zeros(N, 2)
    col = zeros(Int, N);  row = zeros(Int, N)
    for c in 1:nx, r in 1:ny
        i = (c - 1) * ny + r
        coords[i, 1] = (c - 1) * block_km
        coords[i, 2] = (r - 1) * block_km
        col[i] = c;  row[i] = r
    end
    cbd = [(nx - 1) * block_km / 2, (ny - 1) * block_km / 2]
    dist_cbd = [hypot(coords[i,1] - cbd[1], coords[i,2] - cbd[2]) for i in 1:N]
    dmax = maximum(dist_cbd)

    # Land area: smaller blocks downtown, lognormal noise.
    K = area_scale .* exp.(area_gradient .* (dist_cbd ./ dmax) .+
                           area_sd .* randn(rng, N))

    tau = zeros(N, N)
    for i in 1:N, j in 1:N
        d = hypot(coords[i,1] - coords[j,1], coords[i,2] - coords[j,2])
        tau[i, j] = i == j ? tau_own : d / speed
    end
    west = col .<= wall_col
    if divided
        for i in 1:N, j in 1:N
            west[i] != west[j] && (tau[i, j] += wall_penalty)
        end
    end
    return (N = N, coords = coords, K = K, tau = tau, west = west,
            dist_cbd = dist_cbd, col = col, row = row, nx = nx, ny = ny)
end

"""
    plant_fundamentals(geo; scenario, ...)

Draws the fundamentals that the recovery test must get back.

* `log a_j` rises with distance from the CBD -- matching the paper's finding that TOTAL
  productivity falls with distance while its EXOGENOUS component rises (Section 7.1);
* `log b_i` likewise rises with distance;
* `log phi_i` (density of development) falls with distance.

Shocks are independent standard normals, so `a` and `b` are not mechanically collinear.

`scenario`:
* `:mixed` -- every block has `a_j > 0` and `b_i > 0`, so every block is incompletely
  specialised, `q_i = Q_i`, and the single observed floor price is unambiguous. This is
  the benchmark case of supplement S.2.6.
* `:specialised` -- a central 3x3 office district gets `b_i = 0` (theta = 1, H_R = 0) and
  two outer patches get `a_j = 0` (theta = 0, H_M = 0). This exercises the `IcsA`/`IcsB`
  index sets that `:mixed` never touches.

`floorspace_per_worker` fixes the units of floor space (and therefore the level of the
floor price). It is a pure unit choice: scaling all `L` uniformly leaves
`{H_M, H_R, theta, pi_ij}` unchanged.
"""
function plant_fundamentals(geo; scenario::Symbol = :mixed, seed::Int = 20150615,
                            g_a::Float64 = 0.35, sd_a::Float64 = 0.25,
                            g_b::Float64 = 0.45, sd_b::Float64 = 0.25,
                            g_phi::Float64 = -1.10, sd_phi::Float64 = 0.20,
                            H::Float64 = 1.5e6, floorspace_per_worker::Float64 = 0.05,
                            p::ARSWParams = PARAMS_SECTION6)
    rng = MersenneTwister(seed + 1)
    N = geo.N
    dn = geo.dist_cbd ./ maximum(geo.dist_cbd)

    a   = exp.(g_a   .* dn .+ sd_a   .* randn(rng, N))
    b   = exp.(g_b   .* dn .+ sd_b   .* randn(rng, N))
    phi = exp.(g_phi .* dn .+ sd_phi .* randn(rng, N))

    if scenario === :specialised
        cx = (geo.nx + 1) ÷ 2;  cy = (geo.ny + 1) ÷ 2
        office = (abs.(geo.col .- cx) .<= 1) .& (abs.(geo.row .- cy) .<= 1)
        b[office] .= 0.0                                        # purely commercial (IcsA)
        park = ((geo.col .== 3) .& (abs.(geo.row .- cy) .<= 2)) .|
               ((geo.col .== geo.nx - 2) .& (abs.(geo.row .- cy) .<= 2))
        a[park] .= 0.0                                          # purely residential (IcsB)
    end

    # Put floor space on the requested per-worker scale (a unit choice, see docstring).
    L_target = floorspace_per_worker * H
    Lraw = phi .* (geo.K .^ (1 - p.mu))
    phi .*= L_target / sum(Lraw)

    # theta starting values that are CORRECT by construction for specialised blocks:
    # purely commercial -> 1, purely residential -> 0, mixed -> interior guess.
    theta0 = fill(0.5, N)
    theta0[(a .!= 0) .& (b .== 0)] .= 1.0
    theta0[(b .!= 0) .& (a .== 0)] .= 0.0
    return (a = a, b = b, phi = phi, theta0 = theta0, H = H)
end

"""
    forward_fund(geo, pl, p; Q0)

Assemble the `Fund` bundle for a forward (data-generating) solve from planted
fundamentals, with deliberately naive starting guesses: uniform employment, uniform floor
price, and wages from the labour FOC at that guess. Upstream notes that seeding the
solver with the observed equilibrium speeds it up a lot; starting from uniform values is
the more demanding test, and is what is done here.
"""
function forward_fund(geo, pl, p::ARSWParams; Q0::Float64 = 1.0)
    N = geo.N
    L = pl.phi .* (geo.K .^ (1 - p.mu))
    HM0 = fill(pl.H / N, N);  HM0[pl.a .== 0] .= 0.0
    HR0 = fill(pl.H / N, N);  HR0[pl.b .== 0] .= 0.0
    HM0 .*= pl.H / sum(HM0);  HR0 .*= pl.H / sum(HR0)

    w0 = zeros(N)
    Ito = pl.a .!= 0
    w0[Ito] = p.alpha .* pl.a[Ito] .*
              ((pl.theta0[Ito] .* L[Ito] ./ HM0[Ito]) .^ (1 - p.alpha))
    # vv is overwritten at the top of every solver pass, so this is only a placeholder.
    wbar = mean(w0[Ito])
    vv0  = zeros(N);  vv0[pl.b .!= 0] = wbar .* HR0[pl.b .!= 0]
    return Fund(A = copy(pl.a), B = copy(pl.b), V = copy(pl.phi), K = copy(geo.K),
                Q = fill(Q0, N), HM = HM0, HR = HR0, LD = L, theta = copy(pl.theta0),
                wage = w0, vv = vv0)
end

"Turn a solved forward equilibrium into the 'observed' dataset the inversion will see."
function equilibrium_to_data(geo, eq::Equilibrium)
    return ARSWData(geo.N, copy(eq.rent), copy(eq.HM), copy(eq.HR), copy(geo.K),
                    copy(geo.tau), copy(geo.west), copy(geo.coords))
end

"""
    make_synthetic_data(; endogenous, scenario, divided, p, ...)

Generate the synthetic "observed" data by solving the model FORWARD from planted
fundamentals, and return both, on the normalisations the inversion imposes.

Two exact invariances of the model make the normalisation step legitimate (both are
verified numerically in check C7, not merely asserted):

* scaling every `a_j` by `k` scales `{w, q, Q, Y, vv}` by `k` and leaves
  `{H_M, H_R, theta, pi_ij}` unchanged;
* scaling every `b_i` by `c` scales `Phi` by `c^eps` and leaves the ENTIRE allocation
  unchanged.

So the generator solves once, rescales `a <- a / geomean(A_eq)` and
`b <- b * (H/Phi)^(1/eps)`, and re-solves. The second solve produces the same allocation
with `geomean(A~) = 1` and `Phi = H` -- exactly the normalisations that
`calcal_adj_TD.m` and `cmodexog.m` impose on the RECOVERED objects, so planted and
recovered are then directly comparable without any post-hoc rescaling.

With `endogenous = true` the forward solve uses Algorithm 15 and the planted `{a, b}` are
FUNDAMENTALS; total `{A, B}` come out of Eqs. (20)-(21) at the equilibrium.
"""
function make_synthetic_data(; endogenous::Bool = false, scenario::Symbol = :mixed,
                             divided::Bool = false, p::ARSWParams = PARAMS_SECTION6,
                             geo = nothing, tol::Float64 = 1e-14, verbose::Bool = false,
                             kwargs...)
    g  = geo === nothing ? build_grid_geography(divided = divided) : geo
    pl = plant_fundamentals(g; scenario = scenario, p = p, kwargs...)

    # DAMP_ROBUST, not upstream's 0.5: the cold start is a large shock, and on the divided
    # geography the 50/50 blend cycles rather than converging (see DAMP_ROBUST's docstring).
    solve = (f) -> endogenous ?
        solve_equilibrium_endog(p, f, g.tau; closed = true, tol = tol, maxiter = 200_000,
                                damp_far = DAMP_ROBUST, damp_near = DAMP_ROBUST,
                                verbose = verbose) :
        solve_equilibrium_exog(p, f, g.tau; tol = tol, maxiter = 200_000,
                               damp = DAMP_ROBUST, verbose = verbose)

    f_cold = forward_fund(g, pl, p)
    eq_cold = solve(f_cold)
    eq_cold.converged || @warn "forward cold solve did not converge" iters=eq_cold.iters

    # Impose the inversion's two normalisations on the PLANTED fundamentals.
    #
    # They are NOT independent. Rescaling a by 1/kA scales {w, q, Q} by 1/kA, and since
    # phi_ij carries Q^(-(1-beta)eps) w^eps, that multiplies Phi by kA^(-beta*eps). So the
    # amenity rescaling has to be applied AFTER the productivity one, using the Phi that
    # the productivity-rescaled economy actually produces. (Getting this wrong leaves B
    # off by a factor of kA^beta -- a 0.16 log error that the recovery test catches
    # immediately, which is the point of running the test.)
    a_cur = copy(pl.a);  b_cur = copy(pl.b)
    L     = pl.phi .* (g.K .^ (1 - p.mu))
    eq    = eq_cold
    for _ in 1:4
        kA = geomean_pos(eq.A)                             # -> geomean(A~) = 1
        a_cur ./= kA
        eq = solve(Fund(A = a_cur, B = b_cur, V = pl.phi, K = copy(g.K),
                        Q = eq.rent ./ kA, HM = copy(eq.HM), HR = copy(eq.HR), LD = L,
                        theta = copy(eq.theta), wage = eq.wage ./ kA, vv = eq.vv ./ kA))
        kB = (pl.H / eq.Phi) ^ (1 / p.epsilon)             # -> Phi = H
        b_cur .*= kB
        eq = solve(Fund(A = a_cur, B = b_cur, V = pl.phi, K = copy(g.K),
                        Q = copy(eq.rent), HM = copy(eq.HM), HR = copy(eq.HR), LD = L,
                        theta = copy(eq.theta), wage = copy(eq.wage), vv = copy(eq.vv)))
        (abs(log(geomean_pos(eq.A))) < 1e-13 && abs(log(eq.Phi / pl.H)) < 1e-13) && break
    end
    eq.converged || @warn "forward solve did not converge" iters=eq.iters
    f2 = Fund(A = a_cur, B = b_cur, V = pl.phi, K = copy(g.K), Q = copy(eq.rent),
              HM = copy(eq.HM), HR = copy(eq.HR), LD = L, theta = copy(eq.theta),
              wage = copy(eq.wage), vv = copy(eq.vv))
    a2, b2 = a_cur, b_cur

    data = equilibrium_to_data(g, eq)
    planted = (a = a2, b = b2, phi = pl.phi, H = pl.H,
               A = eq.A, B = eq.B, theta = eq.theta, LD = f2.LD)
    return (data = data, planted = planted, eq = eq, eq_cold = eq_cold, geo = g, fund = f2)
end

# =====================================================================================
# 15. QUANTIFICATION DRIVERS  (calcal_TD.m and cftualprep_TD.m, inlined)
# =====================================================================================

"""
    quantify_sequential(p, data; ...)  --  driver `calcal_TD.m`, Algorithms 4,5,6,7,8

Runs the sequential quantification in upstream's order:
  Alg. 4 wages + productivity -> Alg. 5 amenities -> Alg. 6 rescale -> Alg. 7 income ->
  Alg. 8 floor space, land-use share, density of development.

Algorithm 6 is essential: without it the `{A~, B~}` from Algorithms 4-5 are valid only
for RELATIVE cross-location comparisons and are not on the scale the equilibrium solvers
need (upstream's own warning, toolkit README / codebook p. 283).
"""
function quantify_sequential(p::ARSWParams, data::ARSWData; tol::Float64 = 1e-14,
                             rule::Symbol = :tol, verbose::Bool = false)
    w = solve_wages_and_productivity(p, data; tol = tol, rule = rule, maxiter = 200_000,
                                     verbose = verbose)
    am = solve_amenities(p, data, w.wage)                               # Alg. 5
    rs = rescale_fundamentals(p, data, w.A, am.B)                       # Alg. 6
    ei = expected_income(p, data, rs.wage, rs.B)                        # Alg. 7
    de = solve_density(p, data, rs.A, ei.vv)                            # Alg. 8
    return (A = rs.A, B = rs.B, wage = rs.wage, vv = ei.vv, CMA = am.CMA,
            V = de.V, LD = de.LD, LM = de.LM, LR = de.LR, theta = de.theta,
            A_step4 = w.A, B_step5 = am.B, wage_step4 = w.wage, cprob = w.cprob,
            converged = w.converged, iters = w.iters, gap = w.gap, Phi = ei.Phi)
end

"""
    quantify_simultaneous(p, data; ...)  --  driver `cftualprep_TD.m`, Algorithms 9, 10

The quantification the paper uses to seed counterfactuals.
"""
function quantify_simultaneous(p::ARSWParams, data::ARSWData; tol::Float64 = 1e-14,
                               rule::Symbol = :tol, verbose::Bool = false)
    inv = invert_simultaneous(p, data; tol = tol, rule = rule, maxiter = 200_000,
                              verbose = verbose)
    de  = solve_density(p, data, inv.A, inv.vv)                         # Alg. 10
    return (A = inv.A, B = inv.B, wage = inv.wage, vv = inv.vv, CMA = inv.CMA,
            V = de.V, LD = de.LD, LM = de.LM, LR = de.LR, theta = de.theta,
            ucprob = inv.ucprob, Phi = inv.Phi, HH = inv.HH,
            converged = inv.converged, iters = inv.iters,
            mAgap = inv.mAgap, mBgap = inv.mBgap)
end

"Build the `Fund` bundle that seeds an equilibrium solve from a quantification result."
function fund_from_quantification(data::ARSWData, qn; endog = nothing)
    A = endog === nothing ? qn.A : endog.a
    B = endog === nothing ? qn.B : endog.b
    return Fund(A = A, B = B, V = qn.V, K = copy(data.K), Q = copy(data.Q),
                HM = copy(data.HM), HR = copy(data.HR), LD = qn.LD,
                theta = qn.theta, wage = qn.wage, vv = qn.vv)
end

# =====================================================================================
# 16. VALIDATION
# =====================================================================================

const RESULTS = Vector{NamedTuple{(:id, :label, :value, :target, :pass),
                                  Tuple{String,String,Float64,Float64,Bool}}}()

function report(id, label, value, target; smaller_is_better = true)
    pass = smaller_is_better ? (value <= target) : (value >= target)
    push!(RESULTS, (id = id, label = label, value = value, target = target, pass = pass))
    @printf("  %-4s %-58s %12.4e  (target %s %.0e)  %s\n", id, label, value,
            smaller_is_better ? "<=" : ">=", target, pass ? "PASS" : "FAIL")
    return pass
end

function report_note(id, label, value)
    @printf("  %-4s %-58s %12.4e\n", id, label, value)
end

header(s) = (println(); println("="^100); println(s); println("="^100))

"Print a multi-line note at a fixed indent (Julia triple-quoted strings strip the common
leading whitespace, which otherwise flattens these against the left margin)."
note(txt::AbstractString; indent = "       ") =
    for ln in split(rstrip(txt), '\n'); println(indent, lstrip(ln)); end

"""
    run_recovery_tests(p6, p7)

The primary validation: plant fundamentals, solve forward, invert, compare. Tolerances
and test ids follow TRANSLATION_PLAN.md section 8. All errors reported as
`max_i |log(recovered_i / planted_i)|`, a scale-free relative error.
"""
function run_recovery_tests(p6::ARSWParams, p7::ARSWParams)
    header("RECOVERY TESTS -- does the inversion return the planted fundamentals?")

    # ---------------- R1: forward solve is a genuine equilibrium --------------------
    println("\n-- R1: forward solve residual (exogenous fundamentals, :mixed) --")
    syn = make_synthetic_data(; scenario = :mixed, p = p6)
    data, planted, eq = syn.data, syn.planted, syn.eq
    report("R1", "forward-solve max log gap over {w,q,Q,theta}", syn.eq_cold.path[end, 5], 1e-12)
    @printf("       blocks=%d; COLD solve from uniform guesses took %d iterations;\n",
            data.N, syn.eq_cold.iters)
    @printf("       after renormalisation Phi/H=%.15f, geomean(A~)=%.15f\n",
            eq.Phi / planted.H, geomean_pos(eq.A))
    @printf("       mean floor price=%.4f   mean wage=%.4f   mean theta=%.4f\n",
            mean(data.Q), mean(filter(>(0.0), eq.wage)), mean(eq.theta))

    # ---------------- R2: sequential inversion --------------------------------------
    println("\n-- R2: SEQUENTIAL inversion (Algorithms 4->5->6->7->8) --")
    seq = quantify_sequential(p6, data)
    @printf("       Alg. 4 converged=%s in %d iterations (gap %.3e)\n",
            seq.converged, seq.iters, seq.gap)
    report("R2a", "sequential: adjusted productivity A~", maxlogdiff(seq.A, planted.A), 1e-8)
    report("R2b", "sequential: adjusted amenity B~",      maxlogdiff(seq.B, planted.B), 1e-8)
    report("R2c", "sequential: density of development phi", maxlogdiff(seq.V, planted.phi), 1e-8)
    report("R2d", "sequential: land-use share theta",     maxlogdiff(seq.theta, planted.theta), 1e-8)

    # ---------------- R3: simultaneous inversion ------------------------------------
    println("\n-- R3: SIMULTANEOUS inversion (Algorithms 9->10) --")
    sim = quantify_simultaneous(p6, data)
    @printf("       Alg. 9 converged=%s in %d iterations (H_M gap %.3e, H_R gap %.3e)\n",
            sim.converged, sim.iters, sim.mAgap, sim.mBgap)
    report("R3a", "simultaneous: adjusted productivity A~", maxlogdiff(sim.A, planted.A), 1e-8)
    report("R3b", "simultaneous: adjusted amenity B~",      maxlogdiff(sim.B, planted.B), 1e-8)
    report("R3c", "simultaneous: density of development phi", maxlogdiff(sim.V, planted.phi), 1e-8)

    # ---------------- R4: the two procedures agree ----------------------------------
    println("\n-- R4: sequential vs simultaneous (the cross-check cftualprep_TD.m builds in) --")
    report("R4a", "A~ sequential vs simultaneous", maxlogdiff(seq.A, sim.A), 1e-8)
    report("R4b", "B~ sequential vs simultaneous", maxlogdiff(seq.B, sim.B), 1e-8)
    report("R4c", "phi sequential vs simultaneous", maxlogdiff(seq.V, sim.V), 1e-8)

    # ---------------- R5: specialised blocks ----------------------------------------
    println("\n-- R5: :specialised scenario (exercises the IcsA/IcsB/Iis index sets) --")
    syn_s = make_synthetic_data(; scenario = :specialised, p = p6)
    ds, ps, es = syn_s.data, syn_s.planted, syn_s.eq
    nA = count((ps.a .!= 0) .& (ps.b .== 0));  nB = count((ps.b .!= 0) .& (ps.a .== 0))
    @printf("       completely specialised: %d commercial (theta=1), %d residential (theta=0); mixed %d\n",
            nA, nB, count((ps.a .!= 0) .& (ps.b .!= 0)))
    @printf("       forward cold solve: max log gap %.3e in %d iterations\n",
            syn_s.eq_cold.path[end,5], syn_s.eq_cold.iters)
    # coherence of cdensity's single observed floor price (see solve_density docstring)
    mixed = (ps.a .!= 0) .& (ps.b .!= 0)
    report_note("R5*", "max |q-Q| on incompletely specialised blocks",
                maximum(abs.(es.q[mixed] .- es.Q[mixed])))
    seq_s = quantify_sequential(p6, ds)
    sim_s = quantify_simultaneous(p6, ds)
    report("R5a", "specialised, sequential: A~",   maxlogdiff(seq_s.A, ps.A), 1e-8)
    report("R5b", "specialised, sequential: B~",   maxlogdiff(seq_s.B, ps.B), 1e-8)
    report("R5c", "specialised, sequential: phi",  maxlogdiff(seq_s.V, ps.phi), 1e-8)
    report("R5d", "specialised, simultaneous: A~", maxlogdiff(sim_s.A, ps.A), 1e-8)
    report("R5e", "specialised, simultaneous: B~", maxlogdiff(sim_s.B, ps.B), 1e-8)

    # ---------------- R6: Section 7, fundamentals behind the externalities ----------
    println("\n-- R6: SECTION 7 -- endogenous agglomeration, recover {a~, b~} --")
    @printf("       lambda=%.4f delta=%.4f eta=%.4f rho=%.4f  eps=%.3f nu=%.4f\n",
            p7.lambda, p7.delta, p7.eta, p7.rho, p7.epsilon, p7.kappaeps)
    syn7 = make_synthetic_data(; endogenous = true, scenario = :mixed, p = p7)
    d7, pl7, eq7 = syn7.data, syn7.planted, syn7.eq
    @printf("       forward cold solve (Alg. 15): max log gap %.3e in %d iterations\n",
            syn7.eq_cold.path[end,5], syn7.eq_cold.iters)
    sim7 = quantify_simultaneous(p7, d7)                       # Alg. 9 recovers TOTAL A~,B~
    report("R6a", "Section 7: total adjusted productivity A~", maxlogdiff(sim7.A, pl7.A), 1e-8)
    report("R6b", "Section 7: total adjusted amenity B~",      maxlogdiff(sim7.B, pl7.B), 1e-8)
    dp = decompose_productivity(p7, d7, sim7.A)                # Alg. 12
    dr = decompose_amenity(p7, d7, sim7.B)                     # Alg. 13
    report("R6c", "Section 7: production fundamental a~",  maxlogdiff(dp.a, pl7.a), 1e-8)
    report("R6d", "Section 7: residential fundamental b~", maxlogdiff(dr.b, pl7.b), 1e-8)
    # The economics the decomposition is for: total vs exogenous gradient (paper Sec. 7.1)
    g  = syn7.geo
    slope(y) = (X = hcat(ones(g.N), g.dist_cbd); (X \ log.(y))[2])
    @printf("       gradient per km: total A %+.4f, exogenous a %+.4f (wedge %+.4f)\n",
            slope(sim7.A), slope(dp.a), slope(sim7.A) - slope(dp.a))
    @printf("       gradient per km: total B %+.4f, exogenous b %+.4f (wedge %+.4f)\n",
            slope(sim7.B), slope(dr.b), slope(sim7.B) - slope(dr.b))
    note("""The BASE gradients are a property of the synthetic draw and mean nothing. The
       WEDGE does: agglomeration makes total productivity fall away from the centre
       faster than its fundamental, because Upsilon is highest where employment is
       densest. Berlin's wedge (paper Sec. 7.1) is -1.9%/km (total -1.7, exogenous +0.2);
       here it is the number printed above, from lambda=0.0710 and delta=0.3617 acting on
       a synthetic density gradient. Same mechanism, same order of magnitude; the levels
       are not comparable and are not claimed to be.""")

    # ---------------- R7: open city reproduces the closed city ----------------------
    println("\n-- R7: OPEN city (Alg. 16) reproduces CLOSED city (Alg. 15) at its own Ubar --")
    Ubar_target = reservation_utility(p7, eq7.B, eq7.Q, eq7.wage, d7.tau, pl7.a, pl7.b)  # Alg. 14
    @printf("       Alg. 14 Ubar from ubar.m formula: %.10e ; solver's Ubar: %.10e ; rel.diff %.3e\n",
            Ubar_target, eq7.Ubar, abs(Ubar_target / eq7.Ubar - 1))
    # Start DELIBERATELY away from the answer -- 35%% too few workers and prices 20%% off --
    # so that the H update (the one line that distinguishes ussmodendog.m from smodendog.m)
    # actually has to do work. Seeded at the solution the test would be vacuous.
    f7 = Fund(A = pl7.a, B = pl7.b, V = pl7.phi, K = copy(d7.K), Q = eq7.Q .* 1.2,
              HM = eq7.HM .* 0.65, HR = eq7.HR .* 0.65, LD = pl7.LD,
              theta = copy(pl7.theta), wage = eq7.wage .* 1.2, vv = eq7.vv .* 0.8)
    @printf("       starting from H = %.2f%% of the closed-city H, Q 20%% high\n",
            100 * sum(f7.HM) / eq7.HH)
    eq_open = solve_equilibrium_endog(p7, f7, d7.tau; closed = false, Utarget = eq7.Ubar,
                                      tol = 1e-11, maxiter = 200_000,
                                      damp_far = DAMP_ROBUST, damp_near = DAMP_ROBUST)
    @printf("       open-city solver converged=%s in %d iterations; H recovered to %.10f of target\n",
            eq_open.converged, eq_open.iters, eq_open.HH / eq7.HH)
    report("R7a", "open vs closed city: workplace employment H_M",
           maxlogdiff(eq_open.HM, eq7.HM), 1e-6)
    report("R7b", "open vs closed city: total employment H",
           abs(log(eq_open.HH / eq7.HH)), 1e-6)
    report("R7c", "open vs closed city: floor price Q", maxlogdiff(eq_open.Q, eq7.Q), 1e-6)
    report("R7d", "open vs closed city: Ubar hits its target",
           abs(log(eq_open.Ubar / eq7.Ubar)), 1e-10)

    return (syn = syn, seq = seq, sim = sim, syn7 = syn7, sim7 = sim7, dp = dp, dr = dr,
            syn_s = syn_s, eq_open = eq_open)
end

"""
    run_model_checks(p6, p7, R)

Everything on the validation checklist that is not a fundamentals-recovery test:
market clearing, the zero-shock counterfactual, the lambda=delta=eta=rho=0 nesting,
comparative statics, the two scale invariances, and the inertness of mu.
"""
function run_model_checks(p6::ARSWParams, p7::ARSWParams, R)
    header("MODEL CHECKS")
    data, planted, eq = R.syn.data, R.syn.planted, R.syn.eq
    geo = R.syn.geo
    sim = R.sim

    # ---------------- C1: commuting probabilities ------------------------------------
    println("\n-- C1: commuting probabilities --")
    Iwpl = data.HM .!= 0;  Irsd = data.HR .!= 0
    phi_ij = commuting_phi(p6, data.tau[Irsd, Iwpl], sim.B[Irsd], data.Q[Irsd], sim.wage[Iwpl])
    cb = commuting_block(p6, phi_ij, sum(data.HM))
    report("C1a", "max |sum_j pi_{j|i} - 1| over residences",
           maximum(abs.(vec(sum(cb.pp_iji, dims = 2)) .- 1)), 1e-13)
    report("C1b", "max |sum_i pi_{i|j} - 1| over workplaces",
           maximum(abs.(vec(sum(cb.pp_ijj, dims = 2)) .- 1)), 1e-13)
    report("C1c", "|sum_ij pi_ij - 1|", abs(sum(cb.pp_ij) - 1), 1e-13)

    # ---------------- C2: labour market clearing -------------------------------------
    println("\n-- C2: labour market clearing --")
    report("C2a", "max |sum_i pi_{j|i} H_Ri - H_Mj| / H_Mj",
           maximum(abs.((transpose(cb.pp_iji) * data.HR[Irsd]) ./ data.HM[Iwpl] .- 1)), 1e-8)
    report("C2b", "max |sum_j pi_{i|j} H_Mj - H_Ri| / H_Ri",
           maximum(abs.((transpose(cb.pp_ijj) * data.HM[Iwpl]) ./ data.HR[Irsd] .- 1)), 1e-8)
    report("C2c", "|sum_j H_Mj - sum_i H_Ri| / H",
           abs(sum(data.HM) - sum(data.HR)) / sum(data.HM), 1e-12)
    # smodexog.m uses vv = (pi_{j|i} w) .* H_R; smodendog.m uses vv = pi_{i|j}' (w .* H_M).
    # Algebraically the same; upstream never checks it, and the port relies on it.
    vv_a = (cb.pp_iji * sim.wage[Iwpl]) .* data.HR[Irsd]
    vv_b = transpose(cb.pp_ijj) * (sim.wage[Iwpl] .* data.HM[Iwpl])
    report("C2d", "expincome: the two equivalent vv formulas agree",
           maxlogdiff(vv_a, vv_b), 1e-12)

    # ---------------- C3: land market clearing, Eqs. 18 and 19 -----------------------
    println("\n-- C3: land market clearing (Eqs. 18 and 19) --")
    LM_eq = (((1 - p6.alpha) .* eq.A[Iwpl]) ./ eq.q[Iwpl]) .^ (1 / p6.alpha) .* eq.HM[Iwpl]
    report("C3a", "Eq. 18: ((1-a)A/q)^(1/a) H_M  vs  theta*L",
           maximum(abs.(LM_eq ./ (eq.theta[Iwpl] .* planted.LD[Iwpl]) .- 1)), 1e-8)
    LR_eq = ((1 - p6.beta) .* eq.vv[Irsd]) ./ eq.Q[Irsd]
    report("C3b", "Eq. 19: (1-b) E[w|i]H_R / Q  vs  (1-theta)*L",
           maximum(abs.(LR_eq ./ ((1 .- eq.theta[Irsd]) .* planted.LD[Irsd]) .- 1)), 1e-8)
    report("C3c", "Eq. S.31: LM + LR  vs  L = phi K^(1-mu)",
           maximum(abs.((LM_eq .+ LR_eq) ./ planted.LD .- 1)), 1e-8)

    # ---------------- C4: zero-shock counterfactual ----------------------------------
    println("\n-- C4: zero-shock counterfactual --")
    f0  = fund_from_quantification(data, sim)
    eq0 = solve_equilibrium_exog(p6, f0, data.tau; tol = 1e-13, maxiter = 200_000)
    eqz = solve_equilibrium_exog(p6, f0, copy(data.tau); tol = 1e-13, maxiter = 200_000)
    for (id, lab, x, y) in (("C4a", "H_M", eqz.HM, eq0.HM), ("C4b", "H_R", eqz.HR, eq0.HR),
                            ("C4c", "Q",   eqz.Q,  eq0.Q),  ("C4d", "theta", eqz.theta, eq0.theta))
        report(id, "zero-shock: max |change| in $lab", maximum(abs.(x .- y)), 0.0)
    end
    report("C4e", "zero-shock: |change| in Ubar", abs(eqz.Ubar - eq0.Ubar), 0.0)
    @printf("       baseline re-solve reproduces the data: H_M %.3e, H_R %.3e, Q %.3e (max log diff)\n",
            maxlogdiff(eq0.HM, data.HM), maxlogdiff(eq0.HR, data.HR), maxlogdiff(eq0.rent, data.Q))

    # ---------------- C5: the endogenous solver nests the exogenous one --------------
    println("\n-- C5: lambda = delta = eta = rho = 0  =>  smodendog == smodexog --")
    note("""Seeded at the observed data the two solvers agree trivially (both stop at
       iteration 1), so the comparison is run UNDER A SHOCK -- the wall going up -- where
       each solver has to travel several hundred iterations to a different allocation.""")
    p0 = ARSWParams(alpha = p6.alpha, beta = p6.beta, mu = p6.mu, epsilon = p6.epsilon,
                    kappaeps = p6.kappaeps)     # all four spillover parameters zero
    tau_wall = build_grid_geography(divided = true).tau
    eq_x = solve_equilibrium_exog(p6, f0, tau_wall; tol = 1e-13, maxiter = 200_000,
                                  damp = DAMP_ROBUST)
    eq_n = solve_equilibrium_endog(p0, f0, tau_wall; closed = true, tol = 1e-13,
                                   maxiter = 200_000, damp_far = DAMP_ROBUST,
                                   damp_near = DAMP_ROBUST)
    @printf("       under the wall shock: Alg. 11 took %d iterations, Alg. 15 took %d\n",
            eq_x.iters, eq_n.iters)
    report("C5a", "endog vs exog solver, under shock: wage", maxlogdiff(eq_n.wage, eq_x.wage), 1e-10)
    report("C5b", "endog vs exog solver, under shock: Q",    maxlogdiff(eq_n.Q, eq_x.Q), 1e-10)
    report("C5c", "endog vs exog solver, under shock: H_M",  maxlogdiff(eq_n.HM, eq_x.HM), 1e-10)
    report("C5d", "endog vs exog solver, under shock: theta",maxlogdiff(eq_n.theta, eq_x.theta), 1e-10)
    report("C5e", "endog vs exog solver, under shock: Ubar", abs(log(eq_n.Ubar / eq_x.Ubar)), 1e-10)

    # [D5] The same shock through upstream's smodexog, which freezes residential floor
    # prices at their OBSERVED level inside the commuting probability. If that were a
    # stylistic choice rather than a bug, C5 would still hold. It does not.
    eq_legacy = solve_equilibrium_exog(p6, f0, tau_wall; tol = 1e-13, maxiter = 200_000,
                                       damp = DAMP_ROBUST, legacy_QT = true)
    @printf("""
       [D5] the same shock with upstream's smodexog.m:152 (fixed QT in the commuting term):
            H_M differs by %.3e, H_R by %.3e, Q by %.3e, Ubar by %.3e (max log diff)
            -- i.e. the endogenous solver could NOT reproduce the exogenous one, by a wide
            margin, if the QT/Q_i substitution were left in place. It is a bug, not a choice.
""", maxlogdiff(eq_legacy.HM, eq_x.HM), maxlogdiff(eq_legacy.HR, eq_x.HR),
     maxlogdiff(eq_legacy.Q, eq_x.Q), abs(log(eq_legacy.Ubar / eq_x.Ubar)))

    # ---------------- C9 [D9]: smodexog's Q_e/q_e initialisation ----------------------
    println("\n-- C9 [D9]: smodexog.m's Q_e = q_e = QT initialisation vs smodendog.m's Q_e = Q_i --")
    ds = R.syn_s.data
    sim_s = quantify_simultaneous(p6, ds)
    fs = fund_from_quantification(ds, sim_s)
    nspec = count((sim_s.A .!= 0) .& (sim_s.B .== 0)) + count((sim_s.B .!= 0) .& (sim_s.A .== 0))
    @printf("       :specialised scenario, %d completely specialised blocks, upstream ROUNDING rule\n", nspec)
    e_ok  = solve_equilibrium_exog(p6, fs, ds.tau; rule = :round, maxiter = 1000,
                                   damp = DAMP_ROBUST, legacy_Qe_init = false)
    e_bad = solve_equilibrium_exog(p6, fs, ds.tau; rule = :round, maxiter = 1000,
                                   damp = DAMP_ROBUST, legacy_Qe_init = true)
    @printf("       smodendog.m init (Q_e = Q_i): converged=%-5s after %4d iterations\n",
            e_ok.converged, e_ok.iters)
    @printf("       smodexog.m  init (Q_e = QT ): converged=%-5s after %4d iterations\n",
            e_bad.converged, e_bad.iters)
    report("C9a", "both initialisations converge under the rounding rule",
           (e_ok.converged && e_bad.converged) ? 0.0 : 1.0, 0.0)
    report("C9b", "both reach the SAME allocation: H_M", maxlogdiff(e_bad.HM, e_ok.HM), 1e-10)
    report("C9c", "both reach the SAME allocation: Q",   maxlogdiff(e_bad.Q, e_ok.Q), 1e-10)
    report("C9d", "both reach the SAME allocation: theta", maxlogdiff(e_bad.theta, e_ok.theta), 1e-10)
    note("""Upstream's initialisation leaves Q_e[IcsA] = QT[IcsA] against Q_i[IcsA] = 0 and
       compares the FULL vectors in its stopping rule, which looks like it must fail forever
       on any completely specialised block. It does not: the damping blend drags those
       structurally unused entries onto their own initial values as well, so the rule is
       eventually satisfied. The cost is only the extra iterations above, plus meaningless
       nonzero values in the returned q[IcsB] / Q[IcsA] that `Crent` never reads. Recorded
       because the first guess was that it broke convergence outright -- it does not, and
       measuring beat asserting.""")

    # ---------------- C6: comparative statics -----------------------------------------
    println("\n-- C6: comparative statics, uniform +25% commuting time --")
    tau_up = data.tau .* 1.25
    eq_up  = solve_equilibrium_exog(p6, f0, tau_up; tol = 1e-13, maxiter = 200_000,
                                    damp = DAMP_ROBUST)
    cma_base = commuting_market_access(p6, data.tau[Irsd, Iwpl], eq0.wage[Iwpl])
    cma_up   = commuting_market_access(p6, tau_up[Irsd, Iwpl],  eq_up.wage[Iwpl])
    dcma = mean(cma_up ./ cma_base .- 1) * 100
    dw   = mean(eq_up.wage[Iwpl] ./ eq0.wage[Iwpl] .- 1) * 100
    dU   = (eq_up.Ubar / eq0.Ubar - 1) * 100
    dQ   = mean(eq_up.Q[Irsd] ./ eq0.Q[Irsd] .- 1) * 100
    @printf("       mean %% change: CMA %+.3f, wage %+.3f, Ubar %+.3f, Q %+.3f\n", dcma, dw, dU, dQ)
    report("C6a", "tau up => commuting market access falls (need < 0)", -dcma, 0.0;
           smaller_is_better = false)
    report("C6b", "tau up => wages fall (need < 0)", -dw, 0.0; smaller_is_better = false)
    report("C6c", "tau up => reservation utility falls (need < 0)", -dU, 0.0;
           smaller_is_better = false)
    note("""In a CLOSED city the wage response is second order: total employment and total
       floor space are both fixed, so mean w = alpha*Y/H_M barely moves and the shock
       shows up almost entirely in welfare. The first-order margin is population, which
       a closed city cannot use -- so the same shock is run OPEN, with Ubar pinned.""")
    p0o = ARSWParams(alpha = p6.alpha, beta = p6.beta, mu = p6.mu, epsilon = p6.epsilon,
                     kappaeps = p6.kappaeps)
    f0o = Fund(A = f0.A, B = f0.B, V = f0.V, K = f0.K, Q = f0.Q, HM = f0.HM, HR = f0.HR,
               LD = f0.LD, theta = f0.theta, wage = f0.wage, vv = f0.vv)
    eq_open_up = solve_equilibrium_endog(p0o, f0o, tau_up; closed = false,
                                         Utarget = eq0.Ubar, tol = 1e-11,
                                         maxiter = 200_000, damp_far = DAMP_ROBUST,
                                         damp_near = DAMP_ROBUST)
    dH = (eq_open_up.HH / sum(data.HM) - 1) * 100
    @printf("       open city, Ubar pinned: total employment H %+.3f%% (converged=%s, %d iters)\n",
            dH, eq_open_up.converged, eq_open_up.iters)
    report("C6d", "tau up, open city => total employment falls (need < 0)", -dH, 0.0;
           smaller_is_better = false)

    # ---------------- C7: the two scale invariances -----------------------------------
    println("\n-- C7: scale invariances the normalisations rest on --")
    kA = 3.7
    fA = Fund(A = f0.A .* kA, B = f0.B, V = f0.V, K = f0.K, Q = f0.Q .* kA,
              HM = f0.HM, HR = f0.HR, LD = f0.LD, theta = f0.theta,
              wage = f0.wage .* kA, vv = f0.vv .* kA)
    eqA = solve_equilibrium_exog(p6, fA, data.tau; tol = 1e-13, maxiter = 200_000)
    report("C7a", "a -> k*a leaves H_M unchanged",  maxlogdiff(eqA.HM, eq0.HM), 1e-10)
    report("C7b", "a -> k*a leaves theta unchanged", maxlogdiff(eqA.theta, eq0.theta), 1e-10)
    report("C7c", "a -> k*a scales Q by exactly k",  maxlogdiff(eqA.Q ./ kA, eq0.Q), 1e-10)
    cB = 0.41
    fB = Fund(A = f0.A, B = f0.B .* cB, V = f0.V, K = f0.K, Q = f0.Q, HM = f0.HM,
              HR = f0.HR, LD = f0.LD, theta = f0.theta, wage = f0.wage, vv = f0.vv)
    eqB = solve_equilibrium_exog(p6, fB, data.tau; tol = 1e-13, maxiter = 200_000)
    report("C7d", "b -> c*b leaves H_R unchanged", maxlogdiff(eqB.HR, eq0.HR), 1e-10)
    report("C7e", "b -> c*b leaves Q unchanged",   maxlogdiff(eqB.Q, eq0.Q), 1e-10)
    report("C7f", "b -> c*b scales Phi by c^eps",
           abs(log(eqB.Phi / (eq0.Phi * cB ^ p6.epsilon))), 1e-10)

    # ---------------- C8: mu is inert --------------------------------------------------
    println("\n-- C8 [D3]: mu only relabels phi; no equilibrium object moves --")
    p_mu = ARSWParams(alpha = p6.alpha, beta = p6.beta, mu = 0.25,   # upstream's effective value
                      epsilon = p6.epsilon, kappaeps = p6.kappaeps)
    sim_mu = quantify_simultaneous(p_mu, data)
    report("C8a", "mu: A~ unchanged",  maxlogdiff(sim_mu.A, sim.A), 1e-12)
    report("C8b", "mu: B~ unchanged",  maxlogdiff(sim_mu.B, sim.B), 1e-12)
    report("C8c", "mu: total floor space L unchanged", maxlogdiff(sim_mu.LD, sim.LD), 1e-12)
    report("C8d", "mu: theta unchanged", maxlogdiff(sim_mu.theta, sim.theta), 1e-12)
    report_note("C8*", "mu: phi DOES change (that is the point)",
                maxlogdiff(sim_mu.V, sim.V))
    f_mu  = fund_from_quantification(data, sim_mu)
    eq_mu = solve_equilibrium_exog(p_mu, f_mu, data.tau; tol = 1e-13, maxiter = 200_000)
    report("C8e", "mu: solved equilibrium Q unchanged", maxlogdiff(eq_mu.Q, eq0.Q), 1e-12)

    return (eq0 = eq0, f0 = f0, eq_up = eq_up, cma_base = cma_base)
end

"""
    run_rounding_comparison(p6, R)

How much the departure from upstream's rounding stopping rule buys. Reported because the
plan committed to justifying that departure numerically rather than by assertion.
"""
function run_rounding_comparison(p6::ARSWParams, R)
    header("STOPPING RULE: upstream's rounding test vs a tolerance test")
    data, planted = R.syn.data, R.syn.planted
    for (lab, rule) in (("upstream round(gap*1e4)==0", :round), ("tolerance 1e-14", :tol))
        sq = quantify_sequential(p6, data; rule = rule)
        sm = quantify_simultaneous(p6, data; rule = rule)
        @printf("  %-30s  seq: A~ %.3e  B~ %.3e (%d it)  |  sim: A~ %.3e  B~ %.3e (%d it)\n",
                lab, maxlogdiff(sq.A, planted.A), maxlogdiff(sq.B, planted.B), sq.iters,
                maxlogdiff(sm.A, planted.A), maxlogdiff(sm.B, planted.B), sm.iters)
    end
    note("""
      The rounding rule stops at a fixed ABSOLUTE precision in levels, so how good it is
      depends on the units of the data. It is fine for mapping and for the paper's
      purposes, but it is not good enough to distinguish a correct inversion from a subtly
      wrong one, which is why the recovery tests above use the tolerance rule.""")
end

"""
    run_wall_conditioning(p6)

How the wall degrades the CONDITIONING (not the correctness) of the inversion.

Commuting market clearing, Eq. (S.44), determines transformed wages only up to one
global scale, and the uniqueness proof (Lemmas S.6-S.7) needs the commuting matrix to be
irreducible. A hard wall makes it nearly reducible, so the two halves' wage levels are
tied together only through flows of order `exp(-nu*penalty)`: the slow mode of the damped
iteration has modulus `1 - O(exp(-nu*penalty))` and the iteration count scales like
`exp(+nu*penalty)`. The table below shows that scaling directly. Every row still
converges to machine precision -- which is the point: this is a numerical-conditioning
property of a severed city, not a defect, and it is why ARSW estimate epsilon on the
WEST-ONLY 1986 sample instead of on the divided city as a whole.
"""
function run_wall_conditioning(p6::ARSWParams)
    header("DIAGNOSTIC: what a wall does to the conditioning of the inversion")
    @printf("  %8s %14s %12s %12s %14s\n",
            "penalty", "exp(-nu*pen)", "Alg.4 iters", "Alg.9 iters", "Alg.9 A~ error")
    for pen in (0.0, 15.0, 30.0, 45.0)
        gd = build_grid_geography(divided = pen > 0, wall_penalty = max(pen, 1.0))
        sd = make_synthetic_data(; scenario = :mixed, p = p6, geo = gd)
        w  = solve_wages_and_productivity(p6, sd.data; tol = 1e-14, maxiter = 200_000)
        r  = invert_simultaneous(p6, sd.data; tol = 1e-14, maxiter = 200_000)
        @printf("  %8.0f %14.3e %12d %12d %14.3e\n", pen, exp(-p6.kappaeps * pen),
                w.iters, r.iters, maxlogdiff(r.A, sd.planted.A))
    end
    note("""
      Iteration count tracks 1/exp(-nu*penalty) almost exactly; the recovery error does
      not move. Correctness is unaffected -- only the number of iterations needed.""")
end

"""
    run_epsilon_recovery(p6, R)

STRETCH: Algorithms 1-3. A genuine recovery test for epsilon.

The identification (paper p. 2167): commuting market clearing, Eq. (S.44), determines
transformed wages `omega_j` from `{H_M, H_R, tau}` ALONE -- epsilon does not appear.
Epsilon only fixes the monotonic map `w~_j = omega_j^(1/eps)`, hence
`var(log w~) = (1/eps)^2 var(log omega)`, so one moment on the dispersion of observed
wages pins it down.

The test: Algorithm 1 sees only employment and travel times (no wage data at all); the
single data moment is the variance of log DISTRICT wages computed from the TRUE
equilibrium wages. If the port is right, the search must return the planted epsilon.

Upstream's moment is the variance of log wages across the 12 West-Berlin Bezirke in 1986
(`wageworker1986.csv` + `modbezirk.m`). The synthetic stand-in is a 3x3 partition of the
grid into 5x5 districts. The reunified geography is used rather than the divided one:
with the wall up the commuting matrix is nearly block diagonal, so the relative wage
SCALE of the two halves is pinned only through terms of order `exp(-8.4)` and the omega
fixed point becomes numerically ill-conditioned -- an artefact of the synthetic wall
being a pure travel-time penalty, not a statement about the real 1986 data (where West
Berlin is a self-contained city, which is exactly why upstream estimates on the West-only
sample).
"""
function run_epsilon_recovery(p6::ARSWParams, R)
    header("STRETCH -- Algorithms 1-3: recover epsilon from synthetic data")
    data, geo, eq = R.syn.data, R.syn.geo, R.syn.eq

    # 3 x 3 districts of 5 x 5 blocks -- the synthetic stand-in for modbezirk.m.
    dsize = 5
    district = [(min(2, (geo.col[i] - 1) ÷ dsize) + 1) +
                3 * min(2, (geo.row[i] - 1) ÷ dsize) for i in 1:geo.N]

    # The one data moment: dispersion of log district wages at the TRUE wages.
    wtrue = copy(eq.wage);  wtrue ./= geomean_pos(wtrue)
    wbill = wtrue .* data.HM
    ids = sort(unique(district))
    dw  = [sum(wbill[district .== d]) / sum(data.HM[district .== d]) for d in ids]
    lbw = log.(dw);  lbw .-= mean(lbw)
    var_data = var(lbw)
    @printf("  Sample: %d blocks in %d districts; var(log district wage) = %.8e\n",
            geo.N, length(dw), var_data)
    @printf("  Planted epsilon = %.6f   (Algorithm 1 never sees a wage)\n", p6.epsilon)

    est = estimate_epsilon(p6, data, district, var_data)

    # Sanity: the recovered omega must be the true wages^eps up to scale.
    what = est.omega .^ (1 / p6.epsilon);  what ./= geomean_pos(what)
    @printf("    omega check: max |log(w_hat / w_true)| = %.3e\n", maxlogdiff(what, wtrue))
    @printf("  Brent      estimate: %.10f   (relative error %.3e)\n",
            est.epsilon, abs(est.epsilon / p6.epsilon - 1))
    @printf("  NelderMead estimate: %.10f   (relative error %.3e)\n",
            est.epsilon_neldermead, abs(est.epsilon_neldermead / p6.epsilon - 1))
    @printf("  Objective at the optimum: %.6e ; upstream rounds to 2dp -> %.2f\n",
            est.objective, round(est.epsilon * 100) / 100)
    report("S1a", "epsilon recovery (Brent), relative error",
           abs(est.epsilon / p6.epsilon - 1), 1e-6)
    report("S1b", "epsilon recovery (NelderMead), relative error",
           abs(est.epsilon_neldermead / p6.epsilon - 1), 1e-5)
    report("S1c", "Algorithm 1 recovers true wages up to scale",
           maxlogdiff(what, wtrue), 1e-10)
    return est
end

# =====================================================================================
# 17. COUNTERFACTUALS
# =====================================================================================

"""
    run_counterfactuals(p6, p7, R, C)

Three exercises, all seeded from the simultaneous quantification (as upstream's
`cftualexog_TD.m` / `cftualendog_*_TD.m` do):

1. DIVISION: the wall goes up. Same fundamentals, travel times across the cut raised by
   the wall penalty. The synthetic analogue of 1961.
2. REUNIFICATION: quantify the model on DIVIDED data and take the wall away, the
   analogue of 1989. Because the fundamentals are inverted from divided-city data, this
   is the direction the paper actually exploits.
3. EAST PRODUCTIVITY +10%: upstream's `cftualexog_TD.m` "eastern renewal" counterfactual,
   run both with exogenous fundamentals (Alg. 11) and with endogenous agglomeration
   (Alg. 15), so the amplification from the spillovers is visible.
"""
function run_counterfactuals(p6::ARSWParams, p7::ARSWParams, R, C)
    header("COUNTERFACTUALS")
    data, geo = R.syn.data, R.syn.geo
    f0, eq0 = C.f0, C.eq0
    west = data.west

    pct(x, y) = (x ./ y .- 1) .* 100

    # ---- 1. Division ------------------------------------------------------------------
    println("\n-- 1. DIVISION: the wall goes up (same fundamentals, +45 min across the cut) --")
    geo_div = build_grid_geography(divided = true)
    # Upstream's fixed damping weight of 0.5 (smodexog.m:286) does not converge on this
    # shock; it settles into a period-2 limit cycle. Shown, then done properly. See the
    # docstring of DAMP_ROBUST.
    eq_bad = solve_equilibrium_exog(p6, f0, geo_div.tau; tol = 1e-13, maxiter = 3_000,
                                    damp = 0.5)
    @printf("   upstream damping 0.5 : converged=%s after %d iterations, max log gap stuck at %.3e\n",
            eq_bad.converged, eq_bad.iters, eq_bad.path[end, 5])
    eq_div = solve_equilibrium_exog(p6, f0, geo_div.tau; tol = 1e-13, maxiter = 200_000,
                                    damp = DAMP_ROBUST)
    @printf("   damping %.2f        : converged=%s in %d iterations (max log gap %.3e)\n",
            DAMP_ROBUST, eq_div.converged, eq_div.iters, eq_div.path[end, 5])
    @printf("   Ubar  %+.3f%%\n", (eq_div.Ubar / eq0.Ubar - 1) * 100)
    for (lab, x, y) in (("workplace employment H_M", eq_div.HM, eq0.HM),
                        ("residence employment H_R", eq_div.HR, eq0.HR),
                        ("floor price Q",            eq_div.Q,  eq0.Q),
                        ("wage w",                   eq_div.wage, eq0.wage))
        ch = pct(x, y)
        @printf("   %-26s  mean %+7.3f%%   CBD %+7.3f%%   edge %+7.3f%%\n", lab,
                mean(filter(isfinite, ch)),
                mean(ch[geo.dist_cbd .< 3]), mean(ch[geo.dist_cbd .> 7]))
    end
    # The paper's signature: the gradient breaks at the wall.
    near_wall = abs.(geo.coords[:, 1] .- 6.5) .< 1.6
    @printf("   blocks within 1.6 km of the wall: H_R %+.3f%%, Q %+.3f%% (city-wide Q %+.3f%%)\n",
            mean(pct(eq_div.HR, eq0.HR)[near_wall]),
            mean(pct(eq_div.Q, eq0.Q)[near_wall]), mean(pct(eq_div.Q, eq0.Q)))

    # ---- 2. Reunification -------------------------------------------------------------
    println("\n-- 2. REUNIFICATION: quantify on DIVIDED data, then remove the wall --")
    syn_div = make_synthetic_data(; scenario = :mixed, p = p6, geo = geo_div)
    d_div   = syn_div.data
    sim_div = quantify_simultaneous(p6, d_div)
    @printf("   inversion on divided data: A~ %.3e, B~ %.3e (recovery vs planted)\n",
            maxlogdiff(sim_div.A, syn_div.planted.A), maxlogdiff(sim_div.B, syn_div.planted.B))
    f_div    = fund_from_quantification(d_div, sim_div)
    eq_divb  = solve_equilibrium_exog(p6, f_div, geo_div.tau; tol = 1e-13,
                                      maxiter = 200_000, damp = DAMP_ROBUST)
    eq_reun  = solve_equilibrium_exog(p6, f_div, geo.tau; tol = 1e-13,
                                      maxiter = 200_000, damp = DAMP_ROBUST)
    @printf("   reunification: Ubar %+.3f%%\n", (eq_reun.Ubar / eq_divb.Ubar - 1) * 100)
    chHR = pct(eq_reun.HR, eq_divb.HR);  chQ = pct(eq_reun.Q, eq_divb.Q)
    @printf("   H_R  near wall %+7.3f%%   far from wall %+7.3f%%\n",
            mean(chHR[near_wall]), mean(chHR[.!near_wall]))
    @printf("   Q    near wall %+7.3f%%   far from wall %+7.3f%%\n",
            mean(chQ[near_wall]), mean(chQ[.!near_wall]))
    @printf("   (the paper's reduced form: reunification lifts land prices and density\n")
    @printf("    disproportionately in blocks close to the former wall -- reproduced in sign here)\n")

    # ---- 3. East productivity boost ---------------------------------------------------
    println("\n-- 3. EASTERN RENEWAL: +10% productivity in the East (cftualexog_TD.m) --")
    east = .!west
    f_ae = Fund(A = f0.A .* (1 .+ 0.1 .* east), B = f0.B, V = f0.V, K = f0.K, Q = f0.Q,
                HM = f0.HM, HR = f0.HR, LD = f0.LD, theta = f0.theta,
                wage = f0.wage, vv = f0.vv)
    eq_ae = solve_equilibrium_exog(p6, f_ae, data.tau; tol = 1e-13, maxiter = 200_000,
                                   damp = DAMP_ROBUST)
    @printf("   exogenous fundamentals (Alg. 11): Ubar %+.3f%%\n", (eq_ae.Ubar/eq0.Ubar - 1)*100)
    @printf("      H_M  East %+7.3f%%   West %+7.3f%%\n",
            mean(pct(eq_ae.HM, eq0.HM)[east]), mean(pct(eq_ae.HM, eq0.HM)[west]))
    @printf("      Q    East %+7.3f%%   West %+7.3f%%\n",
            mean(pct(eq_ae.Q, eq0.Q)[east]), mean(pct(eq_ae.Q, eq0.Q)[west]))

    d7, pl7 = R.syn7.data, R.syn7.planted
    sim7 = R.sim7;  dp = R.dp;  dr = R.dr
    f7 = Fund(A = dp.a, B = dr.b, V = sim7.V, K = copy(d7.K), Q = copy(d7.Q),
              HM = copy(d7.HM), HR = copy(d7.HR), LD = sim7.LD, theta = sim7.theta,
              wage = sim7.wage, vv = sim7.vv)
    eq7b = solve_equilibrium_endog(p7, f7, d7.tau; closed = true, tol = 1e-13,
                                   maxiter = 200_000, damp_far = DAMP_ROBUST,
                                   damp_near = DAMP_ROBUST)
    f7a  = Fund(A = dp.a .* (1 .+ 0.1 .* east), B = dr.b, V = f7.V, K = f7.K, Q = f7.Q,
                HM = f7.HM, HR = f7.HR, LD = f7.LD, theta = f7.theta,
                wage = f7.wage, vv = f7.vv)
    eq7a = solve_equilibrium_endog(p7, f7a, d7.tau; closed = true, tol = 1e-13,
                                   maxiter = 200_000, damp_far = DAMP_ROBUST,
                                   damp_near = DAMP_ROBUST)
    @printf("   endogenous agglomeration (Alg. 15): Ubar %+.3f%%\n",
            (eq7a.Ubar / eq7b.Ubar - 1) * 100)
    @printf("      H_M  East %+7.3f%%   West %+7.3f%%\n",
            mean(pct(eq7a.HM, eq7b.HM)[east]), mean(pct(eq7a.HM, eq7b.HM)[west]))
    @printf("      total A East %+7.3f%% (the shock is +10%% to the FUNDAMENTAL a~;\n",
            mean(pct(eq7a.A, eq7b.A)[east]))
    @printf("       the gap to +10%% is the productivity externality, Eq. 20, at work)\n")

    return (eq_div = eq_div, eq_reun = eq_reun, eq_divb = eq_divb, eq_ae = eq_ae,
            geo_div = geo_div, syn_div = syn_div, near_wall = near_wall,
            eq7a = eq7a, eq7b = eq7b)
end

# =====================================================================================
# 18. FIGURES
# =====================================================================================

"Reshape a per-block vector onto the (ny, nx) grid for heatmapping."
grid_of(geo, v) = reshape(v, geo.ny, geo.nx)

function make_figures(p6, R, C, CF)
    MAKE_PLOTS || (println("\n[plots disabled: ARSW_PLOTS=0]"); return)
    header("FIGURES")
    gdir = joinpath(@__DIR__, "graphs")          # [D2] self-locating, no hard-coded root
    isdir(gdir) || mkpath(gdir)
    geo, data, planted = R.syn.geo, R.syn.data, R.syn.planted
    sim, eq0 = R.sim, C.eq0

    # 1. Planted vs recovered fundamentals, side by side on the grid.
    f1 = plot(
        heatmap(grid_of(geo, log.(planted.A)), title = "planted  log A~", c = :viridis),
        heatmap(grid_of(geo, log.(sim.A)),     title = "recovered  log A~", c = :viridis),
        heatmap(grid_of(geo, log.(planted.B)), title = "planted  log B~", c = :plasma),
        heatmap(grid_of(geo, log.(sim.B)),     title = "recovered  log B~", c = :plasma),
        layout = (2, 2), size = (900, 720), aspect_ratio = :equal,
        plot_title = "Recovery test: planted vs inverted fundamentals (synthetic)")
    savefig(f1, joinpath(gdir, "recovery_fundamentals.pdf"))

    # 2. Recovery scatter: the whole point, on one axis.
    f2 = scatter(log.(planted.A), log.(sim.A), label = "log A~", ms = 3, mc = :steelblue,
                 xlabel = "planted", ylabel = "recovered (Algorithm 9)",
                 title = "Inversion recovers the planted fundamentals", legend = :topleft)
    scatter!(f2, log.(planted.B), log.(sim.B), label = "log B~", ms = 3, mc = :darkorange)
    lims = extrema(vcat(log.(planted.A), log.(planted.B)))
    plot!(f2, collect(lims), collect(lims), lc = :black, ls = :dash, label = "45 degrees")
    savefig(f2, joinpath(gdir, "recovery_scatter.pdf"))

    # 3. Convergence paths of the equilibrium solver.
    pth = eq0.path
    f3 = plot(pth[:, 6], pth[:, 1], lw = 2, label = "wages", yscale = :log10,
              xlabel = "iteration", ylabel = "max log gap",
              title = "Algorithm 11 convergence (smodexog)")
    plot!(f3, pth[:, 6], pth[:, 2], lw = 2, ls = :dash, label = "floor price q")
    plot!(f3, pth[:, 6], pth[:, 4], lw = 2, ls = :dot,  label = "land use theta")
    savefig(f3, joinpath(gdir, "convergence_smodexog.pdf"))

    # 4. The wall: percentage changes from division.
    chHR = (CF.eq_div.HR ./ eq0.HR .- 1) .* 100
    chQ  = (CF.eq_div.Q  ./ eq0.Q  .- 1) .* 100
    f4 = plot(
        heatmap(grid_of(geo, chHR), title = "division: % change H_R", c = :balance,
                clims = (-maximum(abs.(chHR)), maximum(abs.(chHR)))),
        heatmap(grid_of(geo, chQ),  title = "division: % change Q",  c = :balance,
                clims = (-maximum(abs.(chQ)), maximum(abs.(chQ)))),
        layout = (1, 2), size = (1000, 420), aspect_ratio = :equal,
        plot_title = "The synthetic wall (severs commuting between columns 7 and 8)")
    savefig(f4, joinpath(gdir, "wall_division.pdf"))

    # 5. Gradients from the CBD, and the wall's gradient break.
    xs = geo.coords[:, 1]
    f5 = scatter(xs, chQ, ms = 3, mc = :firebrick, label = "% change in Q",
                 xlabel = "east-west position (km)", ylabel = "% change from division",
                 title = "Floor-price gradient breaks at the wall")
    vline!(f5, [6.5], lc = :black, ls = :dash, label = "wall")
    savefig(f5, joinpath(gdir, "wall_gradient_break.pdf"))

    # 6. Section 7 decomposition: total vs exogenous productivity and amenity.
    geo7 = R.syn7.geo
    f6 = plot(
        scatter(geo7.dist_cbd, log.(R.sim7.A), ms = 3, label = "total log A~",
                xlabel = "distance from CBD (km)", title = "Productivity"),
        scatter(geo7.dist_cbd, log.(R.dp.a), ms = 3, label = "exogenous log a~", mc = :orange),
        scatter(geo7.dist_cbd, log.(R.sim7.B), ms = 3, label = "total log B~",
                xlabel = "distance from CBD (km)", title = "Amenity", mc = :steelblue),
        scatter(geo7.dist_cbd, log.(R.dr.b), ms = 3, label = "exogenous log b~", mc = :purple),
        layout = (2, 2), size = (950, 700),
        plot_title = "Section 7: externalities vs fundamentals (Eqs. 20-21)")
    savefig(f6, joinpath(gdir, "section7_decomposition.pdf"))

    # 7. Commuting market access over the grid, divided vs reunified.
    Iw = data.HM .!= 0;  Ir = data.HR .!= 0
    cma_r = zeros(geo.N); cma_r[Ir] = commuting_market_access(p6, data.tau[Ir, Iw], eq0.wage[Iw])
    cma_d = zeros(geo.N); cma_d[Ir] = commuting_market_access(p6, CF.geo_div.tau[Ir, Iw],
                                                              CF.eq_div.wage[Iw])
    f7 = plot(
        heatmap(grid_of(geo, log.(cma_r)), title = "log CMA, reunified", c = :viridis),
        heatmap(grid_of(geo, log.(cma_d)), title = "log CMA, divided",  c = :viridis),
        layout = (1, 2), size = (1000, 420), aspect_ratio = :equal,
        plot_title = "Commuting market access, Eq. (S.46)")
    savefig(f7, joinpath(gdir, "commuting_market_access.pdf"))

    println("  7 figures written to $gdir")
end

# =====================================================================================
# 19. MAIN
# =====================================================================================

function main()
    t0 = time()
    header("ARSW (2015), The Economics of Density -- Julia port on SYNTHETIC data")
    println("""
  Ahlfeldt, Redding, Sturm & Wolf (2015), Econometrica 83(6), 2127-2189.
  Ported from Gabriel M. Ahlfeldt's MATLAB toolkit (16 codebook algorithms).

  THIS RUNS ON SYNTHETIC DATA. The three .mat files the toolkit needs are external
  downloads and are not present; they are deliberately not fetched. No empirical result
  of the paper is reproduced or claimed. The paper's parameter values are INPUTS here.

  The validation is a recovery test: plant fundamentals, solve the model forward to
  manufacture 'observed' data, invert, and check the planted values come back. Supplement
  Propositions S.3/S.4 say they must.""")

    p6 = PARAMS_SECTION6
    p7 = PARAMS_SECTION7
    @printf("\n  Section 6 parameters: alpha=%.2f beta=%.2f mu=%.2f eps=%.3f nu=%.4f kappa=%.6f\n",
            p6.alpha, p6.beta, p6.mu, p6.epsilon, p6.kappaeps, kappa_of(p6))
    @printf("  Section 7 parameters: lambda=%.4f delta=%.4f eta=%.4f rho=%.4f eps=%.3f nu=%.4f\n",
            p7.lambda, p7.delta, p7.eta, p7.rho, p7.epsilon, p7.kappaeps)
    println("""
  REMINDER: lambda here is the AGGLOMERATION elasticity (higher lambda RAISES
  productivity). In AllenArkolakis-RES-2022 in this same repository, lambda is the traffic
  CONGESTION elasticity at almost the same magnitude, with the opposite sign of effect.
  Do not transplant a calibrated value between the two.""")

    R  = run_recovery_tests(p6, p7)
    C  = run_model_checks(p6, p7, R)
    run_rounding_comparison(p6, R)
    run_wall_conditioning(p6)
    est = run_epsilon_recovery(p6, R)
    CF = run_counterfactuals(p6, p7, R, C)
    make_figures(p6, R, C, CF)

    header("SUMMARY")
    npass = count(r -> r.pass, RESULTS)
    for r in RESULTS
        r.pass || @printf("  FAILED  %-5s %-58s %.4e (target %.0e)\n",
                          r.id, r.label, r.value, r.target)
    end
    @printf("  %d of %d checks passed.\n", npass, length(RESULTS))
    @printf("  worst relative error across all recovery tests: %.4e\n",
            maximum(r.value for r in RESULTS if startswith(r.id, "R")))
    @printf("  elapsed: %.1f s\n", time() - t0)
    return npass == length(RESULTS)
end

if abspath(PROGRAM_FILE) == @__FILE__
    ok = main()
    exit(ok ? 0 : 1)
end
