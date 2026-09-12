# ============================================================================
# Endogenous Freight Costs on a Transport Network - NO PRICED PASSAGES
#
# A variant of Brockhaus, Hinz & Serfaty (2026), "Navigating Shocks", that keeps
# the Caliendo-Parro system and the nested mode-over-route transport block but
# removes the monopoly toll layer entirely. Every passage is free; the only cost
# of traversing one is the physical cost of the link, which is itself endogenous
# through congestion and the global shipping capacity constraint.
#
# WHAT IS KEPT (identical to the full model)
#   - Multi-country, multi-sector Caliendo-Parro with input-output linkages
#   - Freight cost share pass-through   d = f^chi_s                      (eq 5)
#   - Route congestion   Xi_mr = sum_ijs delta * X_ij,smr            (eq 6, 19)
#   - Global shipping capacity   Psi = sum delta * X                 (eq 7, 20)
#   - Route choice CES within a mode, elasticity eta^R              (eq 14, 15)
#   - Mode choice CES across modes, elasticity eta^M_s              (eq 16, 17)
#   - Route-level costs                                                 (eq 18)
#   - Solution in changes (exact hat algebra, Dekle-Eaton-Kortum 2008)
#
# WHAT IS REMOVED
#   - Monopoly toll setting  phi*_q = argmax Pi_q                       (eq 12)
#   - Toll revenue and its log-share division on shared routes          (eq 13)
#   - The toll authorities, the rent rebate to collector countries, and
#     the re-baselining of the equilibrium at calibrated fees
#   - Consequently: the welfare decomposition into rent and residual channels,
#     and the with-fees / without-fees comparison of Section 6.3
#
# This is, in effect, the paper's own "no chokepoint rents" world (Section 6.3)
# promoted from a robustness exercise to the whole model. Incidence therefore
# falls on the economies whose cargo actually has to travel farther, rather than
# on the collector - which is the paper's point, run as a standalone framework.
#
# ---------------------------------------------------------------------------
# THE COST OF REMOVING THE TOLLS: eta^R IS NO LONGER IDENTIFIED
# ---------------------------------------------------------------------------
# In the full model the route-choice elasticity eta^R (the paper's mu^d) is
# pinned by a revenue moment: the revenue-maximising Suez fee must reproduce the
# Canal Authority's observed toll income. With no authority and no fee, that
# moment does not exist, and nothing else in the model identifies eta^R.
#
# This matters, because eta^R is exactly the parameter the headline number is
# proportional to. Closing a route raises the within-mode composite by
#
#       d_hat  =  (1 - rho_closed)^(-1/eta^R)
#
# so the delivered-cost shock - and with it the welfare loss - scales as
# 1/eta^R. Section 6 below reports that sensitivity rather than burying it: the
# model can tell you the SHAPE and INCIDENCE of a chokepoint disruption without
# the toll moment, but not its magnitude.
#
# Candidate external moments, none of them clean:
#   - the paper's own difference-in-differences trough (about -10% for a fully
#     seaborne good facing a 9,000 km detour) would imply eta^R near 95, but it
#     is a SHORT-RUN estimate and the paper documents that trade recovered
#     within a quarter, so it is the wrong moment for a long-run static model;
#   - observed route shares pin only the RATIO of baseline costs, not eta^R;
#   - freight-rate pass-through during the crisis confounds eta^R with the
#     congestion and capacity elasticities lambda and gamma^G.
# The default below is the paper's calibrated 350, used as a reference value and
# NOT as an identified estimate.
#
# ---------------------------------------------------------------------------
# SCOPE
# ---------------------------------------------------------------------------
# Same stylized 24-country x 8-sector baseline as the full model: the paper's
# real chokepoint geography and Gulf-side exposure shares g_i, its calibrated
# elasticities, and a synthetic but internally consistent input-output core.
# It reproduces mechanisms, not GTAP magnitudes.
#
# Because no passage is priced, the policy experiments available here are
# CLOSURES and EXOGENOUS COST SHOCKS (war-risk premia, drought restrictions),
# not fees. A fee would be a toll, which is precisely what this variant drops.
#
# Run with:  julia brockhaus_hinz_serfaty_no_tolls_model.jl
#            BHS_QUICK=1 julia ...   # baseline + verification only
# ============================================================================

using LinearAlgebra, Statistics, Random, Printf, Plots

Random.seed!(1)
gr()
default(titlefontsize = 10, guidefontsize = 9, tickfontsize = 8, legendfontsize = 8,
        left_margin = 6Plots.mm, bottom_margin = 5Plots.mm, top_margin = 4Plots.mm,
        right_margin = 4Plots.mm, size = (640, 430))

const GRAPHDIR = joinpath(@__DIR__, "graphs_no_tolls")
isdir(GRAPHDIR) || mkpath(GRAPHDIR)

banner(s) = (println(); println("="^78); println(s); println("="^78); flush(stdout))
say(args...) = (println(args...); flush(stdout))

# ============================================================================
# 1. PARAMETERS AND ELASTICITIES  (Section 5.2, Table 4)
# ============================================================================

const NS = 8    # sectors
const NM = 3    # modes: 1 = sea, 2 = air, 3 = other (land / pipeline / services)
const NR = 8    # routes

const SECTORS = ["Agri & food", "Energy", "Minerals & bulk", "Chemicals",
                 "Metals", "Machinery & veh.", "Electronics & pharma", "Services"]

# Freight cost share chi_s (Section 5.3): share of DELIVERED price that is
# route-level shipping cost. Bulk 0.20-0.40; steel and basic chemicals
# 0.05-0.10; manufactures 0.02-0.04; electronics and pharma 0.01. Crude and LNG
# sit near 0.10 (see the note in the full model on the paper's implied value).
const CHI = [0.20, 0.10, 0.30, 0.08, 0.08, 0.03, 0.01, 0.00]

# Trade elasticity theta_s. PLACEHOLDERS: the paper uses Fontagne et al. (2022),
# theta_s = sigma_s - 1, which the paper does not tabulate. The large energy
# value reflects the paper's caveat about high long-run substitution.
const THETA = [8.0, 20.0, 8.0, 6.0, 7.0, 5.0, 5.0, 4.0]

# Mode substitution elasticity eta^M_s (Ko et al. 2025, paper Table 4);
# pooled fallback 2.44 (Tolva 2026) for sectors without an estimate.
const ETA_M = [4.87, 11.08, 2.68, 3.13, 5.21, 3.17, 1.67, 2.44]

# Transport-block elasticities (Section 5.2)
const LAMBDA  = fill(0.135, NS)   # route congestion: alpha_tau * beta_tau
const GAMMA_G = fill(0.44, NS)    # global shipping capacity (impact elasticity)
const GAMMA_M = fill(0.00, NS)    # mode-specific capacity: off, as in the paper

# Route-choice elasticity. NOT IDENTIFIED IN THIS VARIANT - see the header.
# 350 is the value the full model calibrates from the Suez revenue moment; here
# it is a reference value only, and Section 6 sweeps it.
const ETA_R_REF = 350.0

# Ton-km capacity weights per dollar shipped (units arbitrary; only ratios matter)
const TPD = [2.0, 5.0, 8.0, 1.0, 3.0, 0.3, 0.05, 0.0]

# Production structure
const GAMMA_VA  = [0.40, 0.50, 0.45, 0.30, 0.30, 0.32, 0.35, 0.60]
const ALPHA_FIN = [0.12, 0.05, 0.03, 0.04, 0.03, 0.10, 0.08, 0.55]

# Mode shares by sector for international trade (sea, air, other).
# Domestic trade is 100% "other", so it is never exposed to a maritime chokepoint.
const MODESHARE = [0.85 0.02 0.13;   # agri & food
                   0.92 0.00 0.08;   # energy
                   0.90 0.00 0.10;   # minerals & bulk
                   0.80 0.05 0.15;   # chemicals
                   0.85 0.02 0.13;   # metals
                   0.75 0.10 0.15;   # machinery & vehicles
                   0.45 0.40 0.15;   # electronics & pharma
                   0.00 0.00 1.00]   # services

# ============================================================================
# 2. GEOGRAPHY: COUNTRIES AND PASSAGES  (Section 5.1)
# ============================================================================
# Note what is NOT here relative to the full model: no AUTH_NAME, no
# AUTH_ROUTES, no AUTH_COLLECTORS, no ROUTE_AUTHS. Passages are geography only.

const COUNTRY_DATA = [
    ("USA", "United States",       19485.0,  39.8,  -98.6, :AMERICAS, 0.00),
    ("CHN", "China",               12310.0,  35.0,  105.0, :EASIA,    0.00),
    ("JPN", "Japan",                4872.0,  36.0,  138.0, :EASIA,    0.00),
    ("DEU", "Germany",              3693.0,  51.0,   10.0, :EUR,      0.00),
    ("IND", "India",                2652.0,  21.0,   78.0, :SASIA,    0.00),
    ("FRA", "France",               2595.0,  46.0,    2.0, :EUR,      0.00),
    ("ITA", "Italy",                1950.0,  42.0,   12.0, :EUR,      0.00),
    ("KOR", "Korea",                1623.0,  36.0,  128.0, :EASIA,    0.00),
    ("NLD", "Netherlands",           831.0,  52.0,    5.5, :EUR,      0.00),
    ("SAU", "Saudi Arabia",          689.0,  24.0,   45.0, :GULF,     0.80),
    ("IRN", "Iran",                  445.0,  32.0,   53.0, :GULF,     0.95),
    ("ARE", "United Arab Emirates",  385.0,  24.0,   54.0, :GULF,     0.75),
    ("SGP", "Singapore",             341.0,   1.35, 103.8, :EASIA,    0.00),
    ("EGY", "Egypt",                 235.0,  27.0,   30.0, :MED,      0.00),
    ("IRQ", "Iraq",                  193.0,  33.0,   44.0, :GULF,     1.00),
    ("QAT", "Qatar",                 161.0,  25.3,   51.2, :GULF,     1.00),
    ("KWT", "Kuwait",                121.0,  29.3,   47.6, :GULF,     1.00),
    ("OMN", "Oman",                   71.0,  21.0,   57.0, :GULF,     0.00),
    ("PAN", "Panama",                 62.0,   9.0,  -80.0, :AMERICAS, 0.00),
    ("BHR", "Bahrain",                35.0,  26.0,   50.5, :GULF,     1.00),
    ("MLT", "Malta",                  13.0,  35.9,   14.4, :EUR,      0.00),
    ("BRN", "Brunei",                 12.0,   4.5,  114.7, :EASIA,    0.00),
    ("GNQ", "Equatorial Guinea",      12.0,   1.6,   10.5, :WAFR,     0.00),
    ("ROW", "Rest of World",       28000.0,  10.0,   20.0, :ROW,      0.00),
]

const NC = length(COUNTRY_DATA)
const ISO    = [c[1] for c in COUNTRY_DATA]
const CNAME  = [c[2] for c in COUNTRY_DATA]
const VATGT  = [c[3] for c in COUNTRY_DATA]
const LAT    = [c[4] for c in COUNTRY_DATA]
const LON    = [c[5] for c in COUNTRY_DATA]
const REGION = [c[6] for c in COUNTRY_DATA]
const GULFSH = [c[7] for c in COUNTRY_DATA]

cid(iso) = findfirst(==(iso), ISO)

const ROUTE_NAME = ["Suez", "Panama", "Cape of Good Hope", "Direct",
                    "Hormuz-Suez", "Hormuz-Panama", "Hormuz-Cape", "Hormuz-Direct"]
const R_SUEZ, R_PANAMA, R_CAPE, R_DIRECT = 1, 2, 3, 4
const HORMUZ_ROUTES = 5:8
const ROUTE_DISTMULT = [1.15, 1.25, 1.75, 1.10, 1.22, 1.32, 1.80, 1.18]

# Baseline sea-route shares by region pair. These are DATA: with a large eta^R
# the implied baseline cost gaps are a fraction of a percent, so generating them
# from a distance logit would be meaningless. Hat algebra needs only the shares.
const PAIR_ROUTES = Dict{Set{Symbol},Vector{Float64}}(
    #                              Suez  Panama  Cape  Direct
    Set([:EUR, :EASIA])        => [0.85, 0.02, 0.13, 0.00],
    Set([:EUR, :SASIA])        => [0.86, 0.00, 0.12, 0.02],
    Set([:EUR, :GULF])         => [0.88, 0.00, 0.10, 0.02],
    Set([:MED, :EASIA])        => [0.88, 0.00, 0.10, 0.02],
    Set([:MED, :SASIA])        => [0.88, 0.00, 0.10, 0.02],
    Set([:MED, :GULF])         => [0.85, 0.00, 0.10, 0.05],
    Set([:EUR, :MED])          => [0.05, 0.00, 0.00, 0.95],
    Set([:EUR, :AMERICAS])     => [0.02, 0.08, 0.00, 0.90],
    Set([:MED, :AMERICAS])     => [0.05, 0.05, 0.00, 0.90],
    Set([:EASIA, :AMERICAS])   => [0.03, 0.32, 0.05, 0.60],
    Set([:SASIA, :AMERICAS])   => [0.30, 0.05, 0.15, 0.50],
    Set([:GULF, :AMERICAS])    => [0.35, 0.03, 0.27, 0.35],
    Set([:GULF, :EASIA])       => [0.00, 0.00, 0.00, 1.00],
    Set([:GULF, :SASIA])       => [0.00, 0.00, 0.00, 1.00],
    Set([:EASIA, :SASIA])      => [0.00, 0.00, 0.00, 1.00],
    Set([:WAFR, :EUR])         => [0.00, 0.00, 0.12, 0.88],
    Set([:WAFR, :MED])         => [0.00, 0.00, 0.10, 0.90],
    Set([:WAFR, :EASIA])       => [0.08, 0.00, 0.62, 0.30],
    Set([:WAFR, :SASIA])       => [0.02, 0.00, 0.60, 0.38],
    Set([:WAFR, :GULF])        => [0.02, 0.00, 0.55, 0.43],
    Set([:WAFR, :AMERICAS])    => [0.00, 0.05, 0.05, 0.90],
    Set([:ROW, :EUR])          => [0.20, 0.05, 0.08, 0.67],
    Set([:ROW, :MED])          => [0.15, 0.00, 0.05, 0.80],
    Set([:ROW, :EASIA])        => [0.18, 0.06, 0.10, 0.66],
    Set([:ROW, :SASIA])        => [0.18, 0.02, 0.10, 0.70],
    Set([:ROW, :GULF])         => [0.20, 0.02, 0.10, 0.68],
    Set([:ROW, :AMERICAS])     => [0.05, 0.12, 0.05, 0.78],
    Set([:ROW, :WAFR])         => [0.05, 0.02, 0.20, 0.73],
)

function passage_shares(ri::Symbol, rj::Symbol)
    v = get(PAIR_ROUTES, Set([ri, rj]), nothing)
    v === nothing && return [0.0, 0.0, 0.0, 1.0]
    return v
end

function gcdist(lat1, lon1, lat2, lon2)
    phi1, phi2 = deg2rad(lat1), deg2rad(lat2)
    dlam = deg2rad(lon2 - lon1)
    c = sin(phi1) * sin(phi2) + cos(phi1) * cos(phi2) * cos(dlam)
    return 6371.0 * acos(clamp(c, -1.0, 1.0))
end

# ============================================================================
# 3. BASELINE CONSTRUCTION
# ============================================================================
# Simpler than the full model in one respect: with no tolls there is no take to
# strip out, so the sourcing share that reaches producers IS pi. The raw
# assembled equilibrium is the baseline - there is no re-baselining step.

function build_distances()
    D = zeros(NC, NC)
    for i in 1:NC, j in 1:NC
        D[i, j] = i == j ? 400.0 : max(gcdist(LAT[i], LON[i], LAT[j], LON[j]), 400.0)
    end
    dr = zeros(NC, NC, NM, NR)
    for i in 1:NC, j in 1:NC
        for r in 1:NR
            dr[i, j, 1, r] = D[i, j] * ROUTE_DISTMULT[r]
        end
        dr[i, j, 2, R_DIRECT] = D[i, j]
        dr[i, j, 3, R_DIRECT] = D[i, j]
    end
    return D, dr
end

function build_transport_shares()
    mu0 = zeros(NC, NC, NS, NM)
    rho0 = zeros(NC, NC, NM, NR)
    for i in 1:NC, j in 1:NC
        for s in 1:NS
            if i == j
                mu0[i, j, s, 3] = 1.0
            else
                mu0[i, j, s, 1] = MODESHARE[s, 1]
                mu0[i, j, s, 2] = MODESHARE[s, 2]
                mu0[i, j, s, 3] = MODESHARE[s, 3]
            end
        end
        rho0[i, j, 2, R_DIRECT] = 1.0
        rho0[i, j, 3, R_DIRECT] = 1.0
        if i == j
            rho0[i, j, 1, R_DIRECT] = 1.0
        else
            base = passage_shares(REGION[i], REGION[j])
            # Hormuz exposure e_ij = g_i(1-g_j) + g_j(1-g_i): intra-Gulf trade
            # is unexposed. The strait is still geography here - it simply
            # cannot be priced.
            e = GULFSH[i] * (1 - GULFSH[j]) + GULFSH[j] * (1 - GULFSH[i])
            for b in 1:4
                rho0[i, j, 1, b]     = (1 - e) * base[b]
                rho0[i, j, 1, b + 4] = e * base[b]
            end
        end
    end
    return mu0, rho0
end

function initial_T()
    T = ones(NC, NS)
    energy = 2
    hydro = ["SAU", "IRN", "IRQ", "QAT", "KWT", "ARE", "OMN", "BHR", "BRN", "GNQ"]
    for iso in hydro
        i = cid(iso)
        T[i, 1:7] .= 0.05
        T[i, NS] = 0.6
    end
    for (iso, t) in [("SAU", 60.0), ("IRN", 25.0), ("IRQ", 22.0), ("QAT", 35.0),
                     ("KWT", 20.0), ("ARE", 22.0), ("OMN", 8.0), ("BHR", 4.0),
                     ("BRN", 12.0), ("GNQ", 12.0), ("USA", 4.0), ("ROW", 6.0)]
        T[cid(iso), energy] = t
    end
    for (iso, t) in [("CHN", 3.0), ("DEU", 2.5), ("JPN", 2.2), ("KOR", 2.2)]
        T[cid(iso), 5:7] .= t
    end
    T[cid("NLD"), :] .*= 1.4
    T[cid("SGP"), :] .*= 1.6
    T .*= exp.(0.20 .* randn(NC, NS))
    return T
end

function build_io()
    pat = [0.25 0.06 0.04 0.08 0.03 0.05 0.02 0.47;
           0.02 0.30 0.04 0.08 0.04 0.06 0.02 0.44;
           0.02 0.14 0.22 0.07 0.06 0.06 0.02 0.41;
           0.04 0.16 0.05 0.28 0.04 0.06 0.03 0.34;
           0.02 0.14 0.16 0.08 0.22 0.06 0.02 0.30;
           0.02 0.05 0.03 0.08 0.20 0.22 0.09 0.31;
           0.02 0.04 0.02 0.08 0.08 0.12 0.28 0.36;
           0.04 0.05 0.02 0.04 0.02 0.07 0.06 0.70]
    pat ./= sum(pat, dims = 2)
    g = zeros(NC, NS, NS)
    for j in 1:NC, k in 1:NS, s in 1:NS
        g[j, k, s] = (1 - GAMMA_VA[k]) * pat[k, s]
    end
    return g
end

function build_domestic_shares()
    lo, hi = log(minimum(VATGT)), log(maximum(VATGT))
    ds = zeros(NC, NS)
    for j in 1:NC
        z = (log(VATGT[j]) - lo) / (hi - lo)
        for s in 1:NS
            ds[j, s] = s == NS ? 0.95 : clamp(0.42 + 0.50 * z, 0.42, 0.93)
        end
    end
    return ds
end

function expand_rho(rho0)
    r = zeros(NC, NC, NS, NM, NR)
    @inbounds for i in 1:NC, j in 1:NC, s in 1:NS, m in 1:NM, rr in 1:NR
        r[i, j, s, m, rr] = rho0[i, j, m, rr]
    end
    return r
end

"""
    solve_levels(pi, gIO, afin, D, Ytot; ...) -> (Y, X, I)

Power iteration on the scale-indeterminate levels system, normalised so world
value added equals `Ytot`. Used ONLY to build the baseline, where the eigenvector
is what determines relative country sizes. No toll term: every dollar of
delivered spending reaches a producer.
"""
function solve_levels(pi_, gIO, afin, D, Ytot; Y0 = nothing,
                      maxiter = 4000, tol = 1e-10, damp = 0.9)
    Y = Y0 === nothing ? fill(Ytot / (NC * NS), NC, NS) : copy(Y0)
    X = zeros(NC, NS); I = zeros(NC); Ynew = zeros(NC, NS)
    for it in 1:maxiter
        @inbounds for j in 1:NC
            acc = 0.0
            for s in 1:NS
                acc += GAMMA_VA[s] * Y[j, s]
            end
            I[j] = acc + D[j]
        end
        @inbounds for j in 1:NC, s in 1:NS
            acc = afin[j, s] * I[j]
            for k in 1:NS
                acc += gIO[j, k, s] * Y[j, k]
            end
            X[j, s] = acc
        end
        fill!(Ynew, 0.0)
        @inbounds for s in 1:NS, j in 1:NC
            xjs = X[j, s]
            for i in 1:NC
                Ynew[i, s] += pi_[i, j, s] * xjs
            end
        end
        va = 0.0
        @inbounds for i in 1:NC, s in 1:NS
            va += GAMMA_VA[s] * Ynew[i, s]
        end
        Ynew .*= Ytot / va
        err = maximum(abs.(Ynew .- Y)) / (Ytot / (NC * NS))
        Y .= damp .* Ynew .+ (1 - damp) .* Y
        err < tol && break
    end
    @inbounds for j in 1:NC
        acc = 0.0
        for s in 1:NS
            acc += GAMMA_VA[s] * Y[j, s]
        end
        I[j] = acc + D[j]
    end
    @inbounds for j in 1:NC, s in 1:NS
        acc = afin[j, s] * I[j]
        for k in 1:NS
            acc += gIO[j, k, s] * Y[j, k]
        end
        X[j, s] = acc
    end
    return Y, X, I
end

"""
    solve_expenditure(pi, gIO, afin, LI, D; ...) -> (Y, X, I)

Same accounting with labour income `LI[i] = what[i] * VA[i]` taken as GIVEN.
Fixing final demand makes the Leontief map a strict contraction, so the system
has a unique solution and the outer wage loop can be driven by excess demand.
"""
function solve_expenditure(pi_, gIO, afin, LI, D; Y0 = nothing,
                           maxiter = 2000, tol = 1e-12)
    Y = Y0 === nothing ? fill(sum(LI) / (NC * NS), NC, NS) : copy(Y0)
    X = zeros(NC, NS); Ynew = zeros(NC, NS)
    I = LI .+ D
    scale = sum(abs, I) / (NC * NS)
    for it in 1:maxiter
        @inbounds for j in 1:NC, s in 1:NS
            acc = afin[j, s] * I[j]
            for k in 1:NS
                acc += gIO[j, k, s] * Y[j, k]
            end
            X[j, s] = acc
        end
        fill!(Ynew, 0.0)
        @inbounds for s in 1:NS, j in 1:NC
            xjs = X[j, s]
            for i in 1:NC
                Ynew[i, s] += pi_[i, j, s] * xjs
            end
        end
        err = maximum(abs.(Ynew .- Y)) / scale
        Y .= Ynew
        err < tol && break
    end
    @inbounds for j in 1:NC, s in 1:NS
        acc = afin[j, s] * I[j]
        for k in 1:NS
            acc += gIO[j, k, s] * Y[j, k]
        end
        X[j, s] = acc
    end
    return Y, X, I
end

function build_baseline(; verbose = true)
    Dist, distroute = build_distances()
    mu0, rho0 = build_transport_shares()
    rho = expand_rho(rho0)
    mu = mu0
    gIO = build_io()
    ds = build_domestic_shares()
    afin = repeat(ALPHA_FIN', NC, 1)
    D = zeros(NC)
    T = initial_T()
    Ytot = sum(VATGT)

    pi_ = zeros(NC, NC, NS)
    Y = nothing; X = nothing; I = nothing
    for outer_it in 1:80
        @inbounds for j in 1:NC, s in 1:NS
            tot = 0.0
            for i in 1:NC
                if i == j
                    pi_[i, j, s] = 0.0
                else
                    v = T[i, s] * (Dist[i, j] / 1000)^(-1.0)
                    pi_[i, j, s] = v
                    tot += v
                end
            end
            scale = (1 - ds[j, s]) / max(tot, 1e-300)
            for i in 1:NC
                i != j && (pi_[i, j, s] *= scale)
            end
            pi_[j, j, s] = ds[j, s]
        end
        Y, X, I = solve_levels(pi_, gIO, afin, D, Ytot; Y0 = Y)
        VAimp = vec(sum(GAMMA_VA' .* Y, dims = 2))
        ratio = VATGT ./ VAimp
        err = maximum(abs.(log.(ratio)))
        T .*= ratio .^ 0.5
        if err < 5e-3
            verbose && @printf("  country-size fixed point converged in %d iterations (max log gap %.4f)\n", outer_it, err)
            break
        end
    end

    VA = vec(sum(GAMMA_VA' .* Y, dims = 2))
    Xijs = zeros(NC, NC, NS)
    @inbounds for i in 1:NC, j in 1:NC, s in 1:NS
        Xijs[i, j, s] = pi_[i, j, s] * X[j, s]
    end

    Xibar = zeros(NM, NR); Xrs = zeros(NR, NS)
    @inbounds for i in 1:NC, j in 1:NC, s in 1:NS
        x = Xijs[i, j, s]
        x == 0.0 && continue
        for m in 1:NM, r in 1:NR
            f = x * mu[i, j, s, m] * rho[i, j, s, m, r]
            f == 0.0 && continue
            Xibar[m, r] += distroute[i, j, m, r] * TPD[s] * f
            m == 1 && (Xrs[r, s] += f)
        end
    end
    Xibar[3, :] .= 0.0
    Psibar = sum(Xibar[1, :])
    Xiactive = Xibar .> 0
    Xibar[.!Xiactive] .= 1.0

    return (; pi_, mu, rho, gIO, afin, D, Dist, distroute,
            Y, X, I, VA, Xijs, Xrs, Xibar, Xiactive, Psibar, Ytot)
end

# ============================================================================
# 4. TRANSPORT BLOCK  (eq 14-20, in changes)
# ============================================================================

@inline function lse(v::AbstractVector{Float64})
    mx = -Inf
    @inbounds for x in v
        x > mx && (mx = x)
    end
    isfinite(mx) || return -Inf
    acc = 0.0
    @inbounds for x in v
        isfinite(x) && (acc += exp(x - mx))
    end
    return mx + log(acc)
end

mutable struct TransportState
    Xihat::Matrix{Float64}
    Psihat::Float64
    dhat::Array{Float64,3}
    dmode::Array{Float64,4}        # mode-level composite change, for diagnostics
    mu::Array{Float64,4}
    rho::Array{Float64,5}
    Xrs::Matrix{Float64}
    damp::Float64
    resid::Float64
end

TransportState() = TransportState(ones(NM, NR), 1.0, ones(NC, NC, NS),
                                  ones(NC, NC, NS, NM),
                                  zeros(NC, NC, NS, NM), zeros(NC, NC, NS, NM, NR),
                                  zeros(NR, NS), 0.4, Inf)

"""
    transport_pass!(st, base, etaR, Xijs, route_open, dbarhat; niter, tol)

Iterate the transport fixed point (18) -> (14)-(17) -> (19)-(20), holding
bilateral sectoral flows `Xijs` fixed. `dbarhat` is an (NS, NM, NR) array of
EXOGENOUS baseline freight-wedge shocks - war-risk premia, drought
restrictions - expressed already in cargo-iceberg terms.

There is no toll wedge in this variant, so the only forces moving route costs
are the congestion and capacity aggregates and whatever `dbarhat` imposes.

The congestion aggregates use an ADAPTIVE step. The map Xi -> Xi has loop gain
of order etaR * chi_s * lambda_s * (1 - rho_r); at large etaR that exceeds one
for high-freight-share sectors and a fixed step oscillates.
"""
function transport_pass!(st::TransportState, base, etaR::Float64, Xijs,
                         route_open::Vector{Bool}, dbarhat::Array{Float64,3};
                         niter::Int = 3, tol::Float64 = 1e-12)
    rho0, mu0 = base.rho, base.mu
    Xibar, Psibar, distroute = base.Xibar, base.Psibar, base.distroute

    g = zeros(NS, NM, NR)
    lnrho0 = zeros(NR); lnrhonew = zeros(NR)
    lnmu0 = zeros(NM); lndm = zeros(NM); tmpm = zeros(NM); tmpr = zeros(NR)
    Xinew = zeros(NM, NR)

    for _ in 1:niter
        lnPsi = log(st.Psihat)
        @inbounds for s in 1:NS, m in 1:NM, r in 1:NR
            sea = (m == 1)
            g[s, m, r] = log(dbarhat[s, m, r]) +
                         (sea ? CHI[s] * (LAMBDA[s] * log(st.Xihat[m, r]) + GAMMA_G[s] * lnPsi) : 0.0)
        end

        fill!(Xinew, 0.0); fill!(st.Xrs, 0.0)
        @inbounds for i in 1:NC, j in 1:NC
            for s in 1:NS
                for m in 1:NM
                    for r in 1:NR
                        p = rho0[i, j, s, m, r]
                        lnrho0[r] = (p > 0 && (m != 1 || route_open[r])) ? log(p) : -Inf
                        tmpr[r] = lnrho0[r] - etaR * g[s, m, r]
                    end
                    L = lse(tmpr)
                    if isfinite(L)
                        lndm[m] = -L / etaR
                        for r in 1:NR
                            lnrhonew[r] = lnrho0[r] - etaR * (g[s, m, r] - lndm[m])
                            st.rho[i, j, s, m, r] = isfinite(lnrhonew[r]) ? exp(lnrhonew[r]) : 0.0
                        end
                    else
                        lndm[m] = Inf
                        for r in 1:NR
                            st.rho[i, j, s, m, r] = 0.0
                        end
                    end
                    st.dmode[i, j, s, m] = isfinite(lndm[m]) ? exp(lndm[m]) : Inf
                    p = mu0[i, j, s, m]
                    lnmu0[m] = (p > 0 && isfinite(lndm[m])) ? log(p) : -Inf
                end
                etaM = ETA_M[s]
                for m in 1:NM
                    tmpm[m] = isfinite(lnmu0[m]) ? lnmu0[m] - etaM * lndm[m] : -Inf
                end
                L = lse(tmpm)
                if isfinite(L)
                    lnds = -L / etaM
                    st.dhat[i, j, s] = exp(lnds)
                    for m in 1:NM
                        st.mu[i, j, s, m] = isfinite(lnmu0[m]) ? exp(lnmu0[m] - etaM * (lndm[m] - lnds)) : 0.0
                    end
                else
                    st.dhat[i, j, s] = Inf
                    for m in 1:NM
                        st.mu[i, j, s, m] = 0.0
                    end
                end
                x = Xijs[i, j, s]
                for m in 1:NM
                    mum = st.mu[i, j, s, m]
                    mum == 0.0 && continue
                    for r in 1:NR
                        rr = st.rho[i, j, s, m, r]
                        rr == 0.0 && continue
                        f = x * mum * rr
                        if m == 1
                            Xinew[m, r] += distroute[i, j, m, r] * TPD[s] * f
                            st.Xrs[r, s] += f
                        elseif m == 2
                            Xinew[m, r] += distroute[i, j, m, r] * TPD[s] * f
                        end
                    end
                end
            end
        end

        Psinew = sum(@view Xinew[1, :]) / Psibar
        resid = 0.0
        d = st.damp
        @inbounds for m in 1:NM, r in 1:NR
            if base.Xiactive[m, r]
                tgt = max(Xinew[m, r] / Xibar[m, r], 1e-3)
                resid = max(resid, abs(log(tgt) - log(st.Xihat[m, r])))
                st.Xihat[m, r] = d * tgt + (1 - d) * st.Xihat[m, r]
            else
                st.Xihat[m, r] = 1.0
            end
        end
        tgtPsi = max(Psinew, 1e-3)
        resid = max(resid, abs(log(tgtPsi) - log(st.Psihat)))
        st.Psihat = d * tgtPsi + (1 - d) * st.Psihat
        st.damp = resid > st.resid ? max(0.5 * st.damp, 0.01) : min(1.06 * st.damp, 0.6)
        st.resid = resid
        resid < tol && break
    end
    return st
end

# ============================================================================
# 5. PRICE BLOCK AND COUNTERFACTUAL SOLVER
# ============================================================================

function solve_prices!(Phat, base, what, dhat; maxiter = 400, tol = 1e-12)
    pi_, gIO = base.pi_, base.gIO
    chat = ones(NC, NS)
    for it in 1:maxiter
        @inbounds for i in 1:NC, k in 1:NS
            acc = GAMMA_VA[k] * log(what[i])
            for s in 1:NS
                acc += gIO[i, k, s] * log(Phat[i, s])
            end
            chat[i, k] = exp(acc)
        end
        err = 0.0
        @inbounds for j in 1:NC, s in 1:NS
            th = THETA[s]; acc = 0.0
            for i in 1:NC
                p = pi_[i, j, s]
                p == 0.0 && continue
                cd = chat[i, s] * dhat[i, j, s]
                isfinite(cd) || continue
                acc += p * cd^(-th)
            end
            new = acc <= 0.0 ? 1e10 : acc^(-1 / th)
            err = max(err, abs(log(new) - log(Phat[j, s])))
            Phat[j, s] = new
        end
        err < tol && break
    end
    return Phat, chat
end

"""
    solve_counterfactual(base, sh; etaR, ...) -> NamedTuple

Exact-hat solver. `sh` carries `route_open :: Vector{Bool}` and
`dbarhat :: Array{Float64,3}` (NS, NM, NR). No `phi_new` field: there are no
tolls to set.

Structure: outer loop on wages driven by excess demand; inside it a transport
pass, the price fixed point, and the expenditure system. Closure: nominal
deficits fixed at baseline, world value added as numeraire.
"""
function solve_counterfactual(base, sh; etaR::Float64 = ETA_R_REF,
                              maxiter = 4000, tol = 1e-9, damp = 0.25, verbose = false)
    what = ones(NC)
    Phat = ones(NC, NS)
    st = TransportState()
    st.mu .= base.mu; st.rho .= base.rho
    Xijs = copy(base.Xijs)
    Y = copy(base.Y); X = copy(base.X); I = copy(base.I)

    conv = false; iters = 0
    step = damp; err_prev = Inf
    for it in 1:maxiter
        iters = it
        transport_pass!(st, base, etaR, Xijs, sh.route_open, sh.dbarhat;
                        niter = it <= 5 ? 200 : 30)
        _, chat = solve_prices!(Phat, base, what, st.dhat; maxiter = it <= 5 ? 400 : 60)

        pinew = similar(base.pi_)
        @inbounds for i in 1:NC, j in 1:NC, s in 1:NS
            p = base.pi_[i, j, s]
            if p == 0.0
                pinew[i, j, s] = 0.0
            else
                cd = chat[i, s] * st.dhat[i, j, s]
                pinew[i, j, s] = isfinite(cd) ? p * (cd / Phat[j, s])^(-THETA[s]) : 0.0
            end
        end

        LI = what .* base.VA
        Y, X, I = solve_expenditure(pinew, base.gIO, base.afin, LI, base.D;
                                    Y0 = Y, maxiter = it <= 5 ? 2000 : 300)

        ED = vec(sum(GAMMA_VA' .* Y, dims = 2)) .- LI
        err = maximum(abs.(ED ./ LI))
        if it > 1
            step = err < err_prev ? min(step * 1.08, 0.9) : max(step * 0.5, 0.02)
        end
        err_prev = err
        what .*= (1 .+ step .* ED ./ LI)
        what .*= sum(base.VA) / sum(what .* base.VA)

        @inbounds for i in 1:NC, j in 1:NC, s in 1:NS
            Xijs[i, j, s] = pinew[i, j, s] * X[j, s]
        end

        if err < tol
            conv = true
            verbose && @printf("  >>>> GE converged in %d iterations (max |ED/LI| = %.2e) <<<<\n", it, err)
            break
        end
    end
    conv || @warn "GE did not converge within maxiter"

    Pfin = [exp(sum(base.afin[j, s] * log(Phat[j, s]) for s in 1:NS)) for j in 1:NC]
    What = (I ./ base.I) ./ Pfin
    dW = 100 .* (What .- 1)
    Wworld = (sum(I) / sum(base.I)) / exp(sum(base.I[j] / sum(base.I) * log(Pfin[j]) for j in 1:NC))

    pinew = similar(base.pi_)
    @inbounds for i in 1:NC, j in 1:NC, s in 1:NS
        pinew[i, j, s] = Xijs[i, j, s] / X[j, s]
    end

    route_flow = vec(sum(st.Xrs, dims = 2))
    route_flow_base = vec(sum(base.Xrs, dims = 2))

    return (; what, Phat, Pfin, dhat = copy(st.dhat), dmode = copy(st.dmode),
            mu = copy(st.mu), rho = copy(st.rho), pi_ = pinew, Y, X, I,
            Xrs = copy(st.Xrs), Xihat = copy(st.Xihat), Psihat = st.Psihat,
            Xijs = copy(Xijs), dW, What, Wworld = 100 * (Wworld - 1),
            route_flow, route_flow_base,
            route_change = 100 .* (route_flow ./ max.(route_flow_base, 1e-12) .- 1),
            converged = conv, iters)
end

# ---- shocks ---------------------------------------------------------------

null_shock() = (route_open = fill(true, NR), dbarhat = ones(NS, NM, NR))

"Close a set of sea routes outright."
function shock_close(routes::Vector{Int})
    ro = fill(true, NR)
    for r in routes
        ro[r] = false
    end
    (route_open = ro, dbarhat = ones(NS, NM, NR))
end

"""
    shock_freight(routes, fhat) -> shock

Raise the FREIGHT WEDGE on `routes` by the factor `fhat` - a war-risk premium, a
drought restriction, a convoy surcharge. This is the closest thing this variant
has to the full model's transit fee, and it is economically distinct: nobody
collects it, so it is a pure resource cost rather than a transfer.

The cargo-value iceberg moves by `fhat^chi_s` (eq 5), which is why a given
premium bites hard on bulk and barely at all on electronics.
"""
function shock_freight(routes::Vector{Int}, fhat::Float64)
    db = ones(NS, NM, NR)
    for s in 1:NS, r in routes
        db[s, 1, r] = fhat^CHI[s]
    end
    (route_open = fill(true, NR), dbarhat = db)
end

# convenience: the Red Sea closure kills Suez AND its Hormuz composite
shock_red_sea() = shock_close([R_SUEZ, R_SUEZ + 4])
shock_hormuz_premium(fhat) = shock_freight(collect(HORMUZ_ROUTES), fhat)
shock_suez_premium(fhat) = shock_freight([R_SUEZ, R_SUEZ + 4], fhat)
shock_panama_premium(fhat) = shock_freight([R_PANAMA, R_PANAMA + 4], fhat)

# ============================================================================
# 6. VERIFICATION
# ============================================================================

check(label, ok) = (@printf("  [%s] %s\n", ok ? " OK " : "FAIL", label); flush(stdout); ok)

function checkv(label, val, tol)
    ok = val < tol
    @printf("  [%s] %-46s  residual %.3e  (tol %.0e)\n", ok ? " OK " : "FAIL", label, val, tol)
    flush(stdout)
    return ok
end

function verify_baseline(b)
    println("\n--- baseline consistency ---")
    ok = true
    ok &= check("sourcing shares sum to one",
                maximum(abs.(vec(sum(b.pi_, dims = 1)) .- 1)) < 1e-9)
    ok &= check("mode shares sum to one",
                maximum(abs.(vec(sum(b.mu, dims = 4)) .- 1)) < 1e-9)
    ok &= check("route shares sum to one",
                maximum(abs.(vec(sum(b.rho, dims = 5)) .- 1)) < 1e-9)
    Yc = zeros(NC, NS)
    @inbounds for s in 1:NS, j in 1:NC, i in 1:NC
        Yc[i, s] += b.pi_[i, j, s] * b.X[j, s]
    end
    ok &= check("goods market clearing  Y = sum_j pi X",
                maximum(abs.(Yc .- b.Y)) / maximum(b.Y) < 1e-7)
    ok &= check("income identity  I = VA + D  (no rents in this variant)",
                maximum(abs.(b.I .- (b.VA .+ b.D))) / maximum(b.I) < 1e-9)
    ok &= check("world VA equals target", abs(sum(b.VA) / b.Ytot - 1) < 1e-6)
    return ok
end

function verify_null_shock(b, etaR)
    println("\n--- null shock reproduces the baseline ---")
    res = solve_counterfactual(b, null_shock(); etaR = etaR, tol = 1e-10, damp = 0.5)
    ok = true
    ok &= checkv("what = 1", maximum(abs.(res.what .- 1)), 1e-6)
    ok &= checkv("Phat = 1", maximum(abs.(res.Phat .- 1)), 1e-6)
    ok &= checkv("dhat = 1", maximum(abs.(res.dhat .- 1)), 1e-8)
    ok &= checkv("route shares unchanged", maximum(abs.(res.rho .- b.rho)), 1e-8)
    ok &= checkv("welfare change zero", maximum(abs.(res.dW)), 1e-4)
    ok &= checkv("baseline Xi reproduced", maximum(abs.(res.Xihat .- 1)), 1e-8)
    return ok
end

"""
Closed-form check of the route nest. With congestion and capacity switched off
and no exogenous freight shock, closing route r0 must leave the surviving route
shares exactly renormalised,

    rho'_r = rho_r / (1 - rho_r0),

and raise the within-mode composite by exactly (1 - rho_r0)^(-1/eta^R). Both are
independent of everything else in the model, so this isolates (14)-(15).
"""
function verify_route_nest_closed_form(b, etaR)
    println("\n--- closed form of the route nest (eq 14-15) ---")
    savel, saveg = copy(LAMBDA), copy(GAMMA_G)
    local res
    try
        LAMBDA .= 0.0; GAMMA_G .= 0.0
        res = solve_counterfactual(b, shock_red_sea(); etaR = etaR, tol = 1e-10, damp = 0.4)
    finally
        LAMBDA .= savel; GAMMA_G .= saveg
    end
    errshare = 0.0; errcomp = 0.0
    @inbounds for i in 1:NC, j in 1:NC, s in 1:NS
        closed = b.rho[i, j, s, 1, R_SUEZ] + b.rho[i, j, s, 1, R_SUEZ + 4]
        open_ = 1 - closed
        open_ < 1e-9 && continue                       # pair entirely cut off
        for r in 1:NR
            (r == R_SUEZ || r == R_SUEZ + 4) && continue
            b.rho[i, j, s, 1, r] == 0.0 && continue
            errshare = max(errshare, abs(res.rho[i, j, s, 1, r] - b.rho[i, j, s, 1, r] / open_))
        end
        errcomp = max(errcomp, abs(res.dmode[i, j, s, 1] - open_^(-1 / etaR)))
    end
    ok = checkv("surviving route shares renormalise exactly", errshare, 1e-10)
    ok &= checkv("composite equals (1 - rho_closed)^(-1/etaR)", errcomp, 1e-10)
    println("  Note this is independent of chi_s: a route CLOSURE costs the same whatever the")
    println("  freight share, because the baseline cost gaps are inferred from observed shares.")
    println("  chi_s scales the congestion and premium channels only.")
    return ok
end

function verify_table7()
    println("\n--- Appendix B, Table 7: multiplicative vs additive pass-through ---")
    phi = 1.5855
    rows = [("Electronics, pharmaceuticals", 0.01, 1.0046, 1.0059),
            ("General manufactures",         0.03, 1.0139, 1.0176),
            ("Steel, basic chemicals",       0.08, 1.0376, 1.0468),
            ("Bulk (grain, ore, fertilizer)",0.30, 1.1483, 1.1756),
            ("Shipping-only limit",          1.00, 1.5855, 1.5855)]
    ok = true
    @printf("  %-32s %6s %10s %10s %10s %10s\n", "sector", "chi", "d=f^chi", "paper", "d=1+chi(.)", "paper")
    for (nm, chi, pm, pa) in rows
        m = phi^chi; a = 1 + chi * (phi - 1)
        @printf("  %-32s %6.2f %10.4f %10.4f %10.4f %10.4f\n", nm, chi, m, pm, a, pa)
        ok &= (abs(m - pm) < 5e-4) && (abs(a - pa) < 5e-4)
    end
    return check("Table 7 reproduced", ok)
end

function verify_limits(b, etaR)
    println("\n--- limiting cases ---")
    ok = true
    savechi = copy(CHI)
    try
        CHI .= 0.0
        r0 = solve_counterfactual(b, shock_suez_premium(1.5); etaR = etaR, tol = 1e-10, damp = 0.5)
        ok &= checkv("chi = 0  =>  a freight premium has zero effect", maximum(abs.(r0.dW)), 1e-6)
    finally
        CHI .= savechi
    end
    savel, saveg = copy(LAMBDA), copy(GAMMA_G)
    local dfull, dnc
    try
        rfull = solve_counterfactual(b, shock_red_sea(); etaR = etaR, tol = 1e-9, damp = 0.35)
        dfull = maximum(abs.(filter(isfinite, rfull.dhat) .- 1))
        LAMBDA .= 0.0; GAMMA_G .= 0.0
        rnc = solve_counterfactual(b, shock_red_sea(); etaR = etaR, tol = 1e-9, damp = 0.35)
        dnc = maximum(abs.(filter(isfinite, rnc.dhat) .- 1))
    finally
        LAMBDA .= savel; GAMMA_G .= saveg
    end
    @printf("  max |dhat - 1| under the closure: with congestion %.5f, without %.5f\n", dfull, dnc)
    ok &= check("congestion amplifies the delivered-cost response", dfull > dnc)
    return ok
end

# ============================================================================
# 7. REPORTING
# ============================================================================

const REPORT_GROUPS = [
    ("Suez-dependent / transshipment", ["MLT", "SGP", "EGY", "ITA"]),
    ("Gulf economies",                 ["QAT", "KWT", "IRQ", "BHR", "ARE", "SAU", "OMN", "IRN"]),
    ("Large trading economies",        ["DEU", "FRA", "NLD", "USA", "CHN", "JPN", "KOR", "IND"]),
    ("Other",                          ["PAN", "GNQ", "BRN", "ROW"]),
]

function table_welfare(results, labels)
    println()
    @printf("%-24s", "")
    for l in labels
        @printf("%13s", l)
    end
    println()
    println("-"^(24 + 13 * length(labels)))
    for (gname, isos) in REPORT_GROUPS
        @printf("%-24s\n", gname)
        for iso in isos
            i = cid(iso)
            @printf("  %-22s", CNAME[i])
            for r in results
                @printf("%13.3f", r.dW[i])
            end
            println()
        end
    end
    @printf("  %-22s", "World")
    for r in results
        @printf("%13.3f", r.Wworld)
    end
    println()
end

function route_pct(res)
    out = Float64[]
    for r in 1:4
        fb = res.route_flow_base[r] + res.route_flow_base[r + 4]
        fn = res.route_flow[r] + res.route_flow[r + 4]
        push!(out, fb > 0 ? 100 * (fn / fb - 1) : 0.0)
    end
    fb = sum(res.route_flow_base[HORMUZ_ROUTES]); fn = sum(res.route_flow[HORMUZ_ROUTES])
    push!(out, fb > 0 ? 100 * (fn / fb - 1) : 0.0)
    return out
end

function table_routes(results, labels)
    rn = ["Suez", "Panama", "Cape of Good Hope", "Direct", "Hormuz transit"]
    println()
    @printf("%-24s", "Route")
    for l in labels
        @printf("%13s", l)
    end
    println()
    println("-"^(24 + 13 * length(labels)))
    pcts = [route_pct(r) for r in results]
    for k in 1:5
        @printf("%-24s", rn[k])
        for p in pcts
            @printf("%13.1f", p[k])
        end
        println()
    end
end

function grouped_bar(M::AbstractMatrix, xlabels, labels; kw...)
    n, g = size(M)
    w = 0.8 / g
    p = plot(; kw...)
    for k in 1:g
        xs = (1:n) .- 0.4 .+ w * (k - 0.5)
        bar!(p, xs, M[:, k]; bar_width = w, label = labels[k], linewidth = 0.3)
    end
    xticks!(p, 1:n, xlabels)
    return p
end

# ============================================================================
# 8. MAIN
# ============================================================================

banner("BHS (2026) WITHOUT PRICED PASSAGES - ENDOGENOUS FREIGHT ONLY")
println("Caliendo-Parro + nested mode-over-route transport block with congestion")
println("and a global capacity constraint. No authorities, no tolls, no rents.")
println("Stylized calibration: $NC countries x $NS sectors x $NM modes x $NR routes")

banner("1. BASELINE")
@time "  build_baseline" base = build_baseline()
flush(stdout)
@printf("  world value added        : \$%.1f tn\n", sum(base.VA) / 1000)
intl = sum(base.Xijs) - sum(base.Xijs[i, i, s] for i in 1:NC, s in 1:NS)
@printf("  international trade      : \$%.1f tn\n", intl / 1000)
seatrade = sum(base.Xrs)
@printf("  sea-borne trade          : \$%.1f tn\n", seatrade / 1000)
@printf("  via Suez                 : \$%.0f bn\n",
        sum(base.Xrs[R_SUEZ, :]) + sum(base.Xrs[R_SUEZ + 4, :]))
@printf("  transiting Hormuz        : \$%.0f bn\n", sum(base.Xrs[HORMUZ_ROUTES, :]))
chibar = sum(CHI[s] * sum(base.Xrs[:, s]) for s in 1:NS) / seatrade
@printf("  trade-weighted chi (sea) : %.3f\n", chibar)
println("\n  No toll revenue appears anywhere: income is value added plus the fixed")
println("  deficit, and every dollar of delivered spending reaches a producer.")
verify_baseline(base)

banner("2. VERIFICATION")
verify_null_shock(base, ETA_R_REF)
verify_route_nest_closed_form(base, ETA_R_REF)
verify_table7()
verify_limits(base, ETA_R_REF)

get(ENV, "BHS_QUICK", "0") == "1" && exit(0)

banner("3. COUNTERFACTUALS")
println("Available experiments are CLOSURES and EXOGENOUS COST SHOCKS. A transit")
println("fee is not available here - a fee is a toll, which this variant drops.")
say("\nSolving Red Sea closure ...")
res_rs = solve_counterfactual(base, shock_red_sea(); etaR = ETA_R_REF, tol = 1e-9)
say("Solving Suez war-risk premium (freight wedge x1.5) ...")
res_sp = solve_counterfactual(base, shock_suez_premium(1.5); etaR = ETA_R_REF, tol = 1e-9)
say("Solving Hormuz war-risk premium (freight wedge x1.5) ...")
res_hp = solve_counterfactual(base, shock_hormuz_premium(1.5); etaR = ETA_R_REF, tol = 1e-9)
say("Solving Panama restriction (freight wedge x1.5) ...")
res_pp = solve_counterfactual(base, shock_panama_premium(1.5); etaR = ETA_R_REF, tol = 1e-9)

allres = [res_rs, res_sp, res_hp, res_pp]
lbls = ["Red Sea shut", "Suez x1.5", "Hormuz x1.5", "Panama x1.5"]

banner("3.1 WELFARE (percent change in real income)")
table_welfare(allres, lbls)

println("\n  Five largest losses and gains per scenario:")
for (lab, r) in zip(lbls, allres)
    ord = sortperm(r.dW)
    @printf("    %-14s losses: %s\n", lab,
            join([@sprintf("%s %.3f", ISO[i], r.dW[i]) for i in ord[1:5]], ", "))
    @printf("    %-14s gains : %s\n", "",
            join([@sprintf("%s %+.3f", ISO[i], r.dW[i]) for i in reverse(ord[end-4:end])], ", "))
end

banner("3.2 SEABORNE TRADE BY ROUTE (percent change)")
table_routes(allres, lbls)

banner("4. WHO BEARS IT, WITHOUT A COLLECTOR")
suez_exp = [sum(base.Xijs[i, j, s] * base.mu[i, j, s, 1] *
                (base.rho[i, j, s, 1, R_SUEZ] + base.rho[i, j, s, 1, R_SUEZ + 4])
                for j in 1:NC, s in 1:NS) / base.I[i] for i in 1:NC]
ord = sortperm(res_rs.dW)
@printf("  %-24s %12s %14s\n", "", "loss (%)", "Suez exp. (%)")
println("  " * "-"^52)
for i in ord[1:8]
    @printf("  %-24s %12.3f %14.1f\n", CNAME[i], res_rs.dW[i], 100 * suez_exp[i])
end
rho_corr = cor(suez_exp, res_rs.dW)
@printf("\n  correlation(Suez exposure, welfare change) = %.3f\n", rho_corr)
println("  In the FULL model this correlation is weak, because incidence is driven by")
println("  who collects the rent rather than by whose cargo is rerouted. With no")
println("  collector, exposure predicts incidence again - that contrast is the whole")
println("  content of the paper's Section 6.3, isolated here as a standalone model.")

banner("5. THE PARAMETER THIS VARIANT CANNOT IDENTIFY")
println("Closing a route raises the within-mode composite by (1-rho_closed)^(-1/etaR),")
println("so the delivered-cost shock and the welfare loss scale as 1/etaR. With no")
println("toll revenue to match, nothing pins it down. The sweep below is the honest")
println("statement of what the model can and cannot say.\n")
etagrid = [50.0, 80.0, 120.0, 180.0, 250.0, 350.0, 500.0, 800.0]
mlt = cid("MLT")
sweep = Float64[]; sweep_world = Float64[]; sweep_cost = Float64[]
for e in etagrid
    r = solve_counterfactual(base, shock_red_sea(); etaR = e, tol = 1e-9)
    push!(sweep, r.dW[mlt]); push!(sweep_world, r.Wworld)
    push!(sweep_cost, maximum(abs.(filter(isfinite, r.dhat) .- 1)))
end
# the same sweep with the congestion and capacity channels switched off, which
# is where the closed form (1 - rho_closed)^(-1/etaR) applies exactly
sweep_nc = Float64[]
let savel = copy(LAMBDA), saveg = copy(GAMMA_G)
    try
        LAMBDA .= 0.0; GAMMA_G .= 0.0
        for e in etagrid
            r = solve_counterfactual(base, shock_red_sea(); etaR = e, tol = 1e-9)
            push!(sweep_nc, r.Wworld)
        end
    finally
        LAMBDA .= savel; GAMMA_G .= saveg
    end
end

prod_c = abs.(sweep_world) .* etagrid
prod_nc = abs.(sweep_nc) .* etagrid
@printf("  %8s %14s %14s %14s %16s\n", "eta^R", "World dW (%)", "  no congest.", "Malta dW (%)", "max |dhat-1|")
for (k, e) in enumerate(etagrid)
    @printf("  %8.0f %14.4f %14.4f %14.4f %16.5f\n",
            e, sweep_world[k], sweep_nc[k], sweep[k], sweep_cost[k])
end

spread(v) = maximum(v) / minimum(v) - 1
@printf("\n  world loss at eta^R = 50 is %.1fx the loss at eta^R = 800 (16x range in eta^R)\n",
        sweep_world[1] / sweep_world[end])
@printf("  |world loss| x eta^R varies by %3.0f%% with congestion OFF  -> the closed form's\n",
        100 * spread(prod_nc))
@printf("                        and by %3.0f%% with congestion ON\n", 100 * spread(prod_c))
println("\n  Reading: with congestion off the loss is proportional to 1/eta^R exactly, as")
println("  (1 - rho_closed)^(-1/etaR) implies. Switching congestion on puts a FLOOR under")
println("  the cost - the volume rerouted onto the Cape is essentially independent of")
println("  eta^R, so the congestion component does not shrink as route choice gets more")
println("  elastic. The model is therefore less hostage to the unidentified parameter")
println("  than the closed form alone suggests, but a 16x range in eta^R still moves the")
println("  headline number by a factor of nearly three.")

banner("6. FIGURES")

p1 = plot(etagrid, abs.(sweep_world), lw = 2.5, marker = :circle, ms = 4,
          label = "world, with congestion", xscale = :log10, yscale = :log10,
          legend = :bottomleft,
          xlabel = "route-choice elasticity eta^R (log scale)",
          ylabel = "real income loss (%, absolute, log scale)",
          title = "What the missing toll moment costs you")
plot!(p1, etagrid, abs.(sweep_nc), lw = 2.5, ls = :dash, marker = :diamond, ms = 4,
      label = "world, no congestion (exactly 1/eta^R)")
plot!(p1, etagrid, abs.(sweep), lw = 2.5, marker = :square, ms = 4, label = "Malta")
vline!(p1, [ETA_R_REF], ls = :dot, lc = :grey, label = "paper's toll-calibrated 350")
savefig(p1, joinpath(GRAPHDIR, "nt_eta_sensitivity.pdf"))

function plot_extremes(res, ttl, fname)
    o = sortperm(res.dW)
    sel = vcat(o[1:5], reverse(o[end-4:end]))
    vals = res.dW[sel]
    bar(1:10, vals, orientation = :h, yticks = (1:10, ISO[sel]), legend = false,
        color = [v < 0 ? :firebrick : :steelblue for v in vals],
        xlabel = "change in real income (%)", title = ttl, yflip = true)
    savefig(joinpath(GRAPHDIR, fname))
end
plot_extremes(res_rs, "Red Sea closure, no tolls: largest losses and gains", "nt_extremes_redsea.pdf")
plot_extremes(res_hp, "Hormuz freight premium x1.5: largest losses and gains", "nt_extremes_hormuz.pdf")

M = hcat(route_pct(res_rs), route_pct(res_sp), route_pct(res_hp), route_pct(res_pp))
p3 = grouped_bar(M, ["Suez", "Panama", "Cape", "Direct", "Hormuz"], lbls;
                 ylabel = "change in flow (%)", legend = :bottomleft,
                 title = "Seaborne trade by route")
savefig(p3, joinpath(GRAPHDIR, "nt_route_reallocation.pdf"))

p4 = scatter(100 .* suez_exp, res_rs.dW, ms = 5, legend = false,
             xlabel = "Suez-routed exports / income (%)",
             ylabel = "change in real income (%)",
             title = @sprintf("Without a collector, exposure predicts incidence (r = %.2f)", rho_corr),
             series_annotations = text.(ISO, 5, :bottom))
savefig(p4, joinpath(GRAPHDIR, "nt_exposure_scatter.pdf"))

secexp = [sum(base.Xrs[:, s]) for s in 1:NS]
p5 = bar(1:NS, CHI, legend = false, xticks = (1:NS, SECTORS), xrotation = 35,
         ylabel = "freight cost share chi_s",
         title = "Pass-through by sector: what makes a passage shock bite")
savefig(p5, joinpath(GRAPHDIR, "nt_freight_shares.pdf"))

println("  6 figures written to $(GRAPHDIR)")

banner("7. CHECKS ON THE MECHANISM")
check("closure costs fall on route-lengthened economies, not a collector",
      argmin(res_rs.dW) != cid("EGY"))
check("greater Suez exposure predicts a larger loss",
      rho_corr < -0.3)
check("world real income falls under the closure", res_rs.Wworld < 0)
check("Suez traffic falls to zero, Cape and Panama rise",
      route_pct(res_rs)[1] < -99.9 && route_pct(res_rs)[3] > 0 && route_pct(res_rs)[2] > 0)
check("a Hormuz premium concentrates on the Gulf",
      res_hp.dW[cid("QAT")] < -0.05 && abs(res_hp.dW[cid("DEU")]) < 0.05)
check("large traders lose little from any single passage shock",
      all(abs(r.dW[cid("DEU")]) < 0.2 for r in allres))
check("with congestion off, closure cost is proportional to 1/eta^R (within 10%)",
      spread(prod_nc) < 0.10)
check("congestion puts a floor under the cost, damping the eta^R sensitivity",
      spread(prod_c) > spread(prod_nc))

banner("DONE")
println("This variant answers: how much does a chokepoint disruption cost, and who")
println("pays, when nobody owns the passage. It cannot answer what the passage is")
println("worth to its owner - that question needs the toll layer of the full model.")
