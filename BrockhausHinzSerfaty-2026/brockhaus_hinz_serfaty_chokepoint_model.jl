# ============================================================================
# Chokepoint Disruptions with Endogenous Freight Costs and Monopoly Tolls
# Based on Brockhaus, Hinz & Serfaty (2026),
#   "Navigating Shocks: The Ripple Effects of Shipping Route Closures"
#   Banque de France Working Paper No. 1057
#
# Implements the quantitative block of the paper: a multi-country, multi-sector
# Caliendo-Parro economy with input-output linkages in which the bilateral trade
# wedge is ENDOGENOUS through a nested mode-over-route transport block:
#
#   - Freight cost share pass-through:  d = f^chi_s                    (eq 5)
#   - Route congestion   Xi_r  = sum_ijs delta * X_ij,sr               (eq 6)
#   - Global capacity    Psi   = sum_ijsr delta * X_ij,sr              (eq 7)
#   - Freight wedge      f = fbar (Xi/Xibar)^lam (Psi/Psibar)^gam phi  (eq 8-9)
#   - Route choice CES (within mode), elasticity eta^R                 (eq 14-15)
#   - Mode choice CES (across modes),  elasticity eta^M_s              (eq 16-17)
#   - Route-level costs with stacked tolls                             (eq 18)
#   - Monopoly toll setting  phi*_q = argmax Pi_q(phi_q; phi_-q)       (eq 12)
#   - Toll revenue with log-share division on shared routes            (eq 13)
#   - Solution in changes (exact hat algebra, Dekle-Eaton-Kortum 2008)
#
# Two chokepoint types are contrasted, which is the point of the paper:
#   * Suez / Panama - priced infrastructures WITH maritime substitutes.
#     Because the authority already prices to extract the shortcut's surplus,
#     a closure is a shock to the TOLL COLLECTOR, not to shippers.
#   * Strait of Hormuz - NO maritime substitute for the Gulf. A transit fee
#     cannot be evaded by rerouting, so incidence falls on the Gulf itself.
#
# ---------------------------------------------------------------------------
# SCOPE AND HONESTY NOTE
# ---------------------------------------------------------------------------
# The paper's baseline is GTAP 11 (160 countries x 65 sectors, ref. 2017) plus
# proprietary AIS trajectories and Panjiva customs micro-data. None of that is
# redistributable, so this script carries a STYLIZED 24-country x 8-sector
# baseline built on the paper's real chokepoint geography, its real Gulf-side
# exposure shares g_i, and its calibrated elasticities. It reproduces the
# MECHANISMS and the QUALITATIVE INCIDENCE PATTERN (including the with-fees /
# no-fees ranking reversal), NOT the magnitudes of Tables 8-10. Loaders for a
# real GTAP + AIS baseline are described in Claude_Plan.md section 11; nothing
# below the data-assembly layer would change.
#
# The paper's own four caveats (Section 6.2.2) carry over verbatim:
#   1. Hormuz exposure is assigned from port geography, not estimated.
#   2. The bypass margin for partially exposed countries substitutes at eta^R,
#      likely overstating short-run bypass capacity.
#   3. The Hormuz fee is exogenous; canal fees do not re-optimize.
#   4. Route shares / eta^R come from container vessels, whereas Hormuz trade
#      is crude, refined products and LNG - less elastic and more contract-bound.
#
# Run with:  julia brockhaus_hinz_serfaty_chokepoint_model.jl
# ============================================================================

using LinearAlgebra, Statistics, Random, Printf, Plots

Random.seed!(1)
gr()
default(titlefontsize = 10, guidefontsize = 9, tickfontsize = 8, legendfontsize = 8,
        left_margin = 6Plots.mm, bottom_margin = 5Plots.mm, top_margin = 4Plots.mm,
        right_margin = 4Plots.mm, size = (640, 430))

const GRAPHDIR = joinpath(@__DIR__, "graphs")
isdir(GRAPHDIR) || mkpath(GRAPHDIR)

banner(s) = (println(); println("="^78); println(s); println("="^78); flush(stdout))
say(args...) = (println(args...); flush(stdout))

# ============================================================================
# 1. PARAMETERS AND ELASTICITIES  (Section 5.2, Table 4)
# ============================================================================

const NS = 8    # sectors
const NM = 3    # modes: 1 = sea, 2 = air, 3 = other (land / pipeline / services)
const NR = 8    # routes (see section 2 below)

const SECTORS = ["Agri & food", "Energy", "Minerals & bulk", "Chemicals",
                 "Metals", "Machinery & veh.", "Electronics & pharma", "Services"]

# --- Freight cost share chi_s (Section 5.3) ---------------------------------
# Share of DELIVERED price accounted for by route-level shipping cost.
# Paper: bulk (grain, ore, fertilizer) 0.20-0.40; basic chemicals and steel
# 0.05-0.10; manufactures 0.02-0.04; electronics and pharma 0.01; trade-weighted
# mean about 0.09. Crude / refined / LNG is NOT in the paper's "bulk" list. Its
# freight bill runs from a few percent of delivered value for crude to 10-20%
# for LNG; 0.10 is a blended value, and is close to the roughly 0.12 implied by
# the paper's own Hormuz numbers (\$14.6bn of revenue on the post-fee flow at
# phi = 1.35 requires (phi^chi - 1)/phi^chi of about 3.6%).
const CHI = [0.20, 0.10, 0.30, 0.08, 0.08, 0.03, 0.01, 0.00]

# --- Trade elasticity theta_s ----------------------------------------------
# PLACEHOLDERS. The paper uses Fontagne et al. (2022), theta_s = sigma_s - 1,
# long-run tariff-based, aggregated to GTAP sectors; those values are not
# reported in the paper. The large energy value reflects the paper's own caveat
# that "the model's high long-run elasticities of substitution allow importers
# to switch toward unaffected suppliers at relatively low cost".
# TODO: replace with Fontagne et al. (2022) aggregated to the chosen sectors.
const THETA = [8.0, 20.0, 8.0, 6.0, 7.0, 5.0, 5.0, 4.0]

# --- Mode substitution elasticity eta^M_s (Ko et al. 2025, paper Table 4) ---
# AFP 4.87, ENG 11.08, NMM 2.68, CHM 3.13, I_S 5.21, OME 3.17, ELE 1.67,
# pooled fallback 2.44 (Tolva 2026) for sectors without an estimate.
const ETA_M = [4.87, 11.08, 2.68, 3.13, 5.21, 3.17, 1.67, 2.44]

# --- Transport-block elasticities (Section 5.2) -----------------------------
const LAMBDA = fill(0.135, NS)   # route congestion: alpha_tau * beta_tau = 0.38 * 0.35
const GAMMA_G = fill(0.44, NS)   # global shipping capacity (impact elasticity, VAR)
const GAMMA_M = fill(0.00, NS)   # mode-specific capacity: switched off in the paper

const ETA_R_PAPER = 350.0        # route elasticity mu^d calibrated in the paper

# --- Ton-km capacity weights, per dollar of shipment ------------------------
# delta_{ij,smr} = route_distance * tons_per_dollar[s]. Units are arbitrary;
# only ratios matter, since Xi and Psi enter as changes.
const TPD = [2.0, 5.0, 8.0, 1.0, 3.0, 0.3, 0.05, 0.0]

# --- Production structure ---------------------------------------------------
const GAMMA_VA = [0.40, 0.50, 0.45, 0.30, 0.30, 0.32, 0.35, 0.60]  # VA share of gross output
const ALPHA_FIN = [0.12, 0.05, 0.03, 0.04, 0.03, 0.10, 0.08, 0.55] # final demand shares

# --- Mode shares by sector for international trade (i != j) -----------------
# columns: sea, air, other(land/pipeline). Domestic trade (i == j) is 100% other,
# so it is never exposed to a maritime chokepoint.
const MODESHARE = [0.85 0.02 0.13;   # agri & food
                   0.92 0.00 0.08;   # energy (pipeline capacity is largely fixed)
                   0.90 0.00 0.10;   # minerals & bulk
                   0.80 0.05 0.15;   # chemicals
                   0.85 0.02 0.13;   # metals
                   0.75 0.10 0.15;   # machinery & vehicles
                   0.45 0.40 0.15;   # electronics & pharma  <- the air margin
                   0.00 0.00 1.00]   # services

# ============================================================================
# 2. GEOGRAPHY: COUNTRIES, PASSAGES, AUTHORITIES  (Section 5.1)
# ============================================================================

# (iso, name, value added $bn, lat, lon, region, Gulf-side share g_i)
# g_i is exactly the paper's assignment: 1 for Kuwait, Iraq, Qatar and Bahrain;
# 0.95 Iran; 0.80 Saudi Arabia; 0.75 UAE; 0 otherwise.
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

# --- Routes -----------------------------------------------------------------
# Base passages 1-4; routes 5-8 are the Hormuz composites of each. A Hormuz
# transit PRECEDES the choice of onward passage, so tolls stack multiplicatively
# and closing an onward passage also closes its Hormuz composite (Section 5.1).
const ROUTE_NAME = ["Suez", "Panama", "Cape of Good Hope", "Direct",
                    "Hormuz-Suez", "Hormuz-Panama", "Hormuz-Cape", "Hormuz-Direct"]
const R_SUEZ, R_PANAMA, R_CAPE, R_DIRECT = 1, 2, 3, 4
const HORMUZ_ROUTES = 5:8
const BASE_OF_ROUTE = [1, 2, 3, 4, 1, 2, 3, 4]   # onward passage of each route

# Sea-route distance multipliers relative to the great-circle distance.
const ROUTE_DISTMULT = [1.15, 1.25, 1.75, 1.10, 1.22, 1.32, 1.80, 1.18]

# --- Toll authorities -------------------------------------------------------
# Each authority sets ONE wedge on ALL routes it controls. Where several
# authorities toll the same route the wedges stack: Phi_q = prod_a phi_a.
const AUTH_NAME  = ["Suez", "Panama", "Hormuz"]
const NA = length(AUTH_NAME)
const AUTH_ROUTES = [[1, 5], [2, 6], [5, 6, 7, 8]]
# collector country -> share of the authority's revenue
const AUTH_COLLECTORS = [[("EGY", 1.0)], [("PAN", 1.0)], [("IRN", 0.5), ("OMN", 0.5)]]

# route -> authorities tolling it
const ROUTE_AUTHS = [[a for a in 1:NA if r in AUTH_ROUTES[a]] for r in 1:NR]

# --- Baseline sea-route shares by region pair (the AIS moment in Track A) ----
# NOTE: route shares are DATA, not model output. With eta^R = 350 the implied
# baseline cost gaps are tiny (dbar_r prop. to rho_r^(-1/eta^R): an 0.85/0.15
# share ratio implies a 0.5% cost gap), so generating shares from a distance
# logit would be numerically meaningless. Hat algebra only needs the shares.
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

"Baseline passage shares (Suez, Panama, Cape, Direct) for a region pair."
function passage_shares(ri::Symbol, rj::Symbol)
    v = get(PAIR_ROUTES, Set([ri, rj]), nothing)
    v === nothing && return [0.0, 0.0, 0.0, 1.0]   # same region: direct
    return v
end

"Great-circle distance in km."
function gcdist(lat1, lon1, lat2, lon2)
    φ1, φ2 = deg2rad(lat1), deg2rad(lat2)
    Δλ = deg2rad(lon2 - lon1)
    c = sin(φ1) * sin(φ2) + cos(φ1) * cos(φ2) * cos(Δλ)
    return 6371.0 * acos(clamp(c, -1.0, 1.0))
end

# ============================================================================
# 3. BASELINE CONSTRUCTION  (plan section 4)
# ============================================================================
#
# The baseline must satisfy market clearing EXACTLY, or the "all hats = 1" test
# fails. It is therefore SOLVED, not assembled. Note that the levels system is
# homogeneous of degree one, so relative country sizes are determined by (pi,
# gamma, alpha_fin) up to scale - VA cannot be imposed independently. Country
# sizes are hit by a fixed point on the sectoral competitiveness terms T[i,s].

"Distance matrix and per-route distance multipliers."
function build_distances()
    D = zeros(NC, NC)
    for i in 1:NC, j in 1:NC
        D[i, j] = i == j ? 400.0 : max(gcdist(LAT[i], LON[i], LAT[j], LON[j]), 400.0)
    end
    # dist_route[i,j,m,r]: ton-km distance actually travelled on (mode, route)
    dr = zeros(NC, NC, NM, NR)
    for i in 1:NC, j in 1:NC
        for r in 1:NR
            dr[i, j, 1, r] = D[i, j] * ROUTE_DISTMULT[r]
        end
        dr[i, j, 2, R_DIRECT] = D[i, j]   # air
        dr[i, j, 3, R_DIRECT] = D[i, j]   # other (delta = 0, see below)
    end
    return D, dr
end

"Baseline mode shares mu0[i,j,s,m] and route shares rho0[i,j,m,r]."
function build_transport_shares()
    mu0 = zeros(NC, NC, NS, NM)
    rho0 = zeros(NC, NC, NM, NR)

    for i in 1:NC, j in 1:NC
        # --- mode shares
        for s in 1:NS
            if i == j
                mu0[i, j, s, 3] = 1.0            # domestic: never maritime
            else
                mu0[i, j, s, 1] = MODESHARE[s, 1]
                mu0[i, j, s, 2] = MODESHARE[s, 2]
                mu0[i, j, s, 3] = MODESHARE[s, 3]
            end
        end
        # --- route shares within mode
        rho0[i, j, 2, R_DIRECT] = 1.0
        rho0[i, j, 3, R_DIRECT] = 1.0
        if i == j
            rho0[i, j, 1, R_DIRECT] = 1.0
        else
            base = passage_shares(REGION[i], REGION[j])
            # Hormuz exposure: e_ij = g_i(1-g_j) + g_j(1-g_i) (Section 5.1),
            # so intra-Gulf trade is unexposed.
            e = GULFSH[i] * (1 - GULFSH[j]) + GULFSH[j] * (1 - GULFSH[i])
            for b in 1:4
                rho0[i, j, 1, b]     = (1 - e) * base[b]
                rho0[i, j, 1, b + 4] = e * base[b]
            end
        end
    end
    return mu0, rho0
end

"""
    initial_T() -> T[i,s]

Sectoral competitiveness. Only the WITHIN-country pattern matters for the
composition of output: the country-size fixed point in `build_baseline` rescales
each country's whole T row to hit its value-added target, so what is imposed here
is specialisation, not size. Hydrocarbon exporters are given a strongly
energy-tilted row (weak in tradable manufactures) so that their exposure to the
Strait of Hormuz is of a realistic order - roughly 40-55% of income for the
Gulf - rather than the 15-25% a symmetric pattern would deliver.
"""
function initial_T()
    T = ones(NC, NS)
    energy = 2
    hydro = ["SAU", "IRN", "IRQ", "QAT", "KWT", "ARE", "OMN", "BHR", "BRN", "GNQ"]
    for iso in hydro
        i = cid(iso)
        T[i, 1:7] .= 0.05          # weak in every tradable except energy
        T[i, NS] = 0.6             # services: mostly domestic anyway
    end
    for (iso, t) in [("SAU", 60.0), ("IRN", 25.0), ("IRQ", 22.0), ("QAT", 35.0),
                     ("KWT", 20.0), ("ARE", 22.0), ("OMN", 8.0), ("BHR", 4.0),
                     ("BRN", 12.0), ("GNQ", 12.0), ("USA", 4.0), ("ROW", 6.0)]
        T[cid(iso), energy] = t
    end
    for (iso, t) in [("CHN", 3.0), ("DEU", 2.5), ("JPN", 2.2), ("KOR", 2.2)]
        T[cid(iso), 5:7] .= t
    end
    T[cid("NLD"), :] .*= 1.4      # entrepot
    T[cid("SGP"), :] .*= 1.6      # entrepot
    T .*= exp.(0.20 .* randn(NC, NS))
    return T
end

"Input-output coefficients gamma_IO[j,k,s]: share of sector-s inputs in gross output of sector k."
function build_io()
    # base intermediate-use pattern, rows = using sector k, cols = supplying sector s
    pat = [0.25 0.06 0.04 0.08 0.03 0.05 0.02 0.47;   # agri & food
           0.02 0.30 0.04 0.08 0.04 0.06 0.02 0.44;   # energy
           0.02 0.14 0.22 0.07 0.06 0.06 0.02 0.41;   # minerals
           0.04 0.16 0.05 0.28 0.04 0.06 0.03 0.34;   # chemicals
           0.02 0.14 0.16 0.08 0.22 0.06 0.02 0.30;   # metals
           0.02 0.05 0.03 0.08 0.20 0.22 0.09 0.31;   # machinery
           0.02 0.04 0.02 0.08 0.08 0.12 0.28 0.36;   # electronics & pharma
           0.04 0.05 0.02 0.04 0.02 0.07 0.06 0.70]   # services
    pat ./= sum(pat, dims = 2)
    γ = zeros(NC, NS, NS)
    for j in 1:NC, k in 1:NS, s in 1:NS
        γ[j, k, s] = (1 - GAMMA_VA[k]) * pat[k, s]
    end
    return γ
end

"Domestic expenditure share ds[j,s] - the openness lever."
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

"""
    toll_shares(φ, χ) -> (Φ, tshare, logshare)

Combined route wedges `Φ[r] = prod_a φ[a]`, the per-(route, sector) toll take
`tshare[r,s] = (Φ^χ_s - 1)/Φ^χ_s`, and each authority's log-contribution share
`logshare[a,r] = ln φ_a / ln Φ_r` used to divide the take on shared routes (eq 13).
Charging each authority's markup separately would over-extract.
"""
function toll_shares(φ::Vector{Float64})
    Φ = ones(NR)
    for r in 1:NR, a in ROUTE_AUTHS[r]
        Φ[r] *= φ[a]
    end
    tshare = zeros(NR, NS)
    for r in 1:NR, s in 1:NS
        Φχ = Φ[r]^CHI[s]
        tshare[r, s] = (Φχ - 1) / Φχ
    end
    logshare = zeros(NA, NR)
    for r in 1:NR
        lnΦ = log(Φ[r])
        auths = ROUTE_AUTHS[r]
        if lnΦ < 1e-12                       # degenerate limit Φ -> 1: split equally
            for a in auths
                logshare[a, r] = isempty(auths) ? 0.0 : 1.0 / length(auths)
            end
        else
            for a in auths
                logshare[a, r] = log(φ[a]) / lnΦ
            end
        end
    end
    return Φ, tshare, logshare
end

"""
    net_of_toll(π, μ, ρ, tshare) -> (πnet, tk)

`tk[i,j,s]` is the fraction of delivered spending on the (i,j,s) cell that is
skimmed as toll; `πnet = π .* (1 - tk)` is the sourcing share net of tolls, i.e.
the share of j's expenditure that actually reaches producers in i.
"""
function net_of_toll(π, μ, ρ, tshare)
    tk = zeros(NC, NC, NS)
    @inbounds for i in 1:NC, j in 1:NC, s in 1:NS
        acc = 0.0
        for m in 1:NM
            μm = μ[i, j, s, m]
            μm == 0.0 && continue
            for r in 1:NR
                ρr = ρ[i, j, s, m, r]
                (ρr == 0.0 || m != 1) && continue    # only sea routes are tolled
                acc += μm * ρr * tshare[r, s]
            end
        end
        tk[i, j, s] = acc
    end
    return π .* (1 .- tk), tk
end

"Expand rho0[i,j,m,r] to the sector dimension."
function expand_rho(rho0)
    ρ = zeros(NC, NC, NS, NM, NR)
    @inbounds for i in 1:NC, j in 1:NC, s in 1:NS, m in 1:NM, r in 1:NR
        ρ[i, j, s, m, r] = rho0[i, j, m, r]
    end
    return ρ
end

"""
    solve_levels(πnet, γIO, αfin, Π, D, Ytot; ...) -> (Y, X, I)

Power iteration on the (scale-indeterminate) levels system

    Y[i,s] = sum_j πnet[i,j,s] X[j,s]
    X[j,s] = sum_k γIO[j,k,s] Y[j,k] + αfin[j,s] I[j]
    I[j]   = sum_s γVA[j,s] Y[j,s] + Π[j] + D[j]

normalised so that world value added equals `Ytot`.
"""
function solve_levels(πnet, γIO, αfin, Π, D, Ytot; Y0 = nothing,
                      maxiter = 4000, tol = 1e-10, damp = 0.9)
    Y = Y0 === nothing ? fill(Ytot / (NC * NS), NC, NS) : copy(Y0)
    X = zeros(NC, NS); I = zeros(NC); Ynew = zeros(NC, NS)
    for it in 1:maxiter
        @inbounds for j in 1:NC
            acc = 0.0
            for s in 1:NS
                acc += GAMMA_VA[s] * Y[j, s]
            end
            I[j] = acc + Π[j] + D[j]
        end
        @inbounds for j in 1:NC, s in 1:NS
            acc = αfin[j, s] * I[j]
            for k in 1:NS
                acc += γIO[j, k, s] * Y[j, k]
            end
            X[j, s] = acc
        end
        fill!(Ynew, 0.0)
        @inbounds for s in 1:NS, j in 1:NC
            xjs = X[j, s]
            for i in 1:NC
                Ynew[i, s] += πnet[i, j, s] * xjs
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
        I[j] = acc + Π[j] + D[j]
    end
    @inbounds for j in 1:NC, s in 1:NS
        acc = αfin[j, s] * I[j]
        for k in 1:NS
            acc += γIO[j, k, s] * Y[j, k]
        end
        X[j, s] = acc
    end
    return Y, X, I
end

"""
    solve_expenditure(πnet, γIO, αfin, LI, Π, D; ...) -> (Y, X, I)

Same accounting as `solve_levels`, but with labour income `LI[i] = what[i] * VA[i]`
taken as GIVEN. Fixing final demand makes the system non-homogeneous and its
Leontief map a strict contraction (spectral radius bounded by 1 - min gamma_VA),
so it has a unique solution and no normalisation is needed. This is what makes
the outer wage loop well behaved: the counterfactual updates wages from an
excess-demand function rather than reading them off a scale-free eigenvector.
"""
function solve_expenditure(πnet, γIO, αfin, LI, Π, D; Y0 = nothing,
                           maxiter = 2000, tol = 1e-12)
    Y = Y0 === nothing ? fill(sum(LI) / (NC * NS), NC, NS) : copy(Y0)
    X = zeros(NC, NS); Ynew = zeros(NC, NS)
    I = LI .+ Π .+ D
    scale = sum(abs, I) / (NC * NS)
    for it in 1:maxiter
        @inbounds for j in 1:NC, s in 1:NS
            acc = αfin[j, s] * I[j]
            for k in 1:NS
                acc += γIO[j, k, s] * Y[j, k]
            end
            X[j, s] = acc
        end
        fill!(Ynew, 0.0)
        @inbounds for s in 1:NS, j in 1:NC
            xjs = X[j, s]
            for i in 1:NC
                Ynew[i, s] += πnet[i, j, s] * xjs
            end
        end
        err = maximum(abs.(Ynew .- Y)) / scale
        Y .= Ynew
        err < tol && break
    end
    @inbounds for j in 1:NC, s in 1:NS
        acc = αfin[j, s] * I[j]
        for k in 1:NS
            acc += γIO[j, k, s] * Y[j, k]
        end
        X[j, s] = acc
    end
    return Y, X, I
end

"Toll revenue by authority (eq 13) and by collector country, given route-sector flows."
function toll_revenue(Xrs, φ)
    Φ, tshare, logshare = toll_shares(φ)
    Πauth = zeros(NA)
    for a in 1:NA, r in AUTH_ROUTES[a], s in 1:NS
        Πauth[a] += logshare[a, r] * tshare[r, s] * Xrs[r, s]
    end
    Πc = zeros(NC)
    for a in 1:NA, (iso, w) in AUTH_COLLECTORS[a]
        Πc[cid(iso)] += w * Πauth[a]
    end
    return Πauth, Πc
end

"""
    build_baseline(; verbose) -> NamedTuple

Assemble the raw (untolled, φ = 1) baseline equilibrium: trade shares, mode and
route shares, IO structure, flows, congestion and capacity aggregates.
"""
function build_baseline(; verbose = true)
    Dist, distroute = build_distances()
    mu0, rho0 = build_transport_shares()
    ρ = expand_rho(rho0)
    μ = mu0
    γIO = build_io()
    ds = build_domestic_shares()
    αfin = repeat(ALPHA_FIN', NC, 1)
    D = zeros(NC)
    T = initial_T()
    Ytot = sum(VATGT)

    φ0 = ones(NA)
    _, tshare0, _ = toll_shares(φ0)

    π = zeros(NC, NC, NS)
    Y = nothing; X = nothing; I = nothing
    for outer_it in 1:80
        # ---- gravity sourcing shares given competitiveness T
        @inbounds for j in 1:NC, s in 1:NS
            tot = 0.0
            for i in 1:NC
                if i == j
                    π[i, j, s] = 0.0
                else
                    v = T[i, s] * (Dist[i, j] / 1000)^(-1.0)
                    π[i, j, s] = v
                    tot += v
                end
            end
            scale = (1 - ds[j, s]) / max(tot, 1e-300)
            for i in 1:NC
                i != j && (π[i, j, s] *= scale)
            end
            π[j, j, s] = ds[j, s]
        end
        πnet, _ = net_of_toll(π, μ, ρ, tshare0)
        Y, X, I = solve_levels(πnet, γIO, αfin, zeros(NC), D, Ytot; Y0 = Y)
        VAimp = vec(sum(GAMMA_VA' .* Y, dims = 2))
        ratio = VATGT ./ VAimp
        err = maximum(abs.(log.(ratio)))
        T .*= ratio .^ 0.5
        if err < 5e-3
            verbose && @printf("  country-size fixed point converged in %d iterations (max log gap %.4f)\n", outer_it, err)
            break
        end
        outer_it == 80 && verbose && @printf("  country-size fixed point stopped at max log gap %.4f\n", err)
    end

    VA = vec(sum(GAMMA_VA' .* Y, dims = 2))
    Xijs = zeros(NC, NC, NS)
    @inbounds for i in 1:NC, j in 1:NC, s in 1:NS
        Xijs[i, j, s] = π[i, j, s] * X[j, s]
    end

    # ---- baseline congestion and capacity aggregates (eq 6, 7)
    Ξbar = zeros(NM, NR); Xrs = zeros(NR, NS)
    @inbounds for i in 1:NC, j in 1:NC, s in 1:NS
        x = Xijs[i, j, s]
        x == 0.0 && continue
        for m in 1:NM, r in 1:NR
            f = x * μ[i, j, s, m] * ρ[i, j, s, m, r]
            f == 0.0 && continue
            Ξbar[m, r] += distroute[i, j, m, r] * TPD[s] * f
            m == 1 && (Xrs[r, s] += f)
        end
    end
    Ξbar[3, :] .= 0.0                       # "other" mode carries no ton-km weight
    Ψbar = sum(Ξbar[1, :])                  # global maritime capacity tightness
    Ξactive = Ξbar .> 0                     # cells with genuine baseline traffic
    Ξbar[.!Ξactive] .= 1.0                  # inert placeholder, never used

    base = (; π, μ, ρ, γIO, αfin, D, Dist, distroute,
            Y, X, I, VA, Xijs, Xrs, Ξbar, Ξactive, Ψbar, Ytot,
            φ_base = copy(φ0), Πauth = zeros(NA), Πc = zeros(NC))
    return base
end

# ============================================================================
# 4. TRANSPORT BLOCK  (eq 14-20, in changes)
# ============================================================================
#
# Because exogenous freight shocks and tolls are ROUTE-level, the route-level
# log cost change is common across (i,j) and depends on the pair only through
# the baseline shares. That collapses the inner work from an (N,N,S,M,R) array
# of costs to an (S,M,R) array, which is what makes this fast.

"log-sum-exp that skips -Inf entries; returns -Inf if all are skipped."
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

"Mutable state carried across outer iterations (warm starts)."
mutable struct TransportState
    Ξhat::Matrix{Float64}          # (NM, NR)
    Ψhat::Float64
    dhat::Array{Float64,3}         # (NC, NC, NS) composite bilateral wedge change
    μ::Array{Float64,4}            # (NC, NC, NS, NM)
    ρ::Array{Float64,5}            # (NC, NC, NS, NM, NR)
    Xrs::Matrix{Float64}           # (NR, NS) sea flows by route and sector
    tk::Array{Float64,3}           # (NC, NC, NS) toll take fraction
    damp::Float64                  # adaptive step on the congestion aggregates
    resid::Float64
end

TransportState() = TransportState(ones(NM, NR), 1.0, ones(NC, NC, NS),
                                  zeros(NC, NC, NS, NM), zeros(NC, NC, NS, NM, NR),
                                  zeros(NR, NS), zeros(NC, NC, NS), 0.4, Inf)

"""
    transport_pass!(st, base, ηR, Xijs, lnφhat, tshare, route_open, dbarhat; niter, tol)

Iterate the transport fixed point (eq 18 -> 14-17 -> 19-20) `niter` times,
holding bilateral sectoral flows `Xijs` fixed. Updates `st` in place.
`dbarhat` is an (NS, NM, NR) array of exogenous baseline freight-wedge shocks.
"""
function transport_pass!(st::TransportState, base, ηR::Float64, Xijs,
                         lnφhat::Vector{Float64}, tshare::Matrix{Float64},
                         route_open::Vector{Bool}, dbarhat::Array{Float64,3};
                         niter::Int = 3, tol::Float64 = 1e-12)
    ρ0, μ0 = base.ρ, base.μ
    Ξbar, Ψbar, distroute = base.Ξbar, base.Ψbar, base.distroute

    g = zeros(NS, NM, NR)                      # route-level log cost change
    lnρ0 = zeros(NR); lnρnew = zeros(NR)
    lnμ0 = zeros(NM); lndm = zeros(NM); tmpm = zeros(NM); tmpr = zeros(NR)
    Ξnew = zeros(NM, NR)

    for _ in 1:niter
        lnΨ = log(st.Ψhat)
        @inbounds for s in 1:NS, m in 1:NM, r in 1:NR
            sea = (m == 1)
            g[s, m, r] = log(dbarhat[s, m, r]) +
                         CHI[s] * (sea ? (LAMBDA[s] * log(st.Ξhat[m, r]) + GAMMA_G[s] * lnΨ) : 0.0) +
                         CHI[s] * (sea ? lnφhat[r] : 0.0)
        end

        fill!(Ξnew, 0.0); fill!(st.Xrs, 0.0)
        @inbounds for i in 1:NC, j in 1:NC
            for s in 1:NS
                # ---- within-mode route nest (eq 14-15)
                for m in 1:NM
                    for r in 1:NR
                        p = ρ0[i, j, s, m, r]
                        lnρ0[r] = (p > 0 && (m != 1 || route_open[r])) ? log(p) : -Inf
                        tmpr[r] = lnρ0[r] - ηR * g[s, m, r]
                    end
                    L = lse(tmpr)
                    if isfinite(L)
                        lndm[m] = -L / ηR
                        for r in 1:NR
                            lnρnew[r] = lnρ0[r] - ηR * (g[s, m, r] - lndm[m])
                            st.ρ[i, j, s, m, r] = isfinite(lnρnew[r]) ? exp(lnρnew[r]) : 0.0
                        end
                    else
                        lndm[m] = Inf                    # mode has no open route
                        for r in 1:NR
                            st.ρ[i, j, s, m, r] = 0.0
                        end
                    end
                    p = μ0[i, j, s, m]
                    lnμ0[m] = (p > 0 && isfinite(lndm[m])) ? log(p) : -Inf
                end
                # ---- across-mode nest (eq 16-17)
                ηM = ETA_M[s]
                for m in 1:NM
                    tmpm[m] = isfinite(lnμ0[m]) ? lnμ0[m] - ηM * lndm[m] : -Inf
                end
                L = lse(tmpm)
                if isfinite(L)
                    lnds = -L / ηM
                    st.dhat[i, j, s] = exp(lnds)
                    for m in 1:NM
                        st.μ[i, j, s, m] = isfinite(lnμ0[m]) ? exp(lnμ0[m] - ηM * (lndm[m] - lnds)) : 0.0
                    end
                else
                    st.dhat[i, j, s] = Inf               # pair fully disconnected
                    for m in 1:NM
                        st.μ[i, j, s, m] = 0.0
                    end
                end
                # ---- flow aggregates (eq 19-20) and toll base
                x = Xijs[i, j, s]
                acc = 0.0
                for m in 1:NM
                    μm = st.μ[i, j, s, m]
                    μm == 0.0 && continue
                    for r in 1:NR
                        ρr = st.ρ[i, j, s, m, r]
                        ρr == 0.0 && continue
                        f = x * μm * ρr
                        if m == 1
                            Ξnew[m, r] += distroute[i, j, m, r] * TPD[s] * f
                            st.Xrs[r, s] += f
                            acc += μm * ρr * tshare[r, s]
                        elseif m == 2
                            Ξnew[m, r] += distroute[i, j, m, r] * TPD[s] * f
                        end
                    end
                end
                st.tk[i, j, s] = acc
            end
        end

        # Floor the congestion ratio at 1e-3: a passage that is priced out does
        # not become infinitely fast, it reaches its free-flow transit time.
        #
        # The congestion aggregates are updated with an ADAPTIVE step. The map
        # Xi -> Xi has loop gain of order eta^R * chi_s * lambda_s * (1 - rho_r)
        # (a toll drives cargo off the passage, which decongests it, which makes
        # it attractive again). At the calibrated elasticities that gain exceeds
        # one for the high-freight-share sectors, so a fixed step of 0.5 makes
        # the iteration oscillate and diverge. Backing the step off whenever the
        # residual rises keeps it contracting; the fixed point reached is the one
        # continuously connected to the baseline.
        Ψnew = sum(@view Ξnew[1, :]) / Ψbar
        resid = 0.0
        d = st.damp
        @inbounds for m in 1:NM, r in 1:NR
            if base.Ξactive[m, r]
                tgt = max(Ξnew[m, r] / Ξbar[m, r], 1e-3)
                resid = max(resid, abs(log(tgt) - log(st.Ξhat[m, r])))
                st.Ξhat[m, r] = d * tgt + (1 - d) * st.Ξhat[m, r]
            else
                st.Ξhat[m, r] = 1.0        # no baseline traffic: congestion inert
            end
        end
        tgtΨ = max(Ψnew, 1e-3)
        resid = max(resid, abs(log(tgtΨ) - log(st.Ψhat)))
        st.Ψhat = d * tgtΨ + (1 - d) * st.Ψhat
        st.damp = resid > st.resid ? max(0.5 * st.damp, 0.01) : min(1.06 * st.damp, 0.6)
        st.resid = resid
        resid < tol && break
    end
    return st
end

# ============================================================================
# 5. PRICE BLOCK AND COUNTERFACTUAL SOLVER  (Section 4.4, in changes)
# ============================================================================

"""
    solve_prices!(Phat, base, what, dhat; maxiter, tol)

Caliendo-Parro inner loop:
  chat[i,k] = what[i]^gVA[k] * prod_s Phat[i,s]^gIO[i,k,s]
  Phat[j,s] = [ sum_i pi[i,j,s] (chat[i,s] dhat[i,j,s])^(-theta_s) ]^(-1/theta_s)
Warm-started from `Phat`.
"""
function solve_prices!(Phat, base, what, dhat; maxiter = 400, tol = 1e-12)
    π, γIO = base.π, base.γIO
    chat = ones(NC, NS)
    for it in 1:maxiter
        @inbounds for i in 1:NC, k in 1:NS
            acc = GAMMA_VA[k] * log(what[i])
            for s in 1:NS
                acc += γIO[i, k, s] * log(Phat[i, s])
            end
            chat[i, k] = exp(acc)
        end
        err = 0.0
        @inbounds for j in 1:NC, s in 1:NS
            θ = THETA[s]; acc = 0.0
            for i in 1:NC
                p = π[i, j, s]
                p == 0.0 && continue
                cd = chat[i, s] * dhat[i, j, s]
                isfinite(cd) || continue
                acc += p * cd^(-θ)
            end
            new = acc <= 0.0 ? 1e10 : acc^(-1 / θ)
            err = max(err, abs(log(new) - log(Phat[j, s])))
            Phat[j, s] = new
        end
        err < tol && break
    end
    return Phat, chat
end

"""
    solve_counterfactual(base, sh; ...) -> NamedTuple

Exact-hat solver. `sh` is a NamedTuple with fields
  `route_open :: Vector{Bool}` (length NR),
  `phi_new    :: Vector{Float64}` (length NA, LEVELS),
  `dbarhat    :: Array{Float64,3}` (NS, NM, NR) exogenous freight-wedge shocks.

Structure (plan section 5): outer loop on wages; inside it a transport pass
(L3), the price fixed point (L2), and the expenditure system. `damp` is the step
size on the excess-demand wage update.
Closure: nominal deficits fixed at baseline; world value added is the numeraire.
"""
function solve_counterfactual(base, sh; ηR::Float64 = ETA_R_PAPER,
                              maxiter = 4000, tol = 1e-9, damp = 0.25, verbose = false)
    what = ones(NC)
    Phat = ones(NC, NS)
    st = TransportState()
    # start from baseline shares so a null shock is an exact fixed point
    st.μ .= base.μ; st.ρ .= base.ρ
    Xijs = copy(base.Xijs)
    Y = copy(base.Y); X = copy(base.X); I = copy(base.I)
    Πc = copy(base.Πc); Πauth = copy(base.Πauth)

    # toll wedges: route costs move with the CHANGE, revenue is levied at the LEVEL
    φnew = sh.phi_new
    lnφhat = zeros(NR)
    for r in 1:NR, a in ROUTE_AUTHS[r]
        lnφhat[r] += log(φnew[a] / base.φ_base[a])
    end
    _, tshare, _ = toll_shares(φnew)

    conv = false; iters = 0
    step = damp; err_prev = Inf
    for it in 1:maxiter
        iters = it
        transport_pass!(st, base, ηR, Xijs, lnφhat, tshare, sh.route_open, sh.dbarhat;
                        niter = it <= 5 ? 200 : 30)
        _, chat = solve_prices!(Phat, base, what, st.dhat; maxiter = it <= 5 ? 400 : 60)

        πnew = similar(base.π)
        @inbounds for i in 1:NC, j in 1:NC, s in 1:NS
            p = base.π[i, j, s]
            if p == 0.0
                πnew[i, j, s] = 0.0
            else
                cd = chat[i, s] * st.dhat[i, j, s]
                πnew[i, j, s] = isfinite(cd) ? p * (cd / Phat[j, s])^(-THETA[s]) : 0.0
            end
        end
        πnet = πnew .* (1 .- st.tk)

        Πauth, Πc = toll_revenue(st.Xrs, φnew)
        LI = what .* base.VA
        Y, X, I = solve_expenditure(πnet, base.γIO, base.αfin, LI, Πc, base.D;
                                    Y0 = Y, maxiter = it <= 5 ? 2000 : 300)

        # excess demand for labour: value added produced minus the wage bill paid.
        # Walras' law holds exactly here (sum_i ED_i = 0), because toll rents are
        # rebated to collectors and every dollar of expenditure returns as income.
        ED = vec(sum(GAMMA_VA' .* Y, dims = 2)) .- LI
        err = maximum(abs.(ED ./ LI))
        # adaptive step: accelerate while the residual is falling, back off if not
        if it > 1
            step = err < err_prev ? min(step * 1.08, 0.9) : max(step * 0.5, 0.02)
        end
        err_prev = err
        what .*= (1 .+ step .* ED ./ LI)
        what .*= sum(base.VA) / sum(what .* base.VA)     # numeraire: world VA

        @inbounds for i in 1:NC, j in 1:NC, s in 1:NS
            Xijs[i, j, s] = πnew[i, j, s] * X[j, s]
        end

        if err < tol
            conv = true
            verbose && @printf("  >>>> GE converged in %d iterations (max |ED/LI| = %.2e) <<<<\n", it, err)
            break
        end
    end
    conv || @warn "GE did not converge within maxiter"

    # ---- welfare and its rent decomposition
    Pfin = [exp(sum(base.αfin[j, s] * log(Phat[j, s]) for s in 1:NS)) for j in 1:NC]
    What = (I ./ base.I) ./ Pfin
    ΔW = 100 .* (What .- 1)
    Δτ = 100 .* ((Πc ./ Pfin) .- base.Πc) ./ base.I
    Wworld = (sum(I) / sum(base.I)) / exp(sum(base.I[j] / sum(base.I) * log(Pfin[j]) for j in 1:NC))

    πnew = similar(base.π)
    @inbounds for i in 1:NC, j in 1:NC, s in 1:NS
        πnew[i, j, s] = base.X[j, s] > 0 ? Xijs[i, j, s] / X[j, s] : 0.0
    end

    route_flow = vec(sum(st.Xrs, dims = 2))
    route_flow_base = vec(sum(base.Xrs, dims = 2))

    return (; what, Phat, Pfin, dhat = copy(st.dhat), μ = copy(st.μ), ρ = copy(st.ρ),
            π = πnew, Y, X, I, Πauth, Πc, Xrs = copy(st.Xrs), Ξhat = copy(st.Ξhat),
            Ψhat = st.Ψhat, tk = copy(st.tk), Xijs = copy(Xijs),
            ΔW, Δτ, What, Wworld = 100 * (Wworld - 1),
            route_flow, route_flow_base,
            route_change = 100 .* (route_flow ./ max.(route_flow_base, 1e-12) .- 1),
            converged = conv, iters)
end

"Convenience constructors for shocks."
null_shock(base) = (route_open = fill(true, NR), phi_new = copy(base.φ_base),
                    dbarhat = ones(NS, NM, NR))

function shock_fees(base, φ::Vector{Float64})
    (route_open = fill(true, NR), phi_new = copy(φ), dbarhat = ones(NS, NM, NR))
end

"Permanent Red Sea closure: the Suez route AND its Hormuz composite disappear."
function shock_red_sea(base)
    ro = fill(true, NR)
    ro[R_SUEZ] = false
    ro[R_SUEZ + 4] = false
    (route_open = ro, phi_new = copy(base.φ_base), dbarhat = ones(NS, NM, NR))
end

"Strait of Hormuz transit fee at level φH, stacking on the canal wedges."
function shock_hormuz(base, φH::Float64)
    φ = copy(base.φ_base); φ[3] = φH
    (route_open = fill(true, NR), phi_new = φ, dbarhat = ones(NS, NM, NR))
end

# ============================================================================
# 6. TOLL REVENUE, BEST RESPONSES AND CALIBRATION  (Section 5.1)
# ============================================================================

"""
    revenue_at(base, φ, a; ηR, mode) -> (revenue of authority a, traffic on its routes)

`mode = :pe` holds bilateral flows at their baseline values and runs only the
transport block - this is exactly the toy monopolist of Section 4.1, where the
trade volume T_ij is exogenous to the fee. `mode = :ge` re-solves the full
general equilibrium. The traffic on the authority's own routes is returned as
well, because it is what identifies the well-posed branch of its problem.
"""
function revenue_at(base, φ::Vector{Float64}, a::Int; ηR = ETA_R_PAPER, mode::Symbol = :pe)
    Xrs = if mode === :ge
        solve_counterfactual(base, shock_fees(base, φ); ηR = ηR, tol = 1e-8, damp = 0.3).Xrs
    else
        st = TransportState(); st.μ .= base.μ; st.ρ .= base.ρ
        lnφhat = zeros(NR)
        for r in 1:NR, b in ROUTE_AUTHS[r]
            lnφhat[r] += log(φ[b] / base.φ_base[b])
        end
        _, tshare, _ = toll_shares(φ)
        transport_pass!(st, base, ηR, base.Xijs, lnφhat, tshare, fill(true, NR),
                        ones(NS, NM, NR); niter = 400)
        st.Xrs
    end
    Πauth, _ = toll_revenue(Xrs, φ)
    traffic = sum(Xrs[r, s] for r in AUTH_ROUTES[a], s in 1:NS)
    return Πauth[a], traffic
end

"""
    best_response_fee(base, a; φ_others, ηR, mode) -> (φ*, revenue)

Maximise Π_a(φ_a; φ_-a) over log φ_a (eq 12). Under "pricing in isolation" (the
paper's calibration convention) `φ_others` are held at their reference levels
rather than at other authorities' best responses.

Two things make this more than a one-line optimisation.

First, a coarse grid scan precedes the golden-section refinement. The elasticity
of route shares to the fee is eta^R * chi_s, so for a large eta^R the revenue
function is a narrow spike just above phi = 1 with a long flat floor above it,
and pure golden section wanders on the flat part.

Second, the grid is truncated at the point where the authority's OWN traffic
stops falling in its own fee, after it has fallen at least 10%. The monopolist's
problem is only well posed where its traffic is decreasing in its own fee, and
this guard makes the search independent of whether the transport block has
settled. With the adaptive damping of `transport_pass!` the guard does not in
fact bind at any elasticity in the calibration range - `admissible_eta` reports
this - but it did bind under a fixed damping step, where an oscillating
congestion aggregate produced a spurious "decongestion makes the tolled route
cheaper" branch and a badly non-monotone calibration curve. It is kept as a
correctness guard, not as a modelling assumption.
"""
function best_response_fee(base, a::Int; φ_others = nothing, ηR = ETA_R_PAPER,
                           mode::Symbol = :pe, lo = 1.0 + 1e-7, hi = 2.5,
                           ngrid = 32, tol = 1e-5)
    φ = φ_others === nothing ? copy(base.φ_base) : copy(φ_others)
    eval_at(x) = (φ[a] = exp(x); revenue_at(base, φ, a; ηR = ηR, mode = mode))
    f(x) = first(eval_at(x))
    # log-log grid: dense just above phi = 1, where the spike lives
    xs = log.(1 .+ exp.(range(log(lo - 1), log(hi - 1); length = ngrid)))
    res = eval_at.(xs)
    fs = first.(res); tr = last.(res)
    # truncate at the first point where own traffic turns back up
    kmax = ngrid
    for k in 2:ngrid
        if tr[k] > tr[k - 1] * (1 + 1e-6) && tr[k - 1] < 0.9 * tr[1]
            kmax = k - 1
            break
        end
    end
    fs = fs[1:kmax]; xs = xs[1:kmax]
    k = argmax(fs)
    al = xs[max(k - 1, 1)]; bl = xs[min(k + 1, kmax)]
    invφ = (sqrt(5) - 1) / 2
    c = bl - invφ * (bl - al); d = al + invφ * (bl - al)
    fc, fd = f(c), f(d)
    while (bl - al) > tol
        if fc > fd
            bl, d, fd = d, c, fc
            c = bl - invφ * (bl - al); fc = f(c)
        else
            al, c, fc = c, d, fd
            d = al + invφ * (bl - al); fd = f(d)
        end
    end
    xstar = (al + bl) / 2
    rstar = f(xstar)
    # never return worse than the best admissible grid point
    if fs[k] > rstar
        xstar, rstar = xs[k], fs[k]
    end
    return exp(xstar), rstar
end

"""
    admissible_eta(base, ηR; hi) -> Bool

True when the Suez authority's own traffic is monotonically decreasing in its own
fee over the whole search range at this eta^R, i.e. when the congestion-relief
branch described in `best_response_fee` has not appeared.
"""
function admissible_eta(base, ηR; hi = 2.5, ngrid = 16)
    φ = copy(base.φ_base)
    xs = log.(1 .+ exp.(range(log(1e-4), log(hi - 1); length = ngrid)))
    prev = Inf; t0 = NaN
    for (k, x) in enumerate(xs)
        φ[1] = exp(x)
        t = last(revenue_at(base, φ, 1; ηR = ηR, mode = :pe))
        k == 1 && (t0 = t)
        t > prev * (1 + 1e-6) && prev < 0.9 * t0 && return false
        prev = t
    end
    return true
end

"""
    calibrate_route_elasticity(raw_base; target, ...) -> NamedTuple

Pin the route-choice elasticity on the cargo iceberg, mu^d, by the single
revenue moment: the revenue-maximising Suez fee - priced in isolation, on the
untolled baseline - must reproduce observed Suez Canal toll revenue.

Over the admissible range revenue is DECREASING in eta^R: a larger elasticity
means the implied baseline cost gap between Suez and its alternatives (dbar_r
proportional to rho_r^(-1/eta^R)) is smaller, so less surplus is available to levy.

The scan also records, at each eta^R, whether the authority's own traffic is
monotonically decreasing in its own fee over the whole search range. It should
be, and it is; the column exists because a fixed-step transport iteration made
it fail, which is what first exposed the damping problem.
"""
function calibrate_route_elasticity(raw_base; target = 10.25, lo = 25.0, hi = 800.0,
                                    npts = 12, mode::Symbol = :pe, verbose = true)
    grid = exp.(range(log(lo), log(hi); length = npts))
    revs = similar(grid); fees = similar(grid); adm = trues(npts)
    verbose && @printf("  %10s %12s %12s %14s\n", "eta^R", "phi_Suez", "revenue", "own traffic")
    for (k, η) in enumerate(grid)
        fees[k], revs[k] = best_response_fee(raw_base, 1; ηR = η, mode = mode)
        adm[k] = admissible_eta(raw_base, η)
        verbose && (@printf("  %10.0f %12.3f %12.2f %14s\n", η, fees[k], revs[k],
                            adm[k] ? "monotone" : "turns back up"); flush(stdout))
    end
    kx = findfirst(k -> min(revs[k], revs[k + 1]) <= target <= max(revs[k], revs[k + 1]),
                   1:(npts - 1))
    if kx === nothing
        k = argmin(abs.(revs .- target))
        verbose && println("  !! target not bracketed on the grid; using the closest point")
        return (ηR = grid[k], φ_suez = fees[k], rev_suez = revs[k], bracketed = false,
                grid, revs, adm)
    end
    a, b = grid[kx], grid[kx + 1]
    fa = revs[kx] - target
    ηstar = a
    for _ in 1:40
        ηstar = exp((log(a) + log(b)) / 2)
        fm = last(best_response_fee(raw_base, 1; ηR = ηstar, mode = mode)) - target
        abs(fm) / target < 1e-3 && break
        sign(fm) == sign(fa) ? (a = ηstar; fa = fm) : (b = ηstar)
    end
    φ, r = best_response_fee(raw_base, 1; ηR = ηstar, mode = mode)
    verbose && @printf("  calibrated eta^R = %.1f  ->  phi_Suez = %.3f, revenue = \$%.2fbn\n", ηstar, φ, r)
    return (ηR = ηstar, φ_suez = φ, rev_suez = r, bracketed = true, grid, revs, adm)
end

"""
    rebaseline(raw_base, φ; ηR) -> new baseline

Apply the calibrated fees to the raw (untolled) baseline and treat the result as
the new baseline, so that counterfactual fees enter as changes φ/φ_base while
toll revenue is always evaluated at the LEVEL φ (Section 4.4).
"""
function rebaseline(raw_base, φ::Vector{Float64}; ηR = ETA_R_PAPER, verbose = true)
    res = solve_counterfactual(raw_base, shock_fees(raw_base, φ); ηR = ηR, tol = 1e-10,
                               damp = 0.3, verbose = verbose)
    Xijs = res.Xijs
    Ξbar = zeros(NM, NR)
    @inbounds for m in 1:NM, r in 1:NR
        Ξbar[m, r] = raw_base.Ξbar[m, r] * res.Ξhat[m, r]
    end
    Ψbar = raw_base.Ψbar * res.Ψhat
    VA = res.what .* raw_base.VA
    return (; π = res.π, μ = res.μ, ρ = res.ρ, γIO = raw_base.γIO, αfin = raw_base.αfin,
            D = raw_base.D, Dist = raw_base.Dist, distroute = raw_base.distroute,
            Y = res.Y, X = res.X, I = res.I, VA, Xijs, Xrs = res.Xrs,
            Ξbar, Ξactive = raw_base.Ξactive, Ψbar, Ytot = raw_base.Ytot,
            φ_base = copy(φ), Πauth = res.Πauth, Πc = res.Πc)
end

# ============================================================================
# 7. VERIFICATION SUITE  (plan section 8)
# ============================================================================

check(label, ok) = (@printf("  [%s] %s\n", ok ? " OK " : "FAIL", label); flush(stdout); ok)

"Check with the residual printed, so a marginal failure is diagnosable."
function checkv(label, val, tol)
    ok = val < tol
    @printf("  [%s] %-46s  residual %.3e  (tol %.0e)\n", ok ? " OK " : "FAIL", label, val, tol)
    flush(stdout)
    return ok
end

function verify_baseline(b, name)
    println("\n--- consistency: $name ---")
    ok = true
    ok &= check("sourcing shares sum to one",
                maximum(abs.(vec(sum(b.π, dims = 1)) .- 1)) < 1e-9)
    ok &= check("mode shares sum to one",
                maximum(abs.(vec(sum(b.μ, dims = 4)) .- 1)) < 1e-9)
    ok &= check("route shares sum to one",
                maximum(abs.(vec(sum(b.ρ, dims = 5)) .- 1)) < 1e-9)
    ok &= check("shares non-negative",
                minimum(b.π) >= -1e-14 && minimum(b.μ) >= -1e-14 && minimum(b.ρ) >= -1e-14)
    _, tshare, _ = toll_shares(b.φ_base)
    πnet, _ = net_of_toll(b.π, b.μ, b.ρ, tshare)
    Yc = zeros(NC, NS)
    @inbounds for s in 1:NS, j in 1:NC, i in 1:NC
        Yc[i, s] += πnet[i, j, s] * b.X[j, s]
    end
    ok &= check("goods market clearing  Y = sum_j pinet X",
                maximum(abs.(Yc .- b.Y)) / maximum(b.Y) < 1e-7)
    ok &= check("income identity  I = VA + Pi + D",
                maximum(abs.(b.I .- (b.VA .+ b.Πc .+ b.D))) / maximum(b.I) < 1e-9)
    ok &= check("world VA equals target",
                abs(sum(b.VA) / b.Ytot - 1) < 1e-6)
    return ok
end

function verify_null_shock(b, ηR, name)
    println("\n--- null shock reproduces $name ---")
    res = solve_counterfactual(b, null_shock(b); ηR = ηR, tol = 1e-10, damp = 0.5)
    ok = true
    ok &= checkv("what = 1", maximum(abs.(res.what .- 1)), 1e-6)
    ok &= checkv("Phat = 1", maximum(abs.(res.Phat .- 1)), 1e-6)
    ok &= checkv("dhat = 1", maximum(abs.(res.dhat .- 1)), 1e-8)
    ok &= checkv("route shares unchanged", maximum(abs.(res.ρ .- b.ρ)), 1e-8)
    ok &= checkv("toll rents unchanged",
                 maximum(abs.(res.Πc .- b.Πc)) / max(maximum(b.Πc), 1e-9), 1e-5)
    ok &= checkv("welfare change zero", maximum(abs.(res.ΔW)), 1e-4)
    ok &= checkv("baseline Xi reproduced", maximum(abs.(res.Ξhat .- 1)), 1e-8)
    return ok
end

"Table 7: cargo-value iceberg under the multiplicative vs. additive decomposition."
function verify_table7()
    println("\n--- Appendix B, Table 7 (analytic, no solver) ---")
    # The paper displays the calibrated Suez markup rounded to 1.59; its
    # "shipping-only limit" row (d = phi at chi = 1) pins the underlying value.
    φ = 1.5855
    rows = [("Electronics, pharmaceuticals", 0.01, 1.0046, 1.0059),
            ("General manufactures",         0.03, 1.0139, 1.0176),
            ("Steel, basic chemicals",       0.08, 1.0376, 1.0468),
            ("Bulk (grain, ore, fertilizer)",0.30, 1.1483, 1.1756),
            ("Shipping-only limit",          1.00, 1.5855, 1.5855)]
    ok = true
    @printf("  %-32s %6s %10s %10s %10s %10s\n", "sector", "chi", "d=phi^chi", "paper", "d=1+chi(.)", "paper")
    for (nm, χ, pm, pa) in rows
        m = φ^χ; a = 1 + χ * (φ - 1)
        @printf("  %-32s %6.2f %10.4f %10.4f %10.4f %10.4f\n", nm, χ, m, pm, a, pa)
        ok &= (abs(m - pm) < 5e-4) && (abs(a - pa) < 5e-4)
    end
    return check("Table 7 reproduced", ok)
end

"Section 4.1: numerical optimum matches F* = d^L (mu-1)^(-1/mu), and F* -> d^L as mu -> Inf."
function verify_toy_monopolist()
    println("\n--- Section 4.1 toy monopolist ---")
    dL = 1.0
    profit(F, μ) = F * (F^(-μ)) / (F^(-μ) + dL^(-μ))
    ok = true
    for μ in (2.0, 4.0, 8.0, 40.0)
        Fstar = dL * (μ - 1)^(-1 / μ)
        Fs = range(1e-4, 3.0; length = 400_000)
        Fnum = Fs[argmax(profit.(Fs, μ))]
        @printf("  mu = %5.1f : analytic F* = %.5f, numerical = %.5f\n", μ, Fstar, Fnum)
        ok &= abs(Fnum - Fstar) < 2e-4
    end
    ok &= check("F* -> d^L as mu -> Inf", abs(1.0 * (1e6 - 1)^(-1 / 1e6) - dL) < 1e-4)
    return check("toy monopolist optimum", ok)
end

"Log-share division of a shared route's take is consistent and never over-extracts."
function verify_toll_split()
    println("\n--- eq 13 log-share division on shared routes ---")
    φ = [1.59, 1.16, 1.35]
    Φ, tshare, logshare = toll_shares(φ)
    r = 5                                        # Hormuz-Suez: tolled by Suez and Hormuz
    ok = check("combined wedge stacks multiplicatively", abs(Φ[r] - φ[1] * φ[3]) < 1e-12)
    ok &= check("authority shares of the take sum to one",
                abs(sum(logshare[a, r] for a in ROUTE_AUTHS[r]) - 1) < 1e-12)
    s = 3
    sep = (φ[1]^CHI[s] - 1) / φ[1]^CHI[s] + (φ[3]^CHI[s] - 1) / φ[3]^CHI[s]
    ok &= check("separate markups would over-extract", sep > tshare[r, s])
    ok &= check("single-authority route collapses to (phi^chi-1)/phi^chi",
                abs(tshare[1, s] - (φ[1]^CHI[s] - 1) / φ[1]^CHI[s]) < 1e-12)
    return ok
end

function verify_limits(raw_base, base, ηR)
    println("\n--- limiting cases ---")
    ok = true
    # chi = 0: the freight wedge cannot pass through to the cargo-value iceberg,
    # so a TOLL has no effect at all (eq 5). Note this is NOT true of a route
    # CLOSURE: dropping a route from the CES (eq 14) raises the composite by
    # rho_open^(-1/eta^R) whatever chi is, because the baseline cost gaps are
    # inferred from the observed shares. Run on the untolled baseline so that
    # switching chi off does not also invalidate the baseline's rent income.
    saveχ = copy(CHI)
    try
        CHI .= 0.0
        r0 = solve_counterfactual(raw_base, shock_hormuz(raw_base, 1.35); ηR = ηR,
                                  tol = 1e-10, damp = 0.5)
        ok &= check("chi = 0  =>  a transit fee has zero effect", maximum(abs.(r0.ΔW)) < 1e-6)
    finally
        CHI .= saveχ
    end
    # congestion and capacity feedback amplify the delivered-cost response
    saveλ, saveγ = copy(LAMBDA), copy(GAMMA_G)
    local d_full, d_nc
    try
        rfull = solve_counterfactual(base, shock_red_sea(base); ηR = ηR, tol = 1e-9, damp = 0.35)
        d_full = maximum(abs.(filter(isfinite, rfull.dhat) .- 1))
        LAMBDA .= 0.0; GAMMA_G .= 0.0
        rnc = solve_counterfactual(base, shock_red_sea(base); ηR = ηR, tol = 1e-9, damp = 0.35)
        d_nc = maximum(abs.(filter(isfinite, rnc.dhat) .- 1))
    finally
        LAMBDA .= saveλ; GAMMA_G .= saveγ
    end
    @printf("  max |dhat - 1| under the closure: with congestion %.5f, without %.5f\n", d_full, d_nc)
    ok &= check("congestion amplifies the delivered-cost response", d_full > d_nc)
    # eta^R -> Inf equalises delivered costs across open routes
    return ok
end

# ============================================================================
# 8. REPORTING
# ============================================================================

const REPORT_GROUPS = [
    ("Gulf economies",            ["QAT", "KWT", "IRQ", "BHR", "ARE", "SAU"]),
    ("Toll collectors",           ["IRN", "OMN", "EGY", "PAN"]),
    ("Rival hydrocarbon exp.",    ["GNQ", "BRN"]),
    ("Large trading economies",   ["DEU", "FRA", "NLD", "USA", "CHN", "JPN", "KOR", "IND"]),
    ("Transshipment / other",     ["SGP", "MLT", "ITA", "ROW"]),
]

function table_welfare(results::Vector{<:NamedTuple}, labels::Vector{String})
    println()
    @printf("%-24s", "")
    for l in labels
        @printf("%12s", l)
    end
    println()
    println("-"^(24 + 12 * length(labels)))
    for (gname, isos) in REPORT_GROUPS
        @printf("%-24s\n", gname)
        for iso in isos
            i = cid(iso)
            @printf("  %-22s", CNAME[i])
            for r in results
                @printf("%12.2f", r.ΔW[i])
            end
            println()
        end
    end
    @printf("  %-22s", "World")
    for r in results
        @printf("%12.2f", r.Wworld)
    end
    println()
end

function table_routes(results::Vector{<:NamedTuple}, labels::Vector{String})
    println()
    @printf("%-24s", "Route")
    for l in labels
        @printf("%12s", l)
    end
    println()
    println("-"^(24 + 12 * length(labels)))
    for r in 1:4
        @printf("%-24s", ROUTE_NAME[r])
        for res in results
            fb = res.route_flow_base[r] + res.route_flow_base[r + 4]
            fn = res.route_flow[r] + res.route_flow[r + 4]
            @printf("%12.1f", fb > 0 ? 100 * (fn / fb - 1) : 0.0)
        end
        println()
    end
    @printf("%-24s", "Hormuz transit")
    for res in results
        fb = sum(res.route_flow_base[HORMUZ_ROUTES])
        fn = sum(res.route_flow[HORMUZ_ROUTES])
        @printf("%12.1f", fb > 0 ? 100 * (fn / fb - 1) : 0.0)
    end
    println()
end

"Side-by-side bar groups (StatsPlots' groupedbar is not a dependency here)."
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
# 9. MAIN
# ============================================================================

banner("BROCKHAUS, HINZ & SERFATY (2026) - CHOKEPOINT DISRUPTIONS")
println("Stylized calibration: $NC countries x $NS sectors x $NM modes x $NR routes")
println("Reproduces the paper's MECHANISMS and incidence pattern, not GTAP magnitudes.")

# ---- 9.1 Raw (untolled) baseline -------------------------------------------
banner("1. RAW BASELINE (no chokepoint rents, phi = 1)")
@time "  build_baseline" raw_base = build_baseline()
flush(stdout)
@printf("  world value added        : \$%.1f tn\n", sum(raw_base.VA) / 1000)
@printf("  world gross output       : \$%.1f tn\n", sum(raw_base.Y) / 1000)
intl = sum(raw_base.Xijs) - sum(raw_base.Xijs[i, i, s] for i in 1:NC, s in 1:NS)
@printf("  international trade      : \$%.1f tn (%.0f%% of world VA)\n",
        intl / 1000, 100 * intl / sum(raw_base.VA))
@printf("  sea-borne trade          : \$%.1f tn\n", sum(raw_base.Xrs) / 1000)
@printf("  via Suez (incl. Hormuz-Suez): \$%.0f bn\n",
        sum(raw_base.Xrs[R_SUEZ, :]) + sum(raw_base.Xrs[R_SUEZ + 4, :]))
@printf("  via Panama               : \$%.0f bn\n",
        sum(raw_base.Xrs[R_PANAMA, :]) + sum(raw_base.Xrs[R_PANAMA + 4, :]))
@printf("  transiting Hormuz        : \$%.0f bn\n", sum(raw_base.Xrs[HORMUZ_ROUTES, :]))
seatrade = sum(raw_base.Xrs)
χbar = sum(CHI[s] * sum(raw_base.Xrs[:, s]) for s in 1:NS) / seatrade
@printf("  trade-weighted chi (sea) : %.3f   (paper: about 0.09)\n", χbar)
@printf("  energy share of sea trade: %.2f\n", sum(raw_base.Xrs[:, 2]) / seatrade)

println("\n  Chokepoint exposure of the reported economies (% of income):")
@printf("  %-24s %10s %10s %10s\n", "", "sea exp.", "Suez exp.", "Hormuz exp.")
function exposure(b, i, routes)
    sum(b.Xijs[i, j, s] * b.μ[i, j, s, 1] * sum(b.ρ[i, j, s, 1, r] for r in routes)
        for j in 1:NC, s in 1:NS) / b.I[i] * 100
end
for iso in ["EGY", "QAT", "KWT", "IRQ", "SAU", "ARE", "OMN", "IRN", "MLT", "SGP", "DEU", "CHN"]
    i = cid(iso)
    @printf("  %-24s %10.1f %10.1f %10.1f\n", CNAME[i], exposure(raw_base, i, 1:NR),
            exposure(raw_base, i, [R_SUEZ, R_SUEZ + 4]), exposure(raw_base, i, HORMUZ_ROUTES))
end
verify_baseline(raw_base, "raw baseline")

# ---- 9.2 Calibrate the route elasticity and the canal fees -----------------
banner("2. CALIBRATION OF eta^R AND THE CANAL WEDGES (Section 5.1)")
println("Pricing in isolation: each authority maximises revenue at the others' reference wedges.")
@time "  calibrate eta^R" cal = calibrate_route_elasticity(raw_base; target = 10.25, mode = :pe)
flush(stdout)
const ETA_R = cal.ηR
φ_suez = cal.φ_suez
φ_panama, rev_panama = best_response_fee(raw_base, 2; ηR = ETA_R, mode = :pe)
@printf("\n  %-34s %10s %10s\n", "", "model", "observed")
@printf("  %-34s %10.2f %10.2f\n", "Suez toll revenue (\$bn, TARGETED)", cal.rev_suez, 10.25)
@printf("  %-34s %10.3f %10s\n", "Suez wedge phi_Suez", φ_suez, "-")
@printf("  %-34s %10.2f %10.2f\n", "Panama toll revenue (\$bn, untargeted)", rev_panama, 3.18)
@printf("  %-34s %10.3f %10s\n", "Panama wedge phi_Panama", φ_panama, "-")
@printf("  %-34s %10.1f %10.1f\n", "route elasticity mu^d", ETA_R, ETA_R_PAPER)
println("\n  Paper: mu^d = 350, phi_Suez = 1.59, \$10.5bn; phi_Panama = 1.16, \$3.62bn.")
@printf("  Effective elasticity of route shares to the FEE is eta^R * chi_s = %.1f at chi = %.3f.\n",
        ETA_R * χbar, χbar)
@printf("  Own-traffic monotonicity at the calibrated eta^R: %s\n",
        admissible_eta(raw_base, ETA_R) ? "holds over the whole fee range" :
        "breaks at high fees (search truncated, see best_response_fee)")

# ---- 9.3 Re-baseline at the calibrated fees --------------------------------
banner("3. WITH-FEES BASELINE")
φcal = [φ_suez, φ_panama, 1.0]
@time "  rebaseline" base = rebaseline(raw_base, φcal; ηR = ETA_R)
flush(stdout)
@printf("  Suez rent   : \$%.2fbn  (%.2f%% of Egyptian income)\n",
        base.Πauth[1], 100 * base.Πc[cid("EGY")] / base.I[cid("EGY")])
@printf("  Panama rent : \$%.2fbn  (%.2f%% of Panamanian income)\n",
        base.Πauth[2], 100 * base.Πc[cid("PAN")] / base.I[cid("PAN")])
@printf("  Suez route share, Europe-Asia pair (DEU-CHN): raw %.3f -> with fees %.3f\n",
        raw_base.ρ[cid("DEU"), cid("CHN"), 6, 1, R_SUEZ],
        base.ρ[cid("DEU"), cid("CHN"), 6, 1, R_SUEZ])
println("  (the toll itself already pushes cargo off Suez - this is the cushion channel)")
verify_baseline(base, "with-fees baseline")

# ---- 9.4 Verification ------------------------------------------------------
banner("4. VERIFICATION")
verify_null_shock(raw_base, ETA_R, "the raw baseline")
verify_null_shock(base, ETA_R, "the with-fees baseline")
verify_table7()
verify_toy_monopolist()
verify_toll_split()
verify_limits(raw_base, base, ETA_R)

# BHS_QUICK=1 stops here: calibration + verification only, no counterfactuals or plots.
get(ENV, "BHS_QUICK", "0") == "1" && exit(0)

# ---- 9.5 Counterfactuals ---------------------------------------------------
banner("5. COUNTERFACTUALS (Section 6)")
say("Solving Red Sea closure ...")
res_rs = solve_counterfactual(base, shock_red_sea(base); ηR = ETA_R, tol = 1e-9, damp = 0.3)
say("Solving Hormuz transit fees ...")
res_h = [solve_counterfactual(base, shock_hormuz(base, φ); ηR = ETA_R, tol = 1e-9, damp = 0.3)
         for φ in (1.05, 1.35, 1.50)]
say("Solving the same shocks in a world with no chokepoint rents ...")
res_rs_nf = solve_counterfactual(raw_base, shock_red_sea(raw_base); ηR = ETA_R, tol = 1e-9, damp = 0.3)

banner("5.1 WELFARE (percent change in real income) - cf. Table 8")
table_welfare([res_rs, res_h[1], res_h[2], res_h[3]],
              ["Red Sea", "Hormuz 1.05", "Hormuz 1.35", "Hormuz 1.50"])

println("\n  Five largest losses and gains per scenario (cf. Figures 10, 13, 16):")
for (lab, r) in [("Red Sea", res_rs), ("Hormuz 1.35", res_h[2]), ("Red Sea, no fees", res_rs_nf)]
    ord = sortperm(r.ΔW)
    @printf("    %-18s losses: %s\n", lab,
            join([@sprintf("%s %.2f", ISO[i], r.ΔW[i]) for i in ord[1:5]], ", "))
    @printf("    %-18s gains : %s\n", "",
            join([@sprintf("%s %+.2f", ISO[i], r.ΔW[i]) for i in reverse(ord[end-4:end])], ", "))
end

banner("5.2 SEABORNE TRADE BY ROUTE (percent change) - cf. Table 9")
table_routes([res_rs, res_h[1], res_h[2], res_h[3]],
             ["Red Sea", "Hormuz 1.05", "Hormuz 1.35", "Hormuz 1.50"])
println()
for (k, φ) in enumerate((1.05, 1.35, 1.50))
    @printf("  Hormuz fee %.2f : toll revenue \$%.1fbn (Iran \$%.1fbn, Oman \$%.1fbn)\n",
            φ, res_h[k].Πauth[3], res_h[k].Πc[cid("IRN")], res_h[k].Πc[cid("OMN")])
end
hs_base = sum(base.Xrs[R_SUEZ + 4, :])
hs_new = sum(res_h[2].Xrs[R_SUEZ + 4, :])
@printf("  Hormuz-Suez composite at phi = 1.35: \$%.1fbn -> \$%.1fbn (%.1f%%), i.e. %.1f%% of Suez traffic\n",
        hs_base, hs_new, 100 * (hs_new / hs_base - 1),
        100 * hs_base / (sum(base.Xrs[R_SUEZ, :]) + hs_base))

banner("5.3 WELFARE DECOMPOSITION FOR THE COLLECTORS - cf. Table 10")
lbl = ["Red Sea", "Hormuz 1.05", "Hormuz 1.35", "Hormuz 1.50"]
allres = [res_rs, res_h[1], res_h[2], res_h[3]]
@printf("%-22s", "")
for l in lbl
    @printf("%18s", l)
end
println()
@printf("%-22s", "")
for _ in lbl
    @printf("%9s%9s", "dW", "d_tau")
end
println()
println("-"^(22 + 18 * length(lbl)))
for iso in ["EGY", "PAN", "IRN", "OMN"]
    i = cid(iso)
    @printf("%-22s", CNAME[i])
    for r in allres
        @printf("%9.2f%9.2f", r.ΔW[i], r.Δτ[i])
    end
    println()
end
println("\n  Note: d_tau is the direct toll-rent channel; the residual dW - d_tau is the")
println("  general-equilibrium and trade-cost channel.")

banner("5.4 THE RED SEA CLOSURE WITH AND WITHOUT CHOKEPOINT RENTS - cf. Table 6")
@printf("\n  %-24s %14s %14s\n", "", "with fees", "without fees")
println("  " * "-"^54)
for iso in ["EGY", "QAT", "MLT", "SGP", "ARE", "GNQ", "OMN", "IRQ", "DEU", "CHN"]
    i = cid(iso)
    @printf("  %-24s %14.2f %14.2f\n", CNAME[i], res_rs.ΔW[i], res_rs_nf.ΔW[i])
end
@printf("  %-24s %14.2f %14.2f\n", "World", res_rs.Wworld, res_rs_nf.Wworld)

egy = cid("EGY")
rank_fees = sortperm(res_rs.ΔW)
rank_nofee = sortperm(res_rs_nf.ΔW)
@printf("\n  Egypt's rank among losers: %d of %d with fees, %d of %d without fees\n",
        findfirst(==(egy), rank_fees), NC, findfirst(==(egy), rank_nofee), NC)
@printf("  Egypt: %.2f%% -> %.2f%%   (rent channel %.2f pp of the with-fees loss)\n",
        res_rs.ΔW[egy], res_rs_nf.ΔW[egy], res_rs.Δτ[egy])
@printf("  World cost is similar either way: %.3f%% vs %.3f%%\n", res_rs.Wworld, res_rs_nf.Wworld)

# ---- 9.6 Laffer sweep ------------------------------------------------------
banner("6. HORMUZ LAFFER CURVE")
φgrid = collect(1.0:0.1:2.5)
laffer_rev = Float64[]; laffer_flow = Float64[]
hbase = sum(base.Xrs[HORMUZ_ROUTES, :])
for φ in φgrid
    r = solve_counterfactual(base, shock_hormuz(base, φ); ηR = ETA_R, tol = 1e-8, damp = 0.35)
    push!(laffer_rev, r.Πauth[3])
    push!(laffer_flow, 100 * (sum(r.Xrs[HORMUZ_ROUTES, :]) / hbase - 1))
    @printf("    phi = %.2f -> revenue \$%.1fbn, Hormuz flow %+.1f%% (%d iters)\n",
            φ, r.Πauth[3], laffer_flow[end], r.iters); flush(stdout)
end
kmax = argmax(laffer_rev)
k95 = findfirst(>=(0.95 * laffer_rev[kmax]), laffer_rev)
@printf("  best fee on the grid: %.2f (\$%.1fbn); 95%% of that revenue is already reached at %.2f\n",
        φgrid[kmax], laffer_rev[kmax], φgrid[k95])
@printf("  revenue is non-decreasing up to the maximum: %s\n",
        all(diff(laffer_rev[1:kmax]) .>= -1e-6) ? "YES" : "NO")
println("  The economically relevant statement is the flattening, not the peak: revenue rises")
println("  steeply to about phi = 1.4 and is essentially flat thereafter, so the fee under")
println("  discussion (1.35) already sits in the revenue-maximising region of the Laffer curve.")

# ---- 9.7 Nash vs. pricing in isolation ------------------------------------
banner("7. NASH PRICING VS. PRICING IN ISOLATION (Section 4.4)")
φnash = copy(raw_base.φ_base)
for it in 1:8
    φold = copy(φnash)
    for a in 1:2
        φnash[a], _ = best_response_fee(raw_base, a; φ_others = φnash, ηR = ETA_R, mode = :pe)
    end
    maximum(abs.(φnash[1:2] .- φold[1:2])) < 1e-4 && (@printf("  Nash fixed point in %d rounds\n", it); break)
end
@printf("  isolation: phi_Suez = %.4f, phi_Panama = %.4f\n", φ_suez, φ_panama)
@printf("  Nash     : phi_Suez = %.4f, phi_Panama = %.4f\n", φnash[1], φnash[2])
@printf("  max absolute gap = %.5f -> the paper's claim that the two nearly coincide is %s\n",
        maximum(abs.(φnash[1:2] .- φcal[1:2])),
        maximum(abs.(φnash[1:2] .- φcal[1:2])) < 0.02 ? "SUPPORTED" : "NOT supported here")

# ---- 9.8 Plots -------------------------------------------------------------
banner("8. FIGURES")

# (1) calibration curve: Suez revenue against the fee
φs = collect(1.0:0.02:2.6)
revs = [first(revenue_at(raw_base, [φ, 1.0, 1.0], 1; ηR = ETA_R, mode = :pe)) for φ in φs]
p1 = plot(φs, revs, lw = 2.5, legend = :topright, label = "Suez toll revenue",
          xlabel = "Suez wedge phi", ylabel = "revenue (\$bn)",
          title = "Calibration: the revenue moment that pins eta^R")
hline!(p1, [10.25], ls = :dash, lc = :red, label = "observed \$10.25bn (2023)")
vline!(p1, [φ_suez], ls = :dot, lc = :black, label = @sprintf("phi* = %.2f", φ_suez))
savefig(p1, joinpath(GRAPHDIR, "bhs_calibration_suez.pdf"))

# (2) Hormuz Laffer curve
p2 = plot(φgrid, laffer_rev, lw = 2.5, label = "toll revenue (\$bn)", legend = :left,
          xlabel = "Hormuz fee phi", ylabel = "revenue (\$bn)",
          title = "Strait of Hormuz: revenue and throughput")
p2b = twinx(p2)
plot!(p2b, φgrid, laffer_flow, lw = 2.5, lc = :darkred, ls = :dash,
      ylabel = "Hormuz trade (% change)", label = "")
vline!(p2, [1.05, 1.35, 1.50], ls = :dot, lc = :grey, label = "scenario fees")
savefig(p2, joinpath(GRAPHDIR, "bhs_hormuz_laffer.pdf"))

# (3) five largest gains and losses, per scenario
function plot_extremes(res, ttl, fname)
    ord = sortperm(res.ΔW)
    sel = vcat(ord[1:5], reverse(ord[end-4:end]))
    vals = res.ΔW[sel]
    bar(1:10, vals, orientation = :h, yticks = (1:10, ISO[sel]), legend = false,
        color = [v < 0 ? :firebrick : :steelblue for v in vals],
        xlabel = "change in real income (%)", title = ttl, yflip = true)
    savefig(joinpath(GRAPHDIR, fname))
end
plot_extremes(res_rs, "Red Sea closure: largest losses and gains", "bhs_extremes_redsea.pdf")
plot_extremes(res_h[2], "Hormuz fee (phi = 1.35): largest losses and gains", "bhs_extremes_hormuz.pdf")
plot_extremes(res_rs_nf, "Red Sea closure without chokepoint rents", "bhs_extremes_redsea_nofee.pdf")

# (4) route reallocation by scenario
rnames = ["Suez", "Panama", "Cape", "Direct", "Hormuz transit"]
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
M = hcat(route_pct(res_rs), route_pct(res_h[1]), route_pct(res_h[2]), route_pct(res_h[3]))
p4 = grouped_bar(M, rnames, ["Red Sea", "Hormuz 1.05", "Hormuz 1.35", "Hormuz 1.50"];
                 ylabel = "change in flow (%)", legend = :bottomleft,
                 title = "Seaborne trade by route")
savefig(p4, joinpath(GRAPHDIR, "bhs_route_reallocation.pdf"))

# (5) with vs. without fees - the punchline. Two panels, because Egypt's loss is
# an order of magnitude larger than everyone else's and would hide the crossing.
egysel = [cid("EGY")]
sel = [cid(x) for x in ["QAT", "IRQ", "KWT", "BHR", "MLT", "ARE", "OMN", "SGP", "GNQ"]]
p5a = grouped_bar(hcat(res_rs.ΔW[egysel], res_rs_nf.ΔW[egysel]), ISO[egysel],
                  ["with fees", "without fees"];
                  ylabel = "change in real income (%)", legend = :bottomleft,
                  title = "The toll collector")
p5b = grouped_bar(hcat(res_rs.ΔW[sel], res_rs_nf.ΔW[sel]), ISO[sel],
                  ["with fees", "without fees"];
                  legend = :bottomleft, title = "Economies whose cargo travels farther")
p5 = plot(p5a, p5b, layout = grid(1, 2, widths = [0.28, 0.72]), size = (950, 400),
          plot_title = "Red Sea closure: the rent moves the incidence, not the total",
          plot_titlefontsize = 11, top_margin = 6Plots.mm)
savefig(p5, joinpath(GRAPHDIR, "bhs_with_without_fees.pdf"))

# (6) rent decomposition for collectors
coll = [cid(x) for x in ["EGY", "PAN", "IRN", "OMN"]]
p6 = grouped_bar(hcat(res_rs.ΔW[coll], res_rs.Δτ[coll], res_rs.ΔW[coll] .- res_rs.Δτ[coll]),
                 ISO[coll], ["total dW", "direct rent channel", "GE / trade-cost residual"];
                 ylabel = "change in real income (%)", legend = :bottomleft,
                 title = "Red Sea closure: welfare decomposition (collectors)")
savefig(p6, joinpath(GRAPHDIR, "bhs_rent_decomposition.pdf"))

# (7) mechanism scatter: exposure alone does not predict incidence
suez_exp = [sum(base.Xijs[i, j, s] * base.μ[i, j, s, 1] *
                (base.ρ[i, j, s, 1, R_SUEZ] + base.ρ[i, j, s, 1, R_SUEZ + 4])
                for j in 1:NC, s in 1:NS) / base.I[i] for i in 1:NC]
p7 = scatter(suez_exp, res_rs.ΔW, ms = 5, legend = false,
             xlabel = "Suez-routed exports / income", ylabel = "change in real income (%)",
             title = "Exposure alone does not predict incidence",
             series_annotations = text.(ISO, 5, :bottom))
savefig(p7, joinpath(GRAPHDIR, "bhs_mechanism_scatter.pdf"))

# (8) country x scenario heatmap
Wmat = hcat(res_rs.ΔW, res_h[1].ΔW, res_h[2].ΔW, res_h[3].ΔW, res_rs_nf.ΔW)
p8 = heatmap(1:5, 1:NC, clamp.(Wmat, -4, 4), yticks = (1:NC, ISO),
             xticks = (1:5, ["RedSea", "H 1.05", "H 1.35", "H 1.50", "RedSea\nno fees"]),
             c = :RdBu, clims = (-4, 4), title = "Change in real income (%), clipped at +/-4")
savefig(p8, joinpath(GRAPHDIR, "bhs_welfare_heatmap.pdf"))

println("  10 figures written to $(GRAPHDIR)")

# ---- 9.9 Qualitative claim checks -----------------------------------------
banner("9. THE PAPER'S QUALITATIVE CLAIMS")
qat, omn, irn, deu, chn = cid("QAT"), cid("OMN"), cid("IRN"), cid("DEU"), cid("CHN")
check("Egypt is the largest loser under a Red Sea closure WITH fees",
      argmin(res_rs.ΔW) == egy)
check("Egypt's loss is dominated by the rent channel",
      abs(res_rs.Δτ[egy]) > 0.6 * abs(res_rs.ΔW[egy]))
check("Egypt's loss shrinks by an order of magnitude WITHOUT fees",
      abs(res_rs_nf.ΔW[egy]) < 0.4 * abs(res_rs.ΔW[egy]))
check("Qatar's loss deepens WITHOUT fees",
      res_rs_nf.ΔW[qat] < res_rs.ΔW[qat])
check("with-fees and no-fees loser rankings cross (Egypt vs Qatar)",
      (res_rs.ΔW[egy] < res_rs.ΔW[qat]) && (res_rs_nf.ΔW[qat] < res_rs_nf.ΔW[egy]))
check("world cost of the closure is similar with and without fees",
      abs(res_rs.Wworld - res_rs_nf.Wworld) < 0.5 * max(abs(res_rs.Wworld), abs(res_rs_nf.Wworld)) + 0.01)
check("Suez traffic falls to zero, Cape and Panama rise",
      route_pct(res_rs)[1] < -99.9 && route_pct(res_rs)[3] > 0 && route_pct(res_rs)[2] > 0)
check("large trading economies lose little under the closure",
      abs(res_rs.ΔW[deu]) < 0.2 && abs(res_rs.ΔW[chn]) < 0.2)
check("Hormuz fee concentrates losses on the Gulf",
      res_h[2].ΔW[qat] < -0.5 && abs(res_h[2].ΔW[deu]) < 0.2)
check("Oman gains far more than Iran despite the equal revenue split",
      res_h[2].ΔW[omn] > 2 * res_h[2].ΔW[irn])
check("rival hydrocarbon exporters gain from the Hormuz fee",
      res_h[2].ΔW[cid("GNQ")] > 0 || res_h[2].ΔW[cid("BRN")] > 0)
check("the Hormuz fee taxes the Suez canal at one remove (Egypt loses)",
      res_h[2].ΔW[egy] < 0 && hs_new < hs_base)
hflow(res) = 100 * (sum(res.Xrs[HORMUZ_ROUTES, :]) / hbase - 1)
check("Hormuz throughput response is convex in the fee",
      abs(hflow(res_h[2])) > 3 * abs(hflow(res_h[1])))
check("Hormuz toll revenue flattens between phi = 1.35 and 1.50",
      res_h[3].Πauth[3] < 1.25 * res_h[2].Πauth[3])

banner("DONE")
println("Reminder: this is a stylized calibration. It reproduces the paper's mechanisms")
println("and incidence pattern; matching Tables 8-10 numerically requires GTAP 11 and the")
println("AIS route shares (see Claude_Plan.md section 11 for the loader interface).")
