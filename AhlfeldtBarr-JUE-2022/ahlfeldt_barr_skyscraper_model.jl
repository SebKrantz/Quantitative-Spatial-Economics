# Stylized Monocentric City with an ENDOGENOUS BUILDING-HEIGHT MARGIN
# Based on Ahlfeldt, Gabriel M. & Jason Barr (2022), "The economics of skyscrapers: A synthesis",
# Journal of Urban Economics 129, 103419.  https://doi.org/10.1016/j.jue.2021.103419
#
# Julia port of the MIT-licensed AB2022 Stata toolkit (see LICENSE-AB2022-toolkit).
# Source: stata_source/AB2022.ado (shipped ado solver) and stata_source/_1_PROGS.do (walkthrough),
# with AB2022-codebook.pdf (SS A.1, "13 equations"; Table 1) as the authoritative documentation.
#
# ---------------------------------------------------------------------------------------------
# THE MODEL
# ---------------------------------------------------------------------------------------------
# A stylized, perfectly-open, ONE-DIMENSIONAL monocentric city: Alonso-Mills-Muth with an
# endogenous VERTICAL margin bolted on. Space is a line x in [-50, 50] km, D = |x|. Two competing
# land uses, commercial (C) and residential (R), plus an agricultural outside option at exogenous
# rent r_a. Developers choose building height to maximize profit per unit land given the
# floor-space price; land goes to whichever use bids the highest land rent; a single city-wide
# wage y and total employment N adjust until commercial labour demand equals residential supply.
#
# There is NO commuting-cost term. Distance enters ONLY through the exponential decay
# exp(-tau^U * D) of the productivity/amenity fundamentals, and one wage applies everywhere.
# tau^C/tau^R are a reduced-form decay of agglomeration/amenity value with distance, NOT a
# commuting-cost primitive. (The Codebook lists a reservation utility Ubar among the primitives,
# but it appears nowhere in the toolkit code; openness is implemented by making N a target object
# rather than fixing utility, so Ubar is absorbed into a^R(x).)
#
# ---------------------------------------------------------------------------------------------
# THE 13 EQUATIONS (Codebook SS A.1). U in {C, R}; Greek in comments, ASCII in code.
# ---------------------------------------------------------------------------------------------
#   Shifters (Table 1):  A~^U(x) = abar^U * atilde^U * N^beta^U * exp(-tau^U * |x|)
#   (1) a^C(x)   = A~^C(x)^(1/(1-alpha^C)) * y^(alpha^C/(alpha^C-1))      commercial floor-price shifter
#   (2) a^R(x)   = A~^R(x)^(1/(1-alpha^R)) * y^(1/(1-alpha^R))            residential floor-price shifter
#   (3) S~^C(x)  = min( (a^C/(c^C(1+theta^C)))^(1/(theta^C-omega^C)), Sbar^C )   constrained height
#   (4) S~^R(x)  = min( (a^R/(c^R(1+theta^R)))^(1/(theta^R-omega^R)), Sbar^R )
#   (5) r^C(x)   = a^C/(1+omega^C) * (S~^C)^(1+omega^C) - c^C * (S~^C)^(1+theta^C)   commercial land rent
#   (6) r^R(x)   = a^R/(1+omega^R) * (S~^R)^(1+omega^R) - c^R * (S~^R)^(1+theta^R)   residential land rent
#   (7) S^C = S~^C, S^R = 0   if r^C >= r^R and r^C >= r_a               highest-and-best use
#   (8) S^R = S~^R, S^C = 0   if r^R >  r^C and r^R >= r_a
#   (9) pbar^C(x) = a^C(x)/(1+omega^C) * S^C(x)^omega^C                  horizontal (per-floor) rent
#  (10) pbar^R(x) = a^R(x)/(1+omega^R) * S^R(x)^omega^R
#  (11) L(x)    = alpha^C/(1-alpha^C) * pbar^C(x)/y * S^C(x)             workplace employment (MRS)
#  (12) n(x)    = S^R(x) / fbar^R(x),  fbar^R = (1-alpha^R)/pbar^R * y   residence employment (Marshallian)
#  (13) int_{-x0}^{x0} L(x)dx = int_{-x1}^{-x0} n(x)dx + int_{x0}^{x1} n(x)dx = N   labour clearing
#
# Profit-maximizing (unconstrained) height, just above Codebook eq. (3):
#       S*^U(x) = (a^U/(c^U(1+theta^U)))^(1/(theta^U-omega^U))
#
# ---------------------------------------------------------------------------------------------
# !! NOTATION COLLISION WARNING -- READ BEFORE COPYING ANY CALIBRATED VALUE OUT OF THIS FILE !!
# ---------------------------------------------------------------------------------------------
# These letters mean DIFFERENT things here than elsewhere in Quantitative-Spatial-Economics/.
# A value transplanted across models will look plausible and be silently wrong:
#
#   theta  HERE: construction-cost elasticity of height,  0.5-0.55
#          ELSEWHERE: Frechet / trade elasticity, 6-10  (AA2022, Redding-JIE-2016, BHS2026)
#   beta   HERE: AGGLOMERATION elasticity,  +0.03
#          ELSEWHERE: AMENITY externality,  -0.3  (AA2022) -- opposite sign AND opposite economics
#   alpha  HERE: non-floor-space input / expenditure share, 0.66-0.85
#          ELSEWHERE: productivity externality, 0.1  (AA2022)
#   omega  HERE: rent elasticity of height, omega^U = omegatilde^U/(1-alpha^U). Not used elsewhere.
#   lambda NOT USED HERE. Elsewhere: traffic-congestion elasticity 0.07-0.09 (AA2022) -- and, in
#          the Ahlfeldt/ARSW2015 toolkit, an AGGLOMERATION elasticity at nearly the same magnitude.
#   sigma, epsilon  NOT USED HERE (no CES demand system, no migration elasticity).
#
# ---------------------------------------------------------------------------------------------
# !! UPSTREAM INCONSISTENCY: THE BID-RENT omega CONVENTION (eqs. 9-10) !!
# ---------------------------------------------------------------------------------------------
# The two Stata sources implement DIFFERENT formulas for the horizontal bid rent:
#
#   stata_source/AB2022.ado:97-98      p_bar_x_C = a_x_C * 1/(1-omega_C) * S_x_C^omega_C
#   stata_source/_1_PROGS.do:132-133   p_bar_x_C = a_x_C * 1/(1+omega_C) * S_x_C^omega_C
#
# AB2022-codebook.pdf eqs. 9-10 read a(x)/(1+omega^U) -- the Codebook agrees with _1_PROGS.do,
# against the .ado.
#
# (1+omega) is the THEORETICALLY CORRECT form, and this is settled without reference to any data.
# The developer maximizes pbar(S)*S - c*S^(1+theta). Under pbar = a/(1+omega)*S^omega, revenue is
# a/(1+omega)*S^(1+omega), so the FOC a*S^omega = c(1+theta)*S^theta gives
# S* = (a/(c(1+theta)))^(1/(theta-omega)) -- EXACTLY the S* both files use. Under (1-omega) the FOC
# would instead give S* = (a(1+omega)/((1-omega)c(1+theta)))^(1/(theta-omega)), which neither file
# implements. So the .ado is internally inconsistent: its pbar does not generate its own S*, and
# pbar*S is not the first term of its own land rent (eq. 5).
#
# Both are implemented here, selected by `bid_rent_convention`; the DEFAULT is :codebook.
#   :codebook -> 1/(1+omega)   Codebook eqs. 9-10, _1_PROGS.do:132-133, and the FOC above.
#   :ado      -> 1/(1-omega)   AB2022.ado:97-98.
# Level effect measured by this port (see validate_inverted / omega_convention_effect below):
# switching :codebook -> :ado raises pbar^C by exactly (1+omega^C)/(1-omega^C) and pbar^R by
# (1+omega^R)/(1-omega^R).
#
# DECISIVE EVIDENCE: the shipped reference_data/INVERTED.csv was generated with (1-omega), i.e. the
# .ado convention -- but with the _1_PROGS/Codebook LAND-RENT bracket. It is a hybrid that neither
# current source file reproduces on its own. reference_data/BASE.csv cannot settle the question
# because it holds NO solved values at all (see validate_base). Details in TRANSLATION_PLAN.md SS4.
#
# Land rent (eqs. 5-6) is NOT offered as a convention: AB2022.ado:69-70 writes c^C*(1+omega^C)
# inside the first bracket where the Codebook and _1_PROGS.do:104-105 write c^C*(1+theta^C). The
# .ado form is inconsistent with its own S* (line 89) and is ~400% off the shipped data. Treated
# here as a plain bug; only the Codebook form is implemented.
#
# ---------------------------------------------------------------------------------------------
# RELATION TO THE OTHER MODELS IN THIS REPOSITORY
# ---------------------------------------------------------------------------------------------
# This is the ONLY model in Quantitative-Spatial-Economics/ with an endogenous building-height
# margin. In AllenArkolakis-RES-2022, FuchsFoongWong-MMN-2026, Santamaria-2022, QRE-HoRaUE-2025,
# Redding-JIE-2016, QSE-ARE-2017 and BrockhausHinzSerfaty-2026, floor/land supply is fixed, an
# exogenous elasticity, or absent, and every location is a point with only a horizontal
# trade/migration problem. The developer's height problem and the construction-cost elasticity
# theta^U are entirely absent from those gravity/trade/traffic models.
#
# Conversely, transport cost is EXOGENOUS here (indeed absent -- see "no commuting cost" above),
# which is precisely the margin those models endogenize.
#
#   ABSJ2027 (Ahlfeldt, Baum-Snow & Jedwab, ReStud forthcoming) is AB2022's direct successor:
#     same vertical margin and 1-D geometry, but imperfectly open (migration closure) and with an
#     explicit urban growth boundary alongside height caps, calibrated to 12,873 cities.
#   ARSW2015 (Berlin blocks) and MRRH2018 (German counties) DROP the vertical margin: floor space
#     there is MEASURED data, not a developer's choice, and they add bilateral commuting instead.
# See Ahlfeldt/CLAUDE.md for the full ladder.
#
# ---------------------------------------------------------------------------------------------
# WHAT IS IMPLEMENTED
# ---------------------------------------------------------------------------------------------
#   Algorithm 1 SOLVER   -> solver!    one deterministic recursive pass given (y, N)
#   Algorithm 2 WAGE     -> wage!      inner loop on the wage
#   Algorithm 3 FINDEQ   -> findeq!    outer loop on employment, nests wage!
#   Algorithm 4 CONV     -> conv!      amenity update toward observed heights
#   Algorithm 5 EMP      -> emp!       amenity scaling toward a target population
#   Algorithm 6 INVERT   -> invert!    drives conv! then emp! (Chicago skyline)
#   Counterfactuals: height limit (_2_ANALYSIS.do:38), commercial subcenter (_4_INVERTEDCOUNTER.do)
#   Validation against reference_data/BASE.csv and reference_data/INVERTED.csv
#
# Run with:  julia ahlfeldt_barr_skyscraper_model.jl

using Statistics, Printf, Random, Plots, DataFrames, CSV

const REFDIR = joinpath(@__DIR__, "reference_data")
const GRAPHDIR = joinpath(@__DIR__, "graphs")

# ============================================================
# Parameters
# ============================================================

"""
    Params

Structural parameters. Defaults are the canonical Ahlfeldt & Barr (2022) values hard-coded in
`AB2022.ado:378-392`. Note that the toolkit uses THREE different parameterizations and they are
not interchangeable:

| context                           | c_C, c_R | r_a | source                |
|-----------------------------------|----------|-----|-----------------------|
| `AB2022.ado` defaults (here)      | 1.4      | 150 | `AB2022.ado:387-389`  |
| `_2_ANALYSIS.do` baseline figures | 1.3      | 50  | `_2_ANALYSIS.do:23`   |
| `_1_PROGS.do` CONV/EMP inversion  | 1.4      | 30  | `_1_PROGS.do:276,288` |

`_1_PROGS.do:17-33` also sets beta_C=0.04, omega=0.1, c=1, a_bar_C=1.5 while building BASE.dta;
those are placeholders overwritten the instant FINDEQ is called, and are not reproduced here.
"""
Base.@kwdef struct Params
    alpha_C::Float64 = 0.85     # alpha^C  commercial share of the non-floor-space input
    alpha_R::Float64 = 0.66     # alpha^R  residential expenditure share on the non-housing good
    beta_C::Float64  = 0.03     # beta^C   commercial AGGLOMERATION elasticity (see collision warning)
    beta_R::Float64  = 0.00     # beta^R   residential agglomeration elasticity (off by default)
    tau_C::Float64   = 0.01     # tau^C    commercial productivity decay, per km
    tau_R::Float64   = 0.005    # tau^R    residential amenity decay, per km
    omega_C::Float64 = 0.03     # omega^C  commercial rent elasticity of height
    omega_R::Float64 = 0.07     # omega^R  residential rent elasticity of height
    theta_C::Float64 = 0.5      # theta^C  commercial CONSTRUCTION-COST elasticity of height
    theta_R::Float64 = 0.55     # theta^R  residential construction-cost elasticity of height
    c_C::Float64     = 1.4      # c^C      commercial baseline construction cost
    c_R::Float64     = 1.4      # c^R      residential baseline construction cost
    a_bar_C::Float64 = 2.0      # abar^C   fundamental commercial productivity
    a_bar_R::Float64 = 1.0      # abar^R   fundamental residential amenity
    r_a::Float64     = 150.0    # r_a      agricultural land rent (a VARIABLE in Stata, not a scalar)
    S_bar_C::Float64 = 999.0    # Sbar^C   commercial height limit (999 ~ unconstrained)
    S_bar_R::Float64 = 999.0    # Sbar^R   residential height limit
    bid_rent_convention::Symbol = :codebook   # :codebook -> 1/(1+omega); :ado -> 1/(1-omega)
    constrained_land_rent::Bool = true        # true = Codebook eq. 5 (CORRECT); false reproduces the toolkit's bug (see solver! step 3)
end

# Parameterizations used by the walkthrough, for reproducibility of its figures.
params_analysis() = Params(c_C=1.3, c_R=1.3, r_a=50.0)        # _2_ANALYSIS.do:23
params_conv()     = Params(c_C=1.4, c_R=1.4, r_a=30.0)        # _1_PROGS.do:276, 288

"Bid-rent denominator in Codebook eqs. 9-10. See the omega-convention block in the header."
bid_denom(omega, convention) = convention === :codebook ? (1 + omega) : (1 - omega)

# ============================================================
# City construction
# ============================================================

"""
    City

The 1-D grid and every endogenous object, mirroring the Stata dataset one column per field.
Missing is represented by `NaN`, matching Stata's `.` for storage (but NOT for comparison --
see `assign_land_use!`).
"""
mutable struct City
    # exogenous grid
    x::Vector{Float64}; D::Vector{Float64}; r_a::Vector{Float64}
    a_rand_C::Vector{Float64}; a_rand_R::Vector{Float64}     # atilde^U, location-specific component
    # target objects
    y::Float64; L::Float64
    # endogenous fields
    A_tilde_x_C::Vector{Float64}; A_tilde_x_R::Vector{Float64}
    a_x_C::Vector{Float64}; a_x_R::Vector{Float64}
    r_x_C::Vector{Float64}; r_x_R::Vector{Float64}
    U::Vector{Float64}
    S_star_x_C::Vector{Float64}; S_star_x_R::Vector{Float64}
    S_x_C::Vector{Float64}; S_x_R::Vector{Float64}; S_x::Vector{Float64}
    p_bar_x_C::Vector{Float64}; p_bar_x_R::Vector{Float64}
    L_x_C::Vector{Float64}; f_bar_x_R::Vector{Float64}; n_x::Vector{Float64}
    URBAN::Vector{Float64}; COM::Vector{Float64}
    # scalars set by solver!
    L_hat_demand::Float64; L_hat_supply::Float64
    sL::Float64; sy::Float64            # scalar sL/sy (AB2022.ado:113-116) -- see note in findeq!
    x0::Float64; x1::Float64
end

"""
    build_city(; large=false, y0=2.5, L0=1e6, r_a0=100.0)

Reproduce `AB2022.ado:332-375` / `_1_PROGS.do:35-78`. `large=false` gives 10,001 points,
x = (n-5001)/100 for n = 1..10001, i.e. [-50, 50] km in 0.01 km steps. Matching this spacing
exactly is what makes a column-by-column comparison against the reference CSVs meaningful rather
than a measurement of interpolation error.
"""
function build_city(; large::Bool=false, y0=2.5, L0=1e6, r_a0=100.0)
    n = large ? 20001 : 10001
    mid = large ? 10001 : 5001
    x = [(i - mid) / 100 for i in 1:n]
    nanv() = fill(NaN, n)
    City(x, abs.(x), fill(r_a0, n), ones(n), ones(n), y0, L0,
         nanv(), nanv(), nanv(), nanv(), nanv(), nanv(), nanv(),
         nanv(), nanv(), nanv(), nanv(), nanv(), nanv(), nanv(),
         nanv(), nanv(), nanv(), nanv(), nanv(),
         NaN, NaN, NaN, NaN, NaN, NaN)
end

# --- Stata missing-value semantics -------------------------------------------------------------
# Stata's sum/egen ignore missing; Stata's min()/max() return the non-missing argument and are
# missing only if all arguments are. Julia's min/max propagate NaN, so these are explicit.
nansum(v) = sum(x for x in v if !isnan(x); init = 0.0)
nanmax2(a, b) = isnan(a) ? b : (isnan(b) ? a : max(a, b))

# ============================================================
# Height / bid-rent solver (Algorithm 1: SOLVER)
# ============================================================

"Smallest x >= 0 whose land use satisfies `pred`; NaN if there is none. (Stata `sum ... if`, r(min).)"
function min_x_where(c::City, pred)
    m = Inf
    @inbounds for i in eachindex(c.x)
        c.x[i] >= 0 && pred(c.U[i]) && (m = min(m, c.x[i]))
    end
    return isfinite(m) ? m : NaN
end

"""
    assign_land_use!(c)

Codebook eqs. 7-8 / `AB2022.ado:73-76` / `_1_PROGS.do:108-111`: highest-and-best use, assigned by
three sequential overwrites in the order agricultural (3), residential (2), commercial (1).

STATA SEMANTICS NOTE. Stata's missing `.` is LARGER than any number in comparisons, whereas Julia's
`NaN` compares false against everything. The two diverge only when exactly one side is missing.
That case cannot arise here: `a^U >= 0` always, and `a^U = 0` yields `r^U = 0` rather than missing,
so `r_x_C`/`r_x_R` are never NaN after the land-rent step. The assertion below enforces that rather
than leaving it to luck.
"""
function assign_land_use!(c::City)
    @assert !any(isnan, c.r_x_C) && !any(isnan, c.r_x_R) "land rent is NaN: Stata/Julia missing-comparison semantics diverge here"
    @inbounds for i in eachindex(c.U)
        rC, rR, ra = c.r_x_C[i], c.r_x_R[i], c.r_a[i]
        u = NaN
        ra > rC && ra > rR && (u = 3.0)          # agricultural rent highest
        rR > rC && rR > ra && (u = 2.0)          # residential rent highest
        rC > rR && rC > ra && (u = 1.0)          # commercial rent highest
        c.U[i] = u
    end
end

"""
    solver!(c, p)

Algorithm 1 (`AB2022.ado:54-124`, `_1_PROGS.do:90-158`). Given guesses of the target objects
(y, N), compute every endogenous object in one deterministic recursive pass. No iteration.
"""
function solver!(c::City, p::Params)
    # -- clear pre-existing values (AB2022.ado:58, _1_PROGS.do:93) ------------------------------
    # UPSTREAM QUIRK, DELIBERATELY REPLICATED: the Stata clear-list has 14 variables and OMITS
    # both S_star_x_C and S_star_x_R, so profit-maximizing heights retain STALE values from
    # earlier iterations at locations whose land use has since changed. In the shipped
    # INVERTED.dta this leaves S_star_x_C non-missing at 1,387 rows where U != 1 (and S_star_x_R
    # at 2,135 rows where U != 2). The reference data depends on this, so S_star is not cleared.
    for v in (c.A_tilde_x_C, c.A_tilde_x_R, c.a_x_C, c.a_x_R, c.r_x_C, c.r_x_R, c.U,
              c.S_x_C, c.S_x_R, c.p_bar_x_C, c.p_bar_x_R, c.L_x_C, c.f_bar_x_R, c.n_x)
        fill!(v, NaN)
    end

    y, L = c.y, c.L

    # -- 1-2. shifters (Table 1) and floor-space price shifters, Codebook eqs. 1-2 --------------
    @. c.A_tilde_x_C = p.a_bar_C * c.a_rand_C * L^p.beta_C * exp(-p.tau_C * c.D)
    @. c.A_tilde_x_R = p.a_bar_R * c.a_rand_R * L^p.beta_R * exp(-p.tau_R * c.D)
    @. c.a_x_C = c.A_tilde_x_C^(1 / (1 - p.alpha_C)) * y^(p.alpha_C / (p.alpha_C - 1))
    @. c.a_x_R = c.A_tilde_x_R^(1 / (1 - p.alpha_R)) * y^(1 / (1 - p.alpha_R))

    # -- 3. land rent, Codebook eqs. 5-6 --------------------------------------------------------
    # r^U = a^U/(1+omega) * (S~^U)^(1+omega) - c^U * (S~^U)^(1+theta),  S*^U = (a^U/(c(1+theta)))^(1/(theta-omega))
    #
    # NOTE the bracket is c*(1+theta), per the Codebook and _1_PROGS.do:104-105. AB2022.ado:69-70
    # writes c*(1+omega) there, which contradicts its own S* and the shipped data. See header.
    #
    # THIRD UPSTREAM INCONSISTENCY, FOUND BY THIS PORT (not in the toolkit's own notes):
    # both Stata files evaluate land rent at the UNCONSTRAINED S*, because they inline the closed
    # form (a/(c(1+theta)))^((1+omega)/(theta-omega)) and S_bar never appears in it. Codebook eq. 5
    # specifies the CONSTRAINED S~ = min(S_bar, S*). The two coincide whenever the height limit
    # does not bind -- which is the case for the baseline (S_bar = 999) and for every column of
    # the reference data, so Tiers A-E are unaffected. They diverge sharply once a cap binds: at
    # the centre with S_bar_C = 20 the toolkit reports land rent 409.26 where Codebook eq. 5 gives
    # 204.90, a 99.7% overstatement. Since land rent drives the land-use allocation (eqs. 7-8),
    # the CBD and urban boundaries of the height-limit counterfactual inherit that error.
    # Default is `true` = Codebook eq. 5, i.e. the CORRECT land rent. Set `constrained_land_rent=false`
    # to reproduce the toolkit's behaviour exactly (needed only for bit-comparison against Stata output;
    # the reference data never binds a cap, so Tiers A-E pass identically either way).
    Sstar_C = @. (c.a_x_C / (p.c_C * (1 + p.theta_C)))^(1 / (p.theta_C - p.omega_C))
    Sstar_R = @. (c.a_x_R / (p.c_R * (1 + p.theta_R)))^(1 / (p.theta_R - p.omega_R))
    Stil_C = p.constrained_land_rent ? min.(p.S_bar_C, Sstar_C) : Sstar_C
    Stil_R = p.constrained_land_rent ? min.(p.S_bar_R, Sstar_R) : Sstar_R
    @. c.r_x_C = c.a_x_C / (1 + p.omega_C) * Stil_C^(1 + p.omega_C) - p.c_C * Stil_C^(1 + p.theta_C)
    @. c.r_x_R = c.a_x_R / (1 + p.omega_R) * Stil_R^(1 + p.omega_R) - p.c_R * Stil_R^(1 + p.theta_R)

    # -- 4. land-use allocation, Codebook eqs. 7-8 ---------------------------------------------
    assign_land_use!(c)

    # -- 5. boundaries -------------------------------------------------------------------------
    # x1 = inner margin of the agricultural zone; x0 = inner margin of any non-commercial zone.
    # UPSTREAM DIVERGENCE: _1_PROGS.do:114-115 takes min(x | U==3, x>=0); AB2022.ado:79-80 takes
    # max(x | U<3, x>=0). These differ by exactly one grid step (0.01 km). The _1_PROGS form is
    # used here. Both are reported diagnostics only -- neither enters any equation.
    c.x1 = min_x_where(c, u -> u == 3)
    c.x0 = min_x_where(c, u -> !(u == 1))

    # -- 6. profit-maximizing and realized height, Codebook eqs. 3-4 ---------------------------
    @inbounds for i in eachindex(c.U)
        if c.U[i] == 1
            c.S_star_x_C[i] = Sstar_C[i]           # NOT cleared between passes -- see above
            c.S_x_C[i] = min(p.S_bar_C, Sstar_C[i])
        elseif c.U[i] == 2
            c.S_star_x_R[i] = Sstar_R[i]
            c.S_x_R[i] = min(p.S_bar_R, Sstar_R[i])
        end
    end

    # -- 7. horizontal bid rent, Codebook eqs. 9-10 (THE omega CONVENTION) ---------------------
    dC = bid_denom(p.omega_C, p.bid_rent_convention)
    dR = bid_denom(p.omega_R, p.bid_rent_convention)
    @. c.p_bar_x_C = c.a_x_C / dC * c.S_x_C^p.omega_C
    @. c.p_bar_x_R = c.a_x_R / dR * c.S_x_R^p.omega_R

    # -- 8. labour demand, Codebook eq. 11 (MRS) + eq. 13 aggregation --------------------------
    @inbounds for i in eachindex(c.U)
        c.U[i] == 1 && (c.L_x_C[i] = p.alpha_C / (1 - p.alpha_C) * c.p_bar_x_C[i] / y * c.S_x_C[i])
    end
    c.L_hat_demand = nansum(c.L_x_C)

    # -- 9. labour supply, Codebook eq. 12 (Marshallian) + eq. 13 aggregation ------------------
    @inbounds for i in eachindex(c.U)
        c.U[i] == 2 && (c.f_bar_x_R[i] = (1 - p.alpha_R) / c.p_bar_x_R[i] * y)
    end
    @. c.n_x = c.S_x_R / c.f_bar_x_R
    c.L_hat_supply = nansum(c.n_x)

    # -- 10. final statistics ------------------------------------------------------------------
    c.sL, c.sy = L, y
    @. c.S_x = nanmax2(c.S_x_C, c.S_x_R)
    @. c.URBAN = Float64(c.U < 3)
    @. c.COM = Float64(c.U == 1)
    return c
end

# ============================================================
# General equilibrium (Algorithms 2-3: WAGE, FINDEQ)
# ============================================================

"""
    wage!(c, p; tol=0.01, maxiter=1000, verbose=false)

Algorithm 2 (`AB2022.ado:131-160`). Inner loop: adjust the wage until aggregate labour demand
equals aggregate supply to within `tol` (1% relative, the upstream value).

UPSTREAM DEFECT FIXED: neither Stata loop has a maxiter guard, so a parameterization that never
reaches tolerance without tripping the `L_hat_demand + L_hat_supply == 0` check hangs Stata
forever. `maxiter` + `@warn` added here.

The `while !(obj <= tol)` form is deliberate: it keeps iterating on NaN, matching Stata, where a
missing objective (0/0) compares as larger than any tolerance.
"""
function wage!(c::City, p::Params; tol=0.01, maxiter::Int=1000, verbose::Bool=false)
    obj = abs(c.L_hat_demand / c.L_hat_supply - 1)     # first objective uses +0 (AB2022.ado:133)
    iter = 0
    while !(obj <= tol)
        iter += 1
        if iter > maxiter
            @warn "wage!: no convergence in $maxiter iterations (objective $(round(obj, digits=5)) > $tol)"
            break
        end
        y_factor = c.L_hat_supply == 0 ? 1.2 :
                   c.L_hat_demand == 0 ? 0.8 :
                   (c.L_hat_demand / c.L_hat_supply)^0.01
        c.y = 0.5 * c.y + 0.5 * c.y * y_factor
        solver!(c, p)
        # subsequent objectives use +0.0001 (AB2022.ado:155) -- faithful to the source
        obj = abs((c.L_hat_demand + 0.0001) / (c.L_hat_supply + 0.0001) - 1)
        verbose && @printf("    wage iter %3d  y = %.6f  obj = %.5f\n", iter, c.y, obj)
    end
    return c
end

"""
    findeq!(c, p; tol=0.01, maxiter=1000, verbose=false)

Algorithm 3 (`AB2022.ado:414-446`, `_1_PROGS.do:202-256`). Outer loop on total employment,
nesting `wage!`.

OFF-BY-ONE, FAITHFULLY PRESERVED: the objective is recomputed BEFORE `L` is updated, and the
`L` update is the last statement in the loop body. So on exit the stored `L` is ONE UPDATE AHEAD
of the endogenous columns it is paired with -- those were produced by the previous `L`. The
scalar `c.sL` holds that previous value (Stata's `scalar sL`, set inside SOLVER). This is not
cosmetic: reproducing the shipped INVERTED.dta requires recovering
`L_prev = 2*(L_stored - 0.25*(L_hat_demand + L_hat_supply))`, without which `A_tilde_x_C` is off
by beta_C*ln(L_stored/L_prev) = 1.4e-4 relative instead of the 8.5e-8 float32 floor.
"""
function findeq!(c::City, p::Params; tol=0.01, maxiter::Int=1000, verbose::Bool=false)
    fill!(c.r_a, p.r_a)              # Stata: `qui replace r_a = `12'` -- a VARIABLE, not a scalar
    solver!(c, p)
    obj_ext = abs(c.L / (0.5 * (c.L_hat_demand + c.L_hat_supply)) - 1)
    iter = 0
    while !(obj_ext <= tol)
        iter += 1
        if iter > maxiter
            @warn "findeq!: no convergence in $maxiter iterations (objective $(round(obj_ext, digits=5)) > $tol)"
            break
        end
        solver!(c, p)
        wage!(c, p; tol=tol, maxiter=maxiter)
        if c.L_hat_demand + c.L_hat_supply == 0
            error("City does not reach critical size. Increase productivity to make the city more attractive.")
        end
        obj_ext = abs(c.L / (0.5 * (c.L_hat_demand + c.L_hat_supply)) - 1)
        c.L = 0.5 * c.L + 0.25 * (c.L_hat_demand + c.L_hat_supply)   # LAST -- see off-by-one note
        verbose && @printf("  outer iter %3d  L = %.2f  y = %.6f  obj = %.5f\n", iter, c.L, c.y, obj_ext)
    end
    return (iters=iter, obj_ext=obj_ext)
end

# ============================================================
# Inversion against an observed skyline (Algorithms 4-6)
# ============================================================

"Stata `round(x, 0.1)`: nearest multiple, halves away from zero."
round_bin(x, w=0.1) = w * round(x / w, RoundNearestTiesAway)

"""
    prepare_chicago(c, path)

Reproduce `_3_INVERSION.do:22-48`: merge the observed Chicago skyline onto the grid by `X`, keep
for each location only the TALLER of the commercial/residential observation, then bin at 0.1 km
and spread the within-bin maximum. Returns `(HEIGHT_C, HEIGHT_R, CONVBIN)`.
"""
function prepare_chicago(c::City, path::AbstractString)
    sky = CSV.read(path, DataFrame)
    n = length(c.x)
    HC, HR = fill(NaN, n), fill(NaN, n)
    idx = Dict(round(Int, xi * 100) => i for (i, xi) in enumerate(c.x))
    for r in eachrow(sky)                                   # merge 1:1 on X
        i = get(idx, round(Int, r.X * 100), 0)
        i == 0 && continue
        HC[i] = ismissing(r.HEIGHT_C) ? NaN : r.HEIGHT_C
        HR[i] = ismissing(r.HEIGHT_R) ? NaN : r.HEIGHT_R
    end
    for i in 1:n                                            # :39-40 keep the taller observation
        if !isnan(HC[i]) && !isnan(HR[i])
            HC[i] >= HR[i] ? (HR[i] = NaN) : (HC[i] = NaN)
        end
    end
    CONVBIN = round_bin.(c.x)                               # :44  (X == x for every master row)
    binmax(v) = (d = Dict{Float64,Float64}();
                 for i in 1:n; isnan(v[i]) || (d[CONVBIN[i]] = nanmax2(get(d, CONVBIN[i], NaN), v[i])); end; d)
    bC = binmax(HC)                                         # :45-46 spread commercial bin maxima
    for i in 1:n
        isnan(HR[i]) && (HC[i] = get(bC, CONVBIN[i], NaN))
    end
    bR = binmax(HR)                                         # :47-48 spread residential bin maxima
    for i in 1:n
        isnan(HC[i]) && (HR[i] = get(bR, CONVBIN[i], NaN))
    end
    return HC, HR, CONVBIN
end

"""
    conv!(c, p, lambda, HC, HR)

Algorithm 4 (`_1_PROGS.do:268-278`). Nudge the location-specific amenity components toward what
would rationalize the observed heights, using `HEIGHT^U / S*^U` as the adjustment factor, then
re-solve. Where either the observed height or the (possibly stale) `S*` is missing, the amenity
is set to a theory-consistent zero.
"""
function conv!(c::City, p::Params, lambda::Float64, HC::Vector{Float64}, HR::Vector{Float64})
    for (a_rand, H, Sstar) in ((c.a_rand_R, HR, c.S_star_x_R), (c.a_rand_C, HC, c.S_star_x_C))
        @inbounds for i in eachindex(a_rand)
            v = (1 - lambda) * a_rand[i] + lambda * (H[i] / Sstar[i]) * a_rand[i]
            a_rand[i] = isnan(v) ? 0.0 : v
        end
    end
    findeq!(c, p)
end

"""
    emp!(c, p, target)

Algorithm 5 (`_1_PROGS.do:284-290`). Scale residential amenities toward a target population.

UPSTREAM QUIRK, REPLICATED: the Stata body sits inside `foreach name in R C` but never uses
`name`, so it executes TWICE per call. That changes the scaling, so it is preserved. `c.sL` (not
`c.L`) is the employment the last solver pass actually used -- see the off-by-one note in findeq!.
"""
function emp!(c::City, p::Params, target::Float64)
    for _ in 1:2
        c.a_rand_R .*= (target / c.sL)^0.01
        findeq!(c, p)
    end
end

"R-squared of `reg HEIGHT_R S_x_R` -- squared correlation over complete cases (Stata drops missing)."
function r2_height(H::Vector{Float64}, S::Vector{Float64})
    k = [i for i in eachindex(H) if !isnan(H[i]) && !isnan(S[i])]
    length(k) < 3 && return 0.0
    h, s = H[k], S[k]
    (std(h) == 0 || std(s) == 0) && return 0.0
    return cor(h, s)^2
end

"""
    invert!(c, p, target, lambda, HC, HR; maxiter=200, verbose=true)

Algorithm 6 (`_1_PROGS.do:296-316`). Iterate `conv!` until the model height gradient correlates
with the data at R^2 >= 0.999, then iterate `emp!` until the population gap is within 1,000.

`_3_INVERSION.do:91` calls `INVERT 1000000 0.05`, and INVERT passes its SECOND argument through as
CONV's convergence parameter, so lambda = 0.05. (Ahlfeldt/AB2022-toolkit/CLAUDE.md states 0.5;
that is a typo in those upstream notes.)

maxiter guards added -- neither upstream `while` has one.
"""
function invert!(c::City, p::Params, target::Float64, lambda::Float64,
                 HC::Vector{Float64}, HR::Vector{Float64}; maxiter::Int=3000, verbose::Bool=true)
    popgap = abs(c.sL - target)          # initialized BEFORE the CONV loop, as upstream
    O, nconv = 0.0, 0
    while O < 0.999
        nconv += 1
        nconv > maxiter && (@warn "invert!: height correlation stalled at $O after $maxiter CONV steps"; break)
        conv!(c, p, lambda, HC, HR)
        O = r2_height(HR, c.S_x_R)
        verbose && (nconv <= 3 || nconv % 10 == 0) && @printf("    CONV %3d  R2 = %.6f\n", nconv, O)
    end
    nemp = 0
    while popgap > 1000
        nemp += 1
        nemp > maxiter && (@warn "invert!: population gap stalled at $popgap after $maxiter EMP steps"; break)
        emp!(c, p, target)
        popgap = abs(c.sL - target)
        verbose && (nemp <= 3 || nemp % 50 == 0) && @printf("    EMP  %4d  population gap = %.1f\n", nemp, popgap)
    end
    return (r2=O, popgap=popgap, n_conv=nconv, n_emp=nemp)
end

# ============================================================
# Validation against the Stata reference data
# ============================================================

"Read a CSV column as Float64 with NaN for missing."
function nanvec(df::DataFrame, col::Symbol)
    v = df[!, col]
    return [ismissing(x) ? NaN : Float64(x) for x in v]
end

"Max absolute and max relative deviation, comparing only where the reference is finite."
function deviation(got::Vector{Float64}, ref::Vector{Float64}; mask=nothing)
    k = [i for i in eachindex(ref) if !isnan(ref[i]) && (mask === nothing || mask[i])]
    isempty(k) && return (n=0, maxabs=0.0, maxrel=0.0, nmiss=0)
    nmiss = count(i -> isnan(got[i]), k)
    k = [i for i in k if !isnan(got[i])]
    isempty(k) && return (n=0, maxabs=NaN, maxrel=NaN, nmiss=nmiss)
    maxabs = maximum(abs(got[i] - ref[i]) for i in k)
    maxrel = maximum(ref[i] == 0 ? 0.0 : abs((got[i] - ref[i]) / ref[i]) for i in k)
    return (n=length(k), maxabs=maxabs, maxrel=maxrel, nmiss=nmiss)
end

# The reference CSVs carry ~8 significant decimal digits, printed from Stata's 4-byte floats, so
# each stored value already carries up to ~1.2e-07 relative quantization (float32 half-ulp 5.96e-08
# plus the decimal printing). A comparison chain touches several such values, so this is the
# irreducible noise floor for any column whose sensitivity to the coarsely-stored y is zero.
const STORAGE_FLOOR = 4 * eps(Float32) / 2      # ~2.4e-07

const BASE_COLS = [:x, :D, :y, :L, :r_a, :a_x_C, :a_x_R, :a_rand_C, :a_rand_R, :A_tilde_x_C,
                   :A_tilde_x_R, :r_x_C, :r_x_R, :U, :S_star_x_C, :S_star_x_R, :S_x_C, :S_x_R,
                   :S_x, :p_bar_x_C, :p_bar_x_R, :L_x_C, :f_bar_x_R, :n_x, :URBAN, :COM,
                   :SHADE, :SHADEU]

"""
    validate_base()

Tier A. Reproduce `reference_data/BASE.csv`, all 28 columns.

FINDING: BASE.csv is NOT a solved equilibrium. `_1_PROGS.do:81` saves BASE.dta BETWEEN generating
the placeholder variables (:48-78) and defining SOLVER (:90), so 21 of its 28 columns are
ENTIRELY missing (0 non-null out of 10,001 each) and only the grid and starting values carry
information. Consequently BASE.csv cannot discriminate the omega convention and cannot supply
land-use boundaries -- both of those move to INVERTED.csv. Note also r_a = 100 here, the `gen`
placeholder from _1_PROGS.do:53, which is none of the three solver parameterizations.
"""
function validate_base()
    df = CSV.read(joinpath(REFDIR, "BASE.csv"), DataFrame)
    c = build_city(large=false, y0=2.5, L0=1e6, r_a0=100.0)   # unsolved, as saved
    got = Dict(:x=>c.x, :D=>c.D, :y=>fill(c.y, length(c.x)), :L=>fill(c.L, length(c.x)),
               :r_a=>c.r_a, :a_rand_C=>c.a_rand_C, :a_rand_R=>c.a_rand_R,
               :a_x_C=>c.a_x_C, :a_x_R=>c.a_x_R, :A_tilde_x_C=>c.A_tilde_x_C,
               :A_tilde_x_R=>c.A_tilde_x_R, :r_x_C=>c.r_x_C, :r_x_R=>c.r_x_R, :U=>c.U,
               :S_star_x_C=>c.S_star_x_C, :S_star_x_R=>c.S_star_x_R, :S_x_C=>c.S_x_C,
               :S_x_R=>c.S_x_R, :S_x=>c.S_x, :p_bar_x_C=>c.p_bar_x_C, :p_bar_x_R=>c.p_bar_x_R,
               :L_x_C=>c.L_x_C, :f_bar_x_R=>c.f_bar_x_R, :n_x=>c.n_x, :URBAN=>c.URBAN,
               :COM=>c.COM, :SHADE=>fill(NaN, length(c.x)), :SHADEU=>fill(NaN, length(c.x)))
    println("\n", "="^78)
    println("TIER A -- reference_data/BASE.csv  (", nrow(df), " rows x ", ncol(df), " cols)")
    println("="^78)
    @printf("%-14s %8s %14s %14s   %s\n", "column", "n_ref", "max abs dev", "max rel dev", "status")
    ok = true
    for col in BASE_COLS
        ref = nanvec(df, col)
        nref = count(!isnan, ref)
        if nref == 0
            allmiss = all(isnan, got[col])
            allmiss || (ok = false)
            @printf("%-14s %8d %14s %14s   %s\n", col, 0, "--", "--",
                    allmiss ? "all-missing both sides OK" : "MISMATCH: port not missing")
        else
            d = deviation(got[col], ref)
            d.maxrel > 1e-7 && (ok = false)
            @printf("%-14s %8d %14.3e %14.3e   %s\n", col, nref, d.maxabs, d.maxrel,
                    d.maxrel <= 1e-7 ? "OK" : "FAIL")
        end
    end
    println("\nPopulated columns: 7 of 28.  The other 21 are all-missing in the reference because")
    println("BASE.dta is saved BEFORE the solver ever runs (_1_PROGS.do:81).")
    println(ok ? ">>>> TIER A PASSED <<<<" : ">>>> TIER A FAILED <<<<")
    return ok
end

const SOLVED_COLS = [:A_tilde_x_C, :A_tilde_x_R, :a_x_C, :a_x_R, :r_x_C, :r_x_R, :U,
                     :S_star_x_C, :S_star_x_R, :S_x_C, :S_x_R, :S_x, :p_bar_x_C, :p_bar_x_R,
                     :L_x_C, :f_bar_x_R, :n_x, :URBAN, :COM]

"""
    validate_inverted(; convention=:ado, quiet=false)

Tier B -- the real test. `reference_data/INVERTED.csv` is the only SOLVED gold standard in the
folder: all 24 solved columns populated over 4,001 rows (|x| <= 20 after `_3_INVERSION.do:119`).

Take its own exogenous state as given -- a_rand_C, a_rand_R, y, r_a = 30, and the CONV/EMP
parameters from `_1_PROGS.do:276` -- run ONE solver pass, and compare every solved column. This
tests all 13 equations simultaneously while side-stepping the path-dependence of a 1%-tolerance
fixed point (see Tier C).

Two replication details are mandatory, both non-obvious:
  (a) the stored L is ONE UPDATE AHEAD of the columns it generated (see findeq!), so the solve
      uses L_prev = 2*(L_stored - 0.25*(L_hat_demand + L_hat_supply));
  (b) S_star is never cleared upstream, so it is compared ONLY within its own land-use zone.

PRECISION FLOOR: Stata stores gen-created variables as 4-byte floats. Round-tripping INVERTED.csv
through Float32 perturbs it by up to 5.7e-08 relative, so no column can do better than ~6e-8.
"""
function validate_inverted(; convention::Symbol=:ado, quiet::Bool=false)
    path = joinpath(REFDIR, "INVERTED.csv")
    df = CSV.read(path, DataFrame)
    n = nrow(df)
    p = Params(c_C=1.4, c_R=1.4, r_a=30.0, bid_rent_convention=convention)

    L_hat_demand = nansum(nanvec(df, :L_x_C))
    L_hat_supply = nansum(nanvec(df, :n_x))
    L_prev = 2 * (df.L[1] - 0.25 * (L_hat_demand + L_hat_supply))     # (a) recover L_prev

    # rebuild the 4,001-point window exactly as saved, with the reference's own exogenous state
    function solve_at(yval)
        cc = City(nanvec(df, :x), nanvec(df, :D), nanvec(df, :r_a),
                  nanvec(df, :a_rand_C), nanvec(df, :a_rand_R), yval, L_prev,
                  (fill(NaN, n) for _ in 1:19)..., NaN, NaN, NaN, NaN, NaN, NaN)
        solver!(cc, p)
        return cc
    end
    c = solve_at(df.y[1])

    # ---- the precision band the reference data itself permits ----------------------------------
    # The CSV stores y with only SIX significant digits (e.g. "2.10701"), and a_x_C scales as
    # y^(alpha_C/(alpha_C-1)) = y^-5.667, which is then raised to 1/(theta_C-omega_C) = 2.128 for
    # height and compounded again for labour demand. Half an ulp of that printed y is therefore
    # amplified by up to ~21x. Re-solving at the two ends of y's own rounding interval bounds how
    # precisely ANY correct implementation could reproduce these columns.
    ndec = (f = split(strip(split(readlines(path)[2], ",")[3]), "."); length(f) > 1 ? length(f[2]) : 0)
    half_ulp = 0.5 * 10.0^(-ndec)
    clo, chi = solve_at(df.y[1] - half_ulp), solve_at(df.y[1] + half_ulp)

    got = Dict(:A_tilde_x_C=>c.A_tilde_x_C, :A_tilde_x_R=>c.A_tilde_x_R, :a_x_C=>c.a_x_C,
               :a_x_R=>c.a_x_R, :r_x_C=>c.r_x_C, :r_x_R=>c.r_x_R, :U=>c.U,
               :S_star_x_C=>c.S_star_x_C, :S_star_x_R=>c.S_star_x_R, :S_x_C=>c.S_x_C,
               :S_x_R=>c.S_x_R, :S_x=>c.S_x, :p_bar_x_C=>c.p_bar_x_C, :p_bar_x_R=>c.p_bar_x_R,
               :L_x_C=>c.L_x_C, :f_bar_x_R=>c.f_bar_x_R, :n_x=>c.n_x, :URBAN=>c.URBAN, :COM=>c.COM)
    refU = nanvec(df, :U)

    if !quiet
        println("\n", "="^94)
        println("TIER B -- reference_data/INVERTED.csv  ($n rows), bid_rent_convention = :$convention")
        println("="^94)
        @printf("  L_stored = %.4f   L_prev (used) = %.6f   ratio = %.8f\n", df.L[1], L_prev, df.L[1]/L_prev)
        @printf("  L_hat_demand = %.4f   L_hat_supply = %.4f\n", L_hat_demand, L_hat_supply)
        @printf("  reference y = %s (%d significant digits) => half-ulp %.1e, the precision band below\n\n",
                df.y[1], ndec + 1, half_ulp)
        @printf("%-14s %7s %13s %13s %13s %7s  %s\n",
                "column", "n_ref", "max abs dev", "max rel dev", "tolerance", "dev/tol", "status")
    end
    worst, worst_ratio, ok = 0.0, 0.0, true
    for col in SOLVED_COLS
        ref = nanvec(df, col)
        # (b) S_star is stale outside its own zone -- compare only where the zone is active
        mask = col === :S_star_x_C ? (refU .== 1) :
               col === :S_star_x_R ? (refU .== 2) : nothing
        d = deviation(got[col], ref; mask=mask)
        lo, hi = getfield(clo, col), getfield(chi, col)
        k = [i for i in eachindex(ref) if !isnan(ref[i]) && ref[i] != 0 && !isnan(lo[i]) &&
             !isnan(hi[i]) && (mask === nothing || mask[i])]
        band = isempty(k) ? 0.0 : maximum(abs((hi[i] - lo[i]) / 2 / ref[i]) for i in k)
        floor_col = max(band, STORAGE_FLOOR)     # y-band, or pure storage noise where band == 0
        ratio = d.maxrel / floor_col
        worst = max(worst, d.maxrel); worst_ratio = max(worst_ratio, ratio)
        pass = d.maxrel <= floor_col
        pass || (ok = false)
        quiet || @printf("%-14s %7d %13.3e %13.3e %13.3e %7.2f  %s\n",
                         col, d.n, d.maxabs, d.maxrel, floor_col, ratio, pass ? "OK" : "FAIL")
    end
    if !quiet
        for col in (:x, :D, :r_a, :a_rand_C, :a_rand_R)
            d = deviation(getfield(c, col), nanvec(df, col))
            @printf("%-14s %7d %13.3e %13.3e %13s %7s  %s\n",
                    col, d.n, d.maxabs, d.maxrel, "-", "-", "OK (exogenous)")
        end
        nU = count(i -> c.U[i] == refU[i], 1:n)
        println()
        @printf("  Land use U identical at %d / %d rows\n", nU, n)
        @printf("  Boundaries: x0 = %.2f km, x1 = %.2f km\n", c.x0, c.x1)
        @printf("  Worst relative deviation: %.3e; worst dev/tolerance ratio: %.2f\n", worst, worst_ratio)
        println("  PASS CRITERION: deviation <= what the reference's own stored precision permits --")
        println("  the y half-ulp band above, floored at %s (storage noise). A ratio well below 1" |> x -> replace(x, "%s" => @sprintf("%.1e", STORAGE_FLOOR)))
        println("  means the port agrees more closely than the reference data can resolve.")
        println(ok ? ">>>> TIER B PASSED <<<<" : ">>>> TIER B FAILED <<<<")
    end
    return (ok=ok, worst=worst, worst_ratio=worst_ratio, city=c, L_hat_demand=L_hat_demand,
            L_hat_supply=L_hat_supply, L_stored=df.L[1],
            nU_match=count(i -> c.U[i] == refU[i], 1:n), n=n)
end


"Tier C: the equilibrium conditions Stata itself tests, evaluated at the reproduced state."
function validate_equilibrium(v)
    println("\n", "="^78); println("TIER C -- general equilibrium at the reproduced state"); println("="^78)
    outer = abs(v.L_stored / (0.5 * (v.L_hat_demand + v.L_hat_supply)) - 1)
    inner = abs(v.L_hat_demand / v.L_hat_supply - 1)
    @printf("  Labour demand  L_hat = %.2f\n", v.L_hat_demand)
    @printf("  Labour supply  N_hat = %.2f\n", v.L_hat_supply)
    @printf("  Demand/supply gap          |L_hat/N_hat - 1|            = %.5f  (< 0.01? %s)\n",
            inner, inner < 0.01 ? "YES" : "NO")
    @printf("  Outer objective            |L/(0.5(L_hat+N_hat)) - 1|   = %.5f  (< 0.01? %s)\n",
            outer, outer < 0.01 ? "YES" : "NO")
    println("\n  NOTE: FINDEQ stops at a 1% relative tolerance with fixed damping, so its terminal")
    println("  (y, N) is not a precise fixed point -- it is path-dependent, and a re-solve from a")
    println("  different start legitimately lands elsewhere inside the 1% band. That is a property")
    println("  of the source, not of this port.")
    return inner < 0.01 && outer < 0.01
end

# ============================================================
# Counterfactuals
# ============================================================

"""
    omega_convention_effect(p)

Tier F. Solve the same city under both bid-rent conventions and report the level shift. The
pointwise ratio on floor-space rent is exactly (1+omega)/(1-omega) by construction; the interest
is in how that propagates through general equilibrium to wage, population and the boundaries.
"""
function omega_convention_effect(pbase::Params)
    out = Dict{Symbol,Any}()
    for conv in (:codebook, :ado)
        p = Params(; (f => getfield(pbase, f) for f in fieldnames(Params) if f !== :bid_rent_convention)...,
                   bid_rent_convention=conv)
        c = build_city()
        findeq!(c, p)
        out[conv] = c
    end
    a, b = out[:codebook], out[:ado]
    println("\n", "="^78); println("TIER F -- effect of the bid-rent omega convention"); println("="^78)
    @printf("%-26s %14s %14s %12s\n", "", ":codebook (1+w)", ":ado (1-w)", "% change")
    pct(u, v) = 100 * (v / u - 1)
    mx(v) = (w = [t for t in v if !isnan(t)]; isempty(w) ? NaN : maximum(w))
    @printf("%-26s %14.4f %14.4f %11.2f%%\n", "max floor rent, commercial", mx(a.p_bar_x_C), mx(b.p_bar_x_C), pct(mx(a.p_bar_x_C), mx(b.p_bar_x_C)))
    @printf("%-26s %14.4f %14.4f %11.2f%%\n", "max floor rent, residential", mx(a.p_bar_x_R), mx(b.p_bar_x_R), pct(mx(a.p_bar_x_R), mx(b.p_bar_x_R)))
    @printf("%-26s %14.4f %14.4f %11.2f%%\n", "wage y", a.sy, b.sy, pct(a.sy, b.sy))
    @printf("%-26s %14.0f %14.0f %11.2f%%\n", "employment N", a.sL, b.sL, pct(a.sL, b.sL))
    @printf("%-26s %14.2f %14.2f %11.2f%%\n", "CBD radius x0 (km)", a.x0, b.x0, pct(a.x0, b.x0))
    @printf("%-26s %14.2f %14.2f %11.2f%%\n", "urban radius x1 (km)", a.x1, b.x1, pct(a.x1, b.x1))
    @printf("\n  Theoretical pointwise ratio on p_bar: (1+w)/(1-w) = %.4f (C), %.4f (R)\n",
            (1 + pbase.omega_C) / (1 - pbase.omega_C), (1 + pbase.omega_R) / (1 - pbase.omega_R))
    return a, b
end

"""
    height_limit_counterfactual(p; S_bar=20.0)

`_2_ANALYSIS.do:38`. A height limit does not merely cap the skyline: constrained developers cannot
bid as much for land, so the CBD and urban footprint shrink and city-wide wage and population fall,
while floor-space price and land-rent LEVELS rise.
"""
function height_limit_counterfactual(pbase::Params; S_bar=20.0)
    base = build_city(); findeq!(base, pbase)
    p_hl = Params(; (f => getfield(pbase, f) for f in fieldnames(Params) if f ∉ (:S_bar_C, :S_bar_R))...,
                  S_bar_C=S_bar, S_bar_R=S_bar)
    hl = build_city(); findeq!(hl, p_hl)
    # same cap, but with the toolkit's (incorrect) unconstrained land rent -- see solver! step 3
    p_up = Params(; (f => getfield(p_hl, f) for f in fieldnames(Params) if f !== :constrained_land_rent)...,
                  constrained_land_rent=false)
    up = build_city(); findeq!(up, p_up)
    println("\n", "="^78); @printf("COUNTERFACTUAL -- height limit of %.0f floors on both uses\n", S_bar); println("="^78)
    @printf("%-26s %14s %14s %12s\n", "", "baseline", "height limit", "% change")
    pct(u, v) = 100 * (v / u - 1)
    mx(v) = (w = [t for t in v if !isnan(t)]; isempty(w) ? NaN : maximum(w))
    for (lab, ba, ha) in (("max commercial height", mx(base.S_x_C), mx(hl.S_x_C)),
                          ("max residential height", mx(base.S_x_R), mx(hl.S_x_R)),
                          ("max floor rent (C)", mx(base.p_bar_x_C), mx(hl.p_bar_x_C)),
                          ("max land rent (C)", mx(base.r_x_C), mx(hl.r_x_C)),
                          ("wage y", base.sy, hl.sy),
                          ("employment N", base.sL, hl.sL),
                          ("CBD radius x0 (km)", base.x0, hl.x0),
                          ("urban radius x1 (km)", base.x1, hl.x1))
        @printf("%-26s %14.4f %14.4f %11.2f%%\n", lab, ba, ha, pct(ba, ha))
    end
    println("\n  The table above uses Codebook eq. 5 (land rent at the CONSTRAINED height), the port's")
    println("  default. The toolkit evaluates land rent at the UNCONSTRAINED S* instead -- see solver! step 3:")
    @printf("    max land rent (C): %10.4f  (Codebook eq. 5) vs %10.4f  (toolkit)\n", mx(hl.r_x_C), mx(up.r_x_C))
    @printf("    employment N:      %10.0f  (Codebook eq. 5) vs %10.0f  (toolkit)\n", hl.sL, up.sL)
    @printf("    CBD radius x0:     %10.2f  (Codebook eq. 5) vs %10.2f  (toolkit)\n", hl.x0, up.x0)
    @printf("    urban radius x1:   %10.2f  (Codebook eq. 5) vs %10.2f  (toolkit)\n", hl.x1, up.x1)
    return base, hl
end

"""
    subcenter_counterfactual(inv, p, CONVBIN; seed=2022)

`_4_INVERTEDCOUNTER.do:27-33`. Starting from the Chicago-inverted fundamentals, inject a positive
commercial-amenity shock between 12 and 15 km from the centre, spread it to 0.1 km location bins,
and re-solve. A new commercial subcenter emerges with its own height gradient; the original CBD
LOSES height (spatial reallocation), but the city overall grows, because it is more productive.

NOT REPRODUCIBLE UPSTREAM: `_4_INVERTEDCOUNTER.do:27` calls `runiform()` with no seed, so the shock
differs on every Stata run. A seed is set here so this port at least reproduces itself; the levels
are therefore not comparable to any particular upstream run.
"""
function subcenter_counterfactual(inv::City, p::Params, CONVBIN::Vector{Float64}; seed::Int=2022)
    before = (sL=inv.sL, sy=inv.sy, x0=inv.x0, x1=inv.x1,
              S_C=[isnan(t) ? 0.0 : t for t in inv.S_x_C])
    rng = MersenneTwister(seed)
    cf = deepcopy(inv)
    @inbounds for i in eachindex(cf.x)
        12 < cf.x[i] < 15 && (cf.a_rand_C[i] = 1.15 * rand(rng))
    end
    binmax = Dict{Float64,Float64}()                              # egen max(a_rand_C), by(CONVBIN)
    for i in eachindex(cf.x)
        binmax[CONVBIN[i]] = max(get(binmax, CONVBIN[i], 0.0), cf.a_rand_C[i])
    end
    for i in eachindex(cf.x)
        cf.a_rand_C[i] = binmax[CONVBIN[i]]
    end
    findeq!(cf, p)

    cbd(v, c) = maximum((isnan(v[i]) ? 0.0 : v[i]) for i in eachindex(c.x) if abs(c.x[i]) <= 5; init=0.0)
    sub(v, c) = maximum((isnan(v[i]) ? 0.0 : v[i]) for i in eachindex(c.x) if 11 <= c.x[i] <= 16; init=0.0)
    println("\n", "="^78)
    println("COUNTERFACTUAL -- commercial subcenter injected at 12-15 km (seeded; upstream is unseeded)")
    println("="^78)
    @printf("%-34s %14s %14s %12s\n", "", "inverted", "w/ subcenter", "% change")
    pct(u, v) = u == 0 ? NaN : 100 * (v / u - 1)
    for (lab, a, b) in (("max commercial height, |x|<=5 km", cbd(before.S_C, inv), cbd(cf.S_x_C, cf)),
                        ("max commercial height, 11-16 km", sub(before.S_C, inv), sub(cf.S_x_C, cf)),
                        ("wage y", before.sy, cf.sy),
                        ("employment N", before.sL, cf.sL),
                        ("urban radius x1 (km)", before.x1, cf.x1))
        @printf("%-34s %14.4f %14.4f %11.2f%%\n", lab, a, b, pct(a, b))
    end
    println("\n  Reallocation vs. growth: the original CBD loses height to the new subcenter, while")
    println("  the city as a whole expands because it has become more productive overall.")
    return cf
end

# ============================================================
# Visualization -- the three canonical gradient plots
# ============================================================

"""
    plot_gradients(c, p; tag="baseline")

`GHEIGHT` / `GBIDRENT` / `GLANDRENT` (`AB2022.ado:167-208`). Saves building height, floor-space
rent and land bid rent to graphs/ as PDF, with the urban area and CBD shaded.
"""
function plot_gradients(c::City, p::Params; tag::AbstractString="baseline")
    mkpath(GRAPHDIR)
    z(v) = [isnan(t) ? 0.0 : t for t in v]                 # Stata's `line` skips missing; area needs 0
    urban, com = z(c.URBAN), z(c.COM)
    shade(top) = (urban .* top, urban .* com .* top)
    made = String[]
    for (fname, yC, yR, ylab, ttl) in
            (("height",   c.S_x_C,     c.S_x_R,     "Building height",  "Building height"),
             ("floorrent", c.p_bar_x_C, c.p_bar_x_R, "Floor space rent", "Floor space rent"),
             ("landrent", c.r_x_C,     c.r_x_R,     "Land bid rent",    "Land bid rent"))
        top = maximum(t for t in vcat(z(yC), z(yR)) if isfinite(t); init=1.0) * 1.05
        su, sc = shade(top)
        plt = plot(c.x, su; seriestype=:sticks, linecolor=RGB(.88,.88,.88), label="Urban area",
                   legend=:topright, framestyle=:box, size=(620, 480), ylims=(0, top),
                   xlabel="Distance from centre (km)", ylabel=ylab, title=ttl,
                   xlims=(-1.6*c.x1, 1.6*c.x1))
        plot!(plt, c.x, sc; seriestype=:sticks, linecolor=RGB(.78,.78,.78), label="CBD")
        plot!(plt, c.x, yC; color=:red, ls=:dash, lw=2, label="Commercial")
        plot!(plt, c.x, yR; color=:blue, ls=:dot, lw=2, label="Residential")
        fname == "landrent" && plot!(plt, c.x, c.r_a; color=:black, lw=1.5, label="Agricultural")
        out = joinpath(GRAPHDIR, "AB2022_$(fname)_$(tag).pdf")
        savefig(plt, out); push!(made, basename(out))
    end
    return made
end

# ============================================================
# Main
# ============================================================

function main()
    gr()
    println("="^78)
    println("Ahlfeldt & Barr (2022), 'The economics of skyscrapers: A synthesis', JUE 129, 103419")
    println("Julia port of the AB2022 Stata toolkit -- endogenous building height")
    println("="^78)
    println("NOTATION WARNING: theta = construction-cost elasticity of height (0.5-0.55), NOT the")
    println("Frechet/trade elasticity (6-10); beta = AGGLOMERATION elasticity (+0.03), NOT the")
    println("amenity externality (-0.3). See the header before reusing any value elsewhere.")

    base_ok = validate_base()

    # Tier B under both conventions -- the decisive evidence on which one the data used.
    v_ado  = validate_inverted(convention=:ado)
    v_code = validate_inverted(convention=:codebook, quiet=true)
    println("\n", "="^78)
    println("WHICH omega CONVENTION GENERATED THE REFERENCE DATA?")
    println("="^78)
    @printf("  worst relative deviation, :ado      (1/(1-omega), AB2022.ado:97-98)      = %.3e\n", v_ado.worst)
    @printf("  worst relative deviation, :codebook (1/(1+omega), _1_PROGS.do:132-133)   = %.3e\n", v_code.worst)
    println("  => reference_data/INVERTED.csv was generated with the :ado convention, 1/(1-omega),")
    println("     even though the Codebook, _1_PROGS.do and the developer's FOC all say 1/(1+omega).")
    println("     The shipped default in this port remains :codebook (theoretically correct).")

    eq_ok = validate_equilibrium(v_ado)

    println("\n", "="^78); println("TIER D -- land-use boundaries"); println("="^78)
    @printf("  U identical at %d / %d rows (%s)\n", v_ado.nU_match, v_ado.n,
            v_ado.nU_match == v_ado.n ? "exact" : "MISMATCH")
    @printf("  CBD radius  x0 = %.2f km\n", v_ado.city.x0)
    @printf("  urban radius x1 = %.2f km\n", v_ado.city.x1)

    # ---- baseline equilibrium under the walkthrough figure parameters --------------------------
    p = params_analysis()
    println("\n", "="^78)
    println("BASELINE EQUILIBRIUM (_2_ANALYSIS.do:23 parameters, default :codebook convention)")
    println("="^78)
    city = build_city()
    res = findeq!(city, p; verbose=false)
    @printf("  converged in %d outer iterations, objective %.5f\n", res.iters, res.obj_ext)
    @printf("  Labour demand: %.2f\n  Labour supply: %.2f\n", city.L_hat_demand, city.L_hat_supply)
    @printf("  Total employment: %.0f\n  Wage: %.4f\n", city.sL, city.sy)
    @printf("  CBD radius x0: %.2f km\n  Urban radius x1: %.2f km\n", city.x0, city.x1)
    figs = plot_gradients(city, p; tag="baseline")
    println("  figures: ", join(figs, ", "))

    omega_convention_effect(p)
    _, hl = height_limit_counterfactual(p; S_bar=20.0)
    plot_gradients(hl, p; tag="heightlimit")

    # ---- Chicago skyline inversion -------------------------------------------------------------
    println("\n", "="^78)
    println("TIER E -- Chicago skyline inversion (_3_INVERSION.do)")
    println("="^78)
    pc = params_conv()
    inv_city = build_city()
    findeq!(inv_city, params_analysis())                    # _3_INVERSION.do:31 initialization
    HC, HR, CONVBIN = prepare_chicago(inv_city, joinpath(REFDIR, "EMPIRICAL_CH_skyline.csv"))
    @printf("  merged skyline: HEIGHT_C non-missing at %d rows, HEIGHT_R at %d rows\n",
            count(!isnan, HC), count(!isnan, HR))
    inv = invert!(inv_city, pc, 1.0e6, 0.05, HC, HR)
    @printf("\n  height correlation R2 = %.5f  (target >= 0.999, %s) after %d CONV steps\n",
            inv.r2, inv.r2 >= 0.999 ? "MET" : "NOT MET", inv.n_conv)
    @printf("  population gap = %.1f  (target <= 1000, %s) after %d EMP steps\n",
            inv.popgap, inv.popgap <= 1000 ? "MET" : "NOT MET", inv.n_emp)
    @printf("  equilibrium: N = %.0f, wage = %.4f, x0 = %.2f km, x1 = %.2f km\n",
            inv_city.sL, inv_city.sy, inv_city.x0, inv_city.x1)

    ref = CSV.read(joinpath(REFDIR, "INVERTED.csv"), DataFrame)
    for (lab, got, rf) in (("a_rand_C", inv_city.a_rand_C, nanvec(ref, :a_rand_C)),
                           ("a_rand_R", inv_city.a_rand_R, nanvec(ref, :a_rand_R)))
        w = [i for i in eachindex(inv_city.x) if abs(inv_city.x[i]) <= 20]
        g = got[w]
        @printf("  %s: zero share %.1f%% (port) vs %.1f%% (reference); mean|>0| %.4f vs %.4f\n",
                lab, 100*count(==(0.0), g)/length(g), 100*count(==(0.0), rf)/length(rf),
                mean(filter(>(0.0), g)), mean(filter(>(0.0), rf)))
    end
    println("\n  NOTE: the recovered amenities are compared IN DISTRIBUTION, not pointwise. The")
    println("  inversion path is as tolerance-limited as the equilibrium (Tier C), so an exact")
    println("  pointwise match is neither expected nor claimed.")
    plot_gradients(inv_city, pc; tag="chicago")

    cf = subcenter_counterfactual(inv_city, params_analysis(), CONVBIN)
    plot_gradients(cf, params_analysis(); tag="subcenter")

    println("\n", "="^78)
    println("SUMMARY: Tier A ", base_ok ? "PASS" : "FAIL",
            " | Tier B ", v_ado.ok ? "PASS" : "FAIL",
            " | Tier C ", eq_ok ? "PASS" : "FAIL",
            " | Tier D ", v_ado.nU_match == v_ado.n ? "PASS" : "FAIL",
            " | Tier E ", (inv.r2 >= 0.999 && inv.popgap <= 1000) ? "PASS" : "FAIL")
    println("Figures written to ", GRAPHDIR)
    println("="^78)
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    main()
end
