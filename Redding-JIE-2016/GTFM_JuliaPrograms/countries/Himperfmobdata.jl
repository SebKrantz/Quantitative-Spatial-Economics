# Monte Carlo for Quantitative Spatial Model;
# Helpman (increasing returns) model;
# Countries specification;

# SJR, November, 2015;

# Helpman counterpart of imperfmobdata.jl. The neoclassical equilibrium is
# solved first, its productivities and amenities are recovered under the
# Helpman model, and the Helpman economy is then solved and opened to trade.
# As in imperfmobdata.jl the transport network is held fixed, so only the
# trade weights (dclosed vs dopen) change. Ported from
# GTFM_MatlabPrograms/countries/Himperfmobdata.m; the figures are omitted.

# Run from this folder:  julia Himperfmobdata.jl

using Random
using Statistics
using LinearAlgebra
using StatsBase: geomean
using SpecialFunctions: gamma # gamma function

include("../graydist.jl")

# ************************
# **** Initialization ****
# ************************

# Set default random number stream
Random.seed!(1)

# *************************
# **** Distance matrix ****
# *************************

N = 11
NN = N * N

# Other latitude-longitude grid
ltd = range(0, stop=4, length=N)'
lgd = range(0, stop=4, length=N)

# Transport weights
tt0 = 7.9
tt1 = 1

# No change in the transport network: the corridor is present in both cases
tau0 = fill(tt0, (N, N))
tau0[6, :] .= tt1
tau0[:, 6] .= tt1

tau1 = fill(tt0, (N, N))
tau1[6, :] .= tt1
tau1[:, 6] .= tt1

dist0 = zeros(NN, NN)
dist1 = zeros(NN, NN)

for z in 1:NN
    seed = falses(N, N)
    seed[z] = true
    temp = graydist(tau0, seed, "quasi-euclidean")
    dist0[z, :] = vec(temp)
    temp = graydist(tau1, seed, "quasi-euclidean")
    dist1[z, :] = vec(temp)
end

for i in 1:size(dist0, 1)
    dist0[i, i] = 1
    dist1[i, i] = 1
end

rdist = dist1 ./ dist0

# Define treatment
treat = zeros(N, N)
treat[6, :] .= 1
treat[:, 6] .= 1
treat = vec(treat)

# Define east and west
Iwest = zeros(N, N)
Ieast = zeros(N, N)
Iwest[:, 1:7] .= 1
Ieast[:, 8:11] .= 1
Iwest = vec(Iwest)
Ieast = vec(Ieast)

# Define distance weights
dopen = ones(size(dist0))
dclosed = zeros(size(dist0))
dclosed[Iwest .== 1, Iwest .== 1] .= 1
dclosed[Ieast .== 1, Ieast .== 1] .= 1

# Trade costs are a power function of effective distance
dist0 = dist0 .^ 0.33
dist1 = dist1 .^ 0.33

# **************************
# **** Parameterization ****
# **************************

# Share of goods in consumption expenditure (1-housing share)
alpha = 0.75
# Elasticity of substitution
Hsigma = 5
# Goods Frechet shape parameter
theta = Hsigma - 1
# Worker Frechet shape parameter
epsilon = 3

param = [alpha, theta, epsilon]

# ***********************
# **** Random shocks ****
# ***********************

# Productivity and amenities are normalized within each country
a = exp.(randn(NN))
a[Iwest .== 1] = a[Iwest .== 1] ./ geomean(a[Iwest .== 1])
a[Ieast .== 1] = a[Ieast .== 1] ./ geomean(a[Ieast .== 1])

b = exp.(randn(NN))
b[Iwest .== 1] = b[Iwest .== 1] ./ geomean(b[Iwest .== 1])
b[Ieast .== 1] = b[Ieast .== 1] ./ geomean(b[Ieast .== 1])

@show mean(a), std(a), maximum(a), minimum(a)
@show mean(b), std(b), maximum(b), minimum(b)

# **************************
# **** Other Parameters ****
# **************************

# Observations
nobs = NN

# Land area
H = fill(100.0, nobs)

# Aggregate labor supply
LL = 153889 # US civilian labor force 2010 (Statistical Abstract, millions)
LLwest = (sum(Iwest) / (sum(Iwest) + sum(Ieast))) * LL
LLeast = (sum(Ieast) / (sum(Iwest) + sum(Ieast))) * LL

# Fixed production cost
F = 1

# ****************************************************************
# **** Closed Economy Solve for Endogenous Neoclassical Model ****
# ****************************************************************

fund = zeros(nobs, 5)
fund[:, 1] = a
fund[:, 2] = b
fund[:, 3] = H
fund[:, 4] = Iwest
fund[:, 5] = Ieast

include("functions/solveLwCtyClosed.jl")
w, L, tradesh, dtradesh, Lconverge, wconverge, xtic =
    solveLwCtyClosed(param, fund, dclosed, dist0, nobs)

@show "Wage and Population System Converged"
@show "Check Wage and Population Convergence"
@show wconverge, Lconverge
@show "Elapsed Time in Seconds"
@show xtic

# **************************************************************
# **** Closed Economy Solve for Unobservables Helpman model ****
# **************************************************************

observe = zeros(nobs, 5)
observe[:, 1] = L
observe[:, 2] = w
observe[:, 3] = H
observe[:, 4] = Iwest
observe[:, 5] = Ieast

include("functions/solveHabCtyClosed.jl")
Ha, Hb, abtradesh, aconverge, bconverge, xtic =
    solveHabCtyClosed(param, observe, dclosed, dist0, nobs)

@show "Productivity and Amenity System Converged"
@show "Check Productivity and Amenity Convergence"
@show aconverge, bconverge
@show "Elapsed Time in Seconds"
@show xtic

# *********************************************************************
# **** Closed Economy Solve for Endogenous variables Helpman model ****
# *********************************************************************

fund = zeros(nobs, 5)
fund[:, 1] = Ha
fund[:, 2] = Hb
fund[:, 3] = H
fund[:, 4] = Iwest
fund[:, 5] = Ieast

include("functions/solveHLwCtyClosed.jl")
w, L, tradesh, dtradesh, Lconverge, wconverge, xtic =
    solveHLwCtyClosed(param, fund, dclosed, dist0, nobs)

@show "Wage and Population System Converged"
@show "Check Wage and Population Convergence"
@show wconverge, Lconverge
@show "Elapsed Time in Seconds"
@show xtic

# Price index
include("functions/Hpindex.jl")
P = Hpindex(param, fund, L, w, dtradesh, nobs)

# Land price
include("functions/landprice.jl")
r = landprice(param, fund, L, w, dclosed, dist0, nobs)

# Expected utility
include("functions/Hexpectut.jl")
EU = Hexpectut(param, fund, L, w, P, r, dist0, nobs)
@show "Expected utility (West, East)"
@show unique(round.(EU .* 10.0^4) ./ 10.0^4)

# Welfare
include("functions/Hwelfare.jl")
welf = Hwelfare(param, fund, L, w, tradesh, dist0, nobs)
@show "Welfare"
welf = round.(welf .* 10.0^4) ./ 10.0^4
@show unique(welf)

# Real wage
include("functions/Hrealw.jl")
realwge = Hrealw(param, fund, L, w, tradesh, dist0, nobs)

# ************************************************
# **** Closed Economy Solve for Unobservables ****
# ************************************************

observe = zeros(nobs, 5)
observe[:, 1] = L
observe[:, 2] = w
observe[:, 3] = H
observe[:, 4] = Iwest
observe[:, 5] = Ieast

a_i, b_i, abtradesh, aconverge, bconverge, xtic =
    solveHabCtyClosed(param, observe, dclosed, dist0, nobs)

@show "Productivity and Amenity System Converged"
@show "Check Productivity and Amenity Convergence"
@show aconverge, bconverge
@show "Elapsed Time in Seconds"
@show xtic

# ***************************************************
# **** CLOSED ECONOMY AGGREGATE TO COUNTRY LEVEL ****
# ***************************************************

NOBSC = 2

WBILL = zeros(NOBSC)
WBILL[1] = sum(w[Iwest .== 1] .* L[Iwest .== 1])
WBILL[2] = sum(w[Ieast .== 1] .* L[Ieast .== 1])

OBSERVE = zeros(NOBSC, 5)
OBSERVE[1, 1] = sum(L[Iwest .== 1])
OBSERVE[2, 1] = sum(L[Ieast .== 1])
OBSERVE[1, 2] = WBILL[1] / OBSERVE[1, 1]
OBSERVE[2, 2] = WBILL[2] / OBSERVE[2, 1]
OBSERVE[1, 3] = sum(H[Iwest .== 1])
OBSERVE[2, 3] = sum(H[Ieast .== 1])
OBSERVE[1, 4] = 1
OBSERVE[2, 4] = 0
OBSERVE[1, 5] = 0
OBSERVE[2, 5] = 1

income = w .* L
trade = tradesh .* income'
TRADE = zeros(NOBSC, NOBSC)
TRADE[1, 1] = sum(trade[Iwest .== 1, Iwest .== 1])
TRADE[2, 2] = sum(trade[Ieast .== 1, Ieast .== 1])
TRADE[1, 2] = sum(trade[Iwest .== 1, Ieast .== 1])
TRADE[2, 1] = sum(trade[Ieast .== 1, Iwest .== 1])
EXPEND = vec(sum(TRADE, dims=2))
TRADESH = TRADE ./ EXPEND
DTRADESH = diag(TRADESH)

# *****************************************************
# **** Open Economy Solve for Endogenous Variables ****
# *****************************************************

include("functions/solveHLwCtyOpen.jl")
Cw, CL, Ctradesh, Cdtradesh, CLconverge, Cwconverge, Cxtic =
    solveHLwCtyOpen(param, fund, dopen, dist1, nobs)

@show "Wage and Population System Converged"
@show "Check Wage and Population Convergence"
@show Cwconverge, CLconverge
@show "Elapsed Time in Seconds"
@show Cxtic

# Counterfactual price index
CP = Hpindex(param, fund, CL, Cw, Cdtradesh, nobs)

# Counterfactual land prices
Cr = landprice(param, fund, CL, Cw, dopen, dist1, nobs)

# Counterfactual expected utility
CEU = Hexpectut(param, fund, CL, Cw, CP, Cr, dist1, nobs)
@show "Expected utility (West, East)"
@show unique(round.(CEU .* 10.0^4) ./ 10.0^4)

# Counterfactual welfare
Cwelf = Hwelfare(param, fund, CL, Cw, Ctradesh, dist1, nobs)
@show "Welfare"
Cwelf = round.(Cwelf .* 10.0^4) ./ 10.0^4
@show unique(Cwelf)

# Counterfactual real wage
Crealwge = Hrealw(param, fund, CL, Cw, Ctradesh, dist1, nobs)

# Welfare gains
include("functions/Hwelfaregains.jl")
welfgain = Hwelfaregains(param, Ctradesh, tradesh, CL, L, nobs)
@show "Welfare Gains"
welfgain = round.(welfgain .* 10.0^4) ./ 10.0^4
@show unique(welfgain)

# Perfectly immobile welfare gains (ACR)
include("functions/acrwelfaregains.jl")
acrwelfgain = acrwelfaregains(param, Ctradesh, tradesh, CL, L, nobs)

# Perfectly mobile welfare gains
include("functions/Hmobwelfaregains.jl")
mobwelfgain = Hmobwelfaregains(param, Ctradesh, tradesh, CL, L, nobs)

# ******************************************
# **** Open Economy Solve Unobservables ****
# ******************************************

Cobserve = zeros(nobs, 5)
Cobserve[:, 1] = CL
Cobserve[:, 2] = Cw
Cobserve[:, 3] = H
Cobserve[:, 4] = Iwest
Cobserve[:, 5] = Ieast

include("functions/solveHabCtyOpen.jl")
Ca_i, Cb_i, Cabtradesh, Caconverge, Cbconverge, Cxtic =
    solveHabCtyOpen(param, Cobserve, dopen, dist1, nobs)

@show "Productivity and Amenity System Converged"
@show "Check Productivity and Amenity Convergence"
@show Caconverge, Cbconverge
@show "Elapsed Time in Seconds"
@show Cxtic

# ************************************
# **** AGGREGATE TO COUNTRY LEVEL ****
# ************************************

CWBILL = zeros(NOBSC)
CWBILL[1] = sum(Cw[Iwest .== 1] .* CL[Iwest .== 1])
CWBILL[2] = sum(Cw[Ieast .== 1] .* CL[Ieast .== 1])

COBSERVE = zeros(NOBSC, 5)
COBSERVE[1, 1] = sum(CL[Iwest .== 1])
COBSERVE[2, 1] = sum(CL[Ieast .== 1])
COBSERVE[1, 2] = CWBILL[1] / COBSERVE[1, 1]
COBSERVE[2, 2] = CWBILL[2] / COBSERVE[2, 1]
COBSERVE[1, 3] = sum(H[Iwest .== 1])
COBSERVE[2, 3] = sum(H[Ieast .== 1])
COBSERVE[1, 4] = 1
COBSERVE[2, 4] = 0
COBSERVE[1, 5] = 0
COBSERVE[2, 5] = 1

Cincome = Cw .* CL
Ctrade = Ctradesh .* Cincome'
CTRADE = zeros(NOBSC, NOBSC)
CTRADE[1, 1] = sum(Ctrade[Iwest .== 1, Iwest .== 1])
CTRADE[2, 2] = sum(Ctrade[Ieast .== 1, Ieast .== 1])
CTRADE[1, 2] = sum(Ctrade[Iwest .== 1, Ieast .== 1])
CTRADE[2, 1] = sum(Ctrade[Ieast .== 1, Iwest .== 1])
CEXPEND = vec(sum(CTRADE, dims=2))
CTRADESH = CTRADE ./ CEXPEND
CDTRADESH = diag(CTRADESH)

# Welfare gains
WELFGAIN = Hwelfaregains(param, CTRADESH, TRADESH, COBSERVE[:, 1], OBSERVE[:, 1], NOBSC)
@show "Aggregate Country Welfare Gains"
WELFGAIN = round.(WELFGAIN .* 10.0^4) ./ 10.0^4
@show unique(WELFGAIN)

# Perfectly immobile welfare gains (ACR)
ACRWELFGAIN = acrwelfaregains(param, CTRADESH, TRADESH, COBSERVE[:, 1], OBSERVE[:, 1], NOBSC)

# Perfectly mobile welfare gains
MOBWELFGAIN = Hmobwelfaregains(param, CTRADESH, TRADESH, COBSERVE[:, 1], OBSERVE[:, 1], NOBSC)

# Compare region and country welfare gains
@show "Region and Country Welfare Gains [Region Country]"
temp = [unique(welfgain[Iwest .== 1]); unique(welfgain[Ieast .== 1])]
@show [temp WELFGAIN]

# Country perfectly immobile welfare gains
@show "Country Perfectly Immobile Welfare Gains"
@show unique(ACRWELFGAIN)

# Country perfectly mobile welfare gains
@show "Region and Country Perfectly Mobile Welfare Gains"
@show unique(MOBWELFGAIN)

# *************************
# **** Geometric Means ****
# *************************

gdtradesh = zeros(NOBSC)
gdtradesh[1] = geomean(dtradesh[Iwest .== 1])
gdtradesh[2] = geomean(dtradesh[Ieast .== 1])

Cgdtradesh = zeros(NOBSC)
Cgdtradesh[1] = geomean(Cdtradesh[Iwest .== 1])
Cgdtradesh[2] = geomean(Cdtradesh[Ieast .== 1])

gL = zeros(NOBSC)
gL[1] = geomean(L[Iwest .== 1])
gL[2] = geomean(L[Ieast .== 1])

gCL = zeros(NOBSC)
gCL[1] = geomean(CL[Iwest .== 1])
gCL[2] = geomean(CL[Ieast .== 1])

test = ((gdtradesh ./ Cgdtradesh) .^ (alpha / theta)) .*
       ((gL ./ gCL) .^ ((1 / epsilon) + (1 - alpha) - (alpha / theta)))

@show "Compare Aggregate and Geometric Mean Domestic Trade Share [Aggregate Geometric]"
@show [1 ./ CDTRADESH gdtradesh ./ Cgdtradesh]

@show "Compare Aggregate and Geometric Mean Labor Supply [Aggregate Geometric]"
@show [ones(NOBSC) gL ./ gCL]

@show "Compare Region Welfare, Test and Aggregate Welfare [region test aggregate]"
@show [temp test ACRWELFGAIN]
