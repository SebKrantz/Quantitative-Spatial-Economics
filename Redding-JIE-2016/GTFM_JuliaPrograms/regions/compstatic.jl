# Comparative statics for Quantitative Spatial Model;
# Constant returns to scale model;
# Regions specification;

# SJR, November, 2015;

# Sweeps the goods Frechet shape parameter theta and the worker Frechet
# shape parameter epsilon over a grid, re-solving the model and the
# transport counterfactual at every point, and reports how the treatment
# effects vary over that grid. Ported from
# GTFM_MatlabPrograms/regions/compstatic.m.

# Run from this folder:  julia compstatic.jl

using Random
using Statistics
using Plots
using LinearAlgebra
using StatsBase: geomean
using SpecialFunctions: gamma # gamma function
using Dierckx # 2D interpolation

# Work in this script's own folder, so that the graphs/ output paths resolve
# however the script is invoked
cd(@__DIR__)

include("../graydist.jl")
include("functions/solveLw.jl")
include("functions/pindex.jl")
include("functions/landprice.jl")
include("functions/welfare.jl")
include("functions/realw.jl")
include("functions/welfaregains.jl")
include("functions/acrwelfaregains.jl")
# note: functions/obwelfaregains.jl defines a second, observation-weighted
# method under the same name; this is the one imperfmobdata.jl uses
include("functions/mobwelfaregains.jl")
include("functions/regress.jl")

# ************************
# **** Initialization ****
# ************************

# Set default random number stream
Random.seed!(1)

# *********************************************
# **** Trade cost matrix (called distance) ****
# *********************************************

N = 11
NN = N * N

# Other latitude-longitude grid
ltd = range(0, stop=4, length=N)'
lgd = range(0, stop=4, length=N)

# Transport weights
tt0 = 7.9
tt1 = 1

tau0 = fill(tt0, (N, N))
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

# Define controls
X = [ones(NN) treat]

# Trade costs are a power function of effective distance
dist0 = dist0 .^ 0.33
dist1 = dist1 .^ 0.33

# **************************
# **** Parameterization ****
# **************************

# Share of goods in consumption expenditure (1-housing share)
alpha = 0.75
# Elasticity of substitution
sigma = 4

# ***********************
# **** Random shocks ****
# ***********************

a = exp.(randn(NN))
a = a ./ geomean(a)

b = exp.(randn(NN))
b = b ./ geomean(b)

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

fund = zeros(nobs, 5)
fund[:, 1] = a
fund[:, 2] = b
fund[:, 3] = H

# ************************************
# **** Loop over parameter values ****
# ************************************

# Parameter grid
# Sigma is 4, so need theta greater than 3
epsgrid = range(3.1, stop=5.1, length=10)
thetagrid = range(3.1, stop=5.1, length=10)
K = length(epsgrid)
KK = K * K

# ndgrid(thetagrid, epsgrid) read out in column-major order: theta varies
# fastest, epsilon slowest
thetavec = repeat(collect(thetagrid), outer=K)
epsvec = repeat(collect(epsgrid), inner=K)

# Matrices to store
dLmat = zeros(KK)
drwmat = zeros(KK)
dwmat = zeros(KK)
drmat = zeros(KK)
dPmat = zeros(KK)
dwelfmat = zeros(KK)
dacrmat = zeros(KK)
dmobmat = zeros(KK)

convmat = zeros(KK, 4)

for p in 1:KK

    epsilon = epsvec[p]
    theta = thetavec[p]
    estparam = [alpha, theta, epsilon]
    @show "Value of parameters are [alpha theta epsilon]", estparam

    # ****************************************
    # **** Solve for Endogenous Variables ****
    # ****************************************

    # Solve for region populations and wages
    w, L, tradesh, dtradesh, Lconverge, wconverge, xtic =
        solveLw(estparam, fund, dist0, nobs)
    convmat[p, 1] = wconverge
    convmat[p, 2] = Lconverge

    # Price index
    P = pindex(estparam, fund, w, dtradesh)

    # Land price
    r = landprice(estparam, fund, L, w)

    # Welfare
    welf = welfare(estparam, fund, L, tradesh)

    # Real wage
    realwage = realw(estparam, fund, L, tradesh)

    # *******************************
    # **** CHANGE IN TRADE COSTS ****
    # *******************************

    # Solve for region populations and wages
    Cw, CL, Ctradesh, Cdtradesh, CLconverge, Cwconverge, Cxtic =
        solveLw(estparam, fund, dist1, nobs)
    convmat[p, 3] = Cwconverge
    convmat[p, 4] = CLconverge

    # Counterfactual price index
    CP = pindex(estparam, fund, Cw, Cdtradesh)

    # Counterfactual land prices
    Cr = landprice(estparam, fund, CL, Cw)

    # Counterfactual welfare
    Cwelf = welfare(estparam, fund, CL, Ctradesh)

    # Counterfactual real wage
    Crealwage = realw(estparam, fund, CL, Ctradesh)

    # Welfare gains
    welfgain = welfaregains(estparam, Ctradesh, tradesh, CL, L)
    welfgain = round.(welfgain .* 10.0^4) ./ 10.0^4

    # Perfectly immobile welfare gains (ACR)
    acrwelfgain = acrwelfaregains(estparam, Ctradesh, tradesh)

    # Perfectly mobile welfare gains
    mobwelfgain = mobwelfaregains(estparam, Ctradesh, tradesh, CL, L)

    # **************************
    # **** Relative Changes ****
    # **************************

    dL = CL ./ L
    ldL = log.(dL)
    dw = Cw ./ w
    ldw = log.(dw)
    dr = Cr ./ r
    ldr = log.(dr)
    dP = CP ./ P
    ldP = log.(dP)
    lacrwelfgain = log.(acrwelfgain)
    lmobwelfgain = log.(mobwelfgain)
    drealw = Crealwage ./ realwage
    ldrealw = log.(drealw)

    # Population treatment
    bL = regress(ldL, X)[1]
    dLmat[p] = bL[2]
    # Wage treatment
    bw = regress(ldw, X)[1]
    dwmat[p] = bw[2]
    # Price index treatment
    bP = regress(ldP, X)[1]
    dPmat[p] = bP[2]
    # Land price treatment
    br = regress(ldr, X)[1]
    drmat[p] = br[2]
    # Real wage treatment
    brw = regress(ldrealw, X)[1]
    drwmat[p] = brw[2]
    # ACR welfare treatment
    bacr = regress(lacrwelfgain, X)[1]
    dacrmat[p] = bacr[2]
    # Mobility welfare treatment
    bmob = regress(lmobwelfgain, X)[1]
    dmobmat[p] = bmob[2]

end

# ****************************************************
# **** Check all parameter combinations converged ****
# ****************************************************

@show "Check all parameter combinations converged"
@show minimum(convmat)

# ********************************************************
# **** Analyze how effects vary with parameter values ****
# ********************************************************

# z is indexed [theta, epsilon]; transpose it so that the interpolated
# surface is indexed [epsilon, theta], which is what contourf(theta, eps, z)
# expects
function my_spline2d(x, y, z, xx, yy)
    spl = Spline2D(collect(vec(x)), collect(vec(y)), collect(transpose(z)))
    xx = collect(vec(xx))
    yy = collect(vec(yy))
    xxe = repeat(xx, inner = length(yy))
    yye = repeat(yy, outer = length(xx))
    zz = reshape(evaluate(spl, xxe, yye), length(yy), length(xx))
    return xx, yy, zz
end

# Fine grid for the contour surfaces
TT = range(minimum(thetagrid), stop=maximum(thetagrid), length=200)
EE = range(minimum(epsgrid), stop=maximum(epsgrid), length=200)

# Each *mat is stored over p with theta varying fastest
dLmat = reshape(dLmat, K, K)
drwmat = reshape(drwmat, K, K)
dPmat = reshape(dPmat, K, K)
dwmat = reshape(dwmat, K, K)
drmat = reshape(drmat, K, K)
dacrmat = reshape(dacrmat, K, K)
dmobmat = reshape(dmobmat, K, K)

XXL, YYL, ZZL = my_spline2d(thetagrid, epsgrid, permutedims(dLmat), TT, EE)
XXw, YYw, ZZw = my_spline2d(thetagrid, epsgrid, permutedims(dwmat), TT, EE)
XXP, YYP, ZZP = my_spline2d(thetagrid, epsgrid, permutedims(dPmat), TT, EE)
XXr, YYr, ZZr = my_spline2d(thetagrid, epsgrid, permutedims(drmat), TT, EE)
XXrw, YYrw, ZZrw = my_spline2d(thetagrid, epsgrid, permutedims(drwmat), TT, EE)
XXacrb, YYacrb, ZZacrb = my_spline2d(thetagrid, epsgrid, permutedims(dacrmat), TT, EE)
XXmobb, YYmobb, ZZmobb = my_spline2d(thetagrid, epsgrid, permutedims(dmobmat), TT, EE)

# Multi-panel figure
p1 = contourf(XXL, YYL, ZZL, levels=10, title="Panel A : Population Treatment")
p2 = contourf(XXw, YYw, ZZw, levels=10, title="Panel B : Wage Treatment")
p3 = contourf(XXP, YYP, ZZP, levels=10, title="Panel C : Price Treatment")
p4 = contourf(XXr, YYr, ZZr, levels=10, title="Panel D : Land Rent Treatment")
p5 = contourf(XXrw, YYrw, ZZrw, levels=10, title="Panel E : Real Wage Treatment")
p6 = contourf(XXacrb, YYacrb, ZZacrb, levels=10, title="Panel F : Incorrect Immobile Welfare")

pall = plot(p1, p2, p3, p4, p5, p6, layout=(3, 2), size=(900, 1000),
            xlabel="Theta", ylabel="Epsilon", titlefontsize=8,
            guidefontsize=8, tickfontsize=6)
savefig(pall, "graphs/compstatic.pdf")

# Mobility welfare bias
pmob = contourf(XXmobb, YYmobb, ZZmobb, levels=10, title="Mobility Welfare Bias",
                xlabel="Theta", ylabel="Epsilon", titlefontsize=8,
                guidefontsize=8, tickfontsize=6)
savefig(pmob, "graphs/compstatic_mobbias.pdf")
