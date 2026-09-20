# Monte Carlo for Quantitative Spatial Model;
# Constant and increasing returns to scale model;
# Regions specification;

# SJR, November, 2015;

# Robustness rerun of imperfmobdata.jl holding the worker Frechet shape
# parameter fixed at epsilon = 5.1. Same specification and the same
# treatment regressions; the figures are omitted, as in the MATLAB original
# (GTFM_MatlabPrograms/regions/imperfmobdata_eps5.m).

# Run from this folder:  julia imperfmobdata_eps5.jl

using Random
using Statistics
using LinearAlgebra
using StatsBase: geomean
using SpecialFunctions: gamma # gamma function

include("../graydist.jl")
include("functions/regress.jl")
include("functions/regstats.jl")

# ************************
# **** Initialization ****
# ************************

# Set default random number stream
Random.seed!(1) 

# *********************************************
# **** Trade cost matrix (called distance) ****
# *********************************************

N = 11
NN = N*N

# Other latitude-longitude grid
ltd = range(0, stop=4, length=N)' 
lgd = range(0, stop=4, length=N)

# Transport weights
tt0 = 7.9
tt1 = 1

tau0 = fill(tt0, (N,N))
tau1 = fill(tt0, (N,N))  
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

for i in 1:size(dist0,1)
    dist0[i,i] = 1
    dist1[i,i] = 1
end

rdist = dist1 ./ dist0

# Define treatment
treat = zeros(N,N)
treat[6, :] .= 1
treat[:, 6] .= 1 
treat = vec(treat)

# Trade costs are a power function of effective distance  
dist0 = dist0.^0.33
dist1 = dist1.^0.33

# ************************** 
# **** Parameterization ****
# **************************

# Share of goods in consumption expenditure (1-housing share)
alpha = 0.75
# Elasticity of substitution
sigma = 4  
Hsigma = 5
# Goods Frechet shape parameter
theta = 4
# Worker Frechet shape parameter  
epsilon = 5.1

param = [alpha, theta, epsilon]

# ***********************
# **** Random shocks ****
# ***********************

a = exp.(randn(NN, 1)) 
a ./= geomean(a)

b = exp.(randn(NN, 1))
b ./= geomean(b)

@show mean(a), std(a), maximum(a), minimum(a)
@show mean(b), std(b), maximum(b), minimum(b)

# **************************
# **** Other Parameters ****
# **************************

# Observations
nobs = NN  

# Land area
H = fill(100.0, nobs) 

# Aggregate labor Supply  
LL = 153889 # US civilian labor force 2010 (Statistical Abstract, millions)

# Fixed production cost
F = 1

# ****************************************  
# **** Solve for Endogenous Variables ****
# ****************************************

fund = zeros(nobs, 5)
fund[:, 1] = a
fund[:, 2] = b
fund[:, 3] = H

# Solve for region populations and wages
include("functions/solveLw.jl")
w, L, tradesh, dtradesh, Lconverge, wconverge, xtic = 
    solveLw(param, fund, dist0, nobs)

@show "Wage and Population System Converged"
@show "Check Wage and Population Convergence"
@show wconverge, Lconverge 
@show "Elapsed Time in Seconds"
@show xtic 


# Price index
include("functions/pindex.jl")
P = pindex(param, fund, w, dtradesh)

# Land price
include("functions/landprice.jl")
r = landprice(param, fund, L, w)

# Expected utility
include("functions/expectut.jl")
EU = expectut(param, fund, L, tradesh)
@show "Expected utility"
@show EU

# Welfare  
include("functions/welfare.jl")
welf = welfare(param, fund, L, tradesh)
@show "Welfare"
welf = round.(welf .* 10.0^4) ./ 10.0^4
@show unique(welf) 

# Real wage
include("functions/realw.jl")
realwage = realw(param, fund, L, tradesh)

# *********************************
# **** Solve for Unobservables ****
# *********************************

observe = zeros(nobs, 5)
observe[:, 1] = L
observe[:, 2] = w
observe[:, 3] = H

# Solve for region productivities and amenities  
include("functions/solveab.jl")
a_i, b_i, abtradesh, aconverge, bconverge, xtic = 
    solveab(param, observe, dist0, nobs)

@show "Productivity and Amenity System Converged"
@show "Check Productivity and Amenity Convergence" 
@show aconverge, bconverge
@show "Elapsed Time in Seconds" 
@show xtic

# *****************************************
# **** Solve for Helpman Unobservables ****
# *****************************************

# Solve for region productivities and amenities
include("functions/solveHab.jl")
Ha_i, Hb_i, Habtradesh, Haconverge, Hbconverge, xtic = 
    solveHab(param, observe, dist0, nobs)
    
@show "Helpman Productivity and Amenity System Converged"
@show "Check Productivity and Amenity Convergence"
@show Haconverge, Hbconverge
@show "Elapsed Time in Seconds"
@show xtic

Hfund = zeros(nobs, 5) 
Hfund[:, 1] = Ha_i
Hfund[:, 2] = Hb_i
Hfund[:, 3] = H

# ************************************************
# **** Solve for Helpman Endogenous Variables ****
# ************************************************

# Solve for region populations and wages
include("functions/solveHLw.jl")
Hw, HL, Htradesh, Hdtradesh, HLconverge, Hwconverge, Hxtic = 
    solveHLw(param, Hfund, dist0, nobs)
    
@show "Wage and Population System Converged"
@show "Check Wage and Population Convergence"  
@show Hwconverge, HLconverge
@show "Elapsed Time in Seconds" 
@show xtic

# Price index
include("functions/Hpindex.jl")
HP = Hpindex(param, Hfund, HL, Hw, Hdtradesh)

# Land price
Hr = landprice(param, Hfund, HL, Hw)  

# Expected utility
include("functions/Hexpectut.jl")
HEU = Hexpectut(param, Hfund, Hw, HP, Hr)
@show "Helpman Expected utility"
@show HEU

# Welfare
include("functions/Hwelfare.jl")
Hwelf = Hwelfare(param, Hfund, HL, Htradesh) 
@show "Helpman Welfare"
Hwelf = round.(Hwelf .* 10.0^4) ./ 10.0^4
@show unique(Hwelf)

# Real wage  
include("functions/Hrealw.jl")
Hrealwage = Hrealw(param, Hfund, HL, Htradesh)

# *******************************
# **** CHANGE IN TRADE COSTS ****
# *******************************

# Solve for region populations and wages
Cw, CL, Ctradesh, Cdtradesh, CLconverge, Cwconverge, Cxtic = 
    solveLw(param, fund, dist1, nobs)
    
@show "Wage and Population System Converged"
@show "Check Wage and Population Convergence"
@show Cwconverge, CLconverge
@show "Elapsed Time in Seconds"
@show Cxtic

# Counterfactual price index
CP = pindex(param, fund, Cw, Cdtradesh)  

# Counterfactual land prices
Cr = landprice(param, fund, CL, Cw)

# Counterfactual expected utility
CEU = expectut(param, fund, CL, Ctradesh)  
@show "Expected utility"
@show CEU

# Counterfactual welfare
Cwelf = welfare(param, fund, CL, Ctradesh)
@show "Welfare"
Cwelf = round.(Cwelf .* 10.0^4) ./ 10.0^4  
@show unique(Cwelf)

# Counterfactual real wage
Crealwage = realw(param, fund, CL, Ctradesh) 

# Welfare gains
include("functions/welfaregains.jl")
welfgain = welfaregains(param, Ctradesh, tradesh, CL, L)
@show "Welfare Gains"
welfgain = round.(welfgain .* 10.0^4) ./ 10.0^4
@show unique(welfgain)

# Perfectly immobile welfare gains (ACR)
include("functions/acrwelfaregains.jl")
acrwelfgain = acrwelfaregains(param, Ctradesh, tradesh)

# Perfectly mobile no preference heterogeneity welfare gains  
include("functions/mobwelfaregains.jl")
mobwelfgain = mobwelfaregains(param, Ctradesh, tradesh, CL, L)


# ***************************************
# **** HELPMAN CHANGE IN TRADE COSTS ****
# ***************************************

# Solve for region populations and wages
CHw, CHL, CHtradesh, CHdtradesh, CHLconverge, CHwconverge, CHxtic = 
    solveHLw(param, Hfund, dist1, nobs)

@show "Wage and Population System Converged"
@show "Check Wage and Population Convergence"
@show CHwconverge, CHLconverge
@show "Elapsed Time in Seconds"
@show CHxtic

# Counterfactual price index 
CHP = Hpindex(param, Hfund, CHL, CHw, CHdtradesh)

# Counterfactual land prices
CHr = landprice(param, Hfund, CHL, CHw)

# Counterfactual expected utility
CHEU = Hexpectut(param, Hfund, CHw, CHP, CHr)
@show "Helpman Expected utility"
@show CHEU

# Counterfactual welfare
CHwelf = Hwelfare(param, Hfund, CHL, CHtradesh)
@show "Helpman Welfare" 
CHwelf = round.(CHwelf .* 10.0^4) ./ 10.0^4
@show unique(CHwelf)

# Counterfactual real wage
CHrealwage = Hrealw(param, Hfund, CHL, CHtradesh)

# Welfare gains
include("functions/Hwelfaregains.jl")
Hwelfgain = Hwelfaregains(param, CHtradesh, tradesh, CHL, L)
@show "Helpman Welfare Gains"
Hwelfgain = round.(Hwelfgain .* 10.0^4) ./ 10.0^4  
@show unique(Hwelfgain)

# Perfectly immobile welfare gains (ACR)
Hacrwelfgain = acrwelfaregains(param, CHtradesh, tradesh)

# Perfectly mobile welfare gains
Hmobwelfgain = mobwelfaregains(param, CHtradesh, tradesh, CHL, L)

# *****************************
# **** Solve Unobservables ****
# *****************************

Cobserve = zeros(nobs, 5)
Cobserve[:, 1] = CL
Cobserve[:, 2] = Cw
Cobserve[:, 3] = H

# Solve for region productivities and amenities
include("functions/solveab.jl")
Ca_i, Cb_i, Cabtradesh, Caconverge, Cbconverge, Cxtic = 
    solveab(param, Cobserve, dist1, nobs)
    
@show "Productivity and Amenity System Converged" 
@show "Check Productivity and Amenity Convergence"
@show Caconverge, Cbconverge
@show "Elapsed Time in Seconds"
@show Cxtic 

# *************************************
# **** Solve Helpman Unobservables ****  
# *************************************

CHobserve = zeros(nobs, 5) 
CHobserve[:, 1] = CHL
CHobserve[:, 2] = CHw
CHobserve[:, 3] = H

# Solve for region productivities and amenities
CHa_i, CHb_i, CHabtradesh, CHaconverge, CHbconverge, CHxtic = 
    solveHab(param, CHobserve, dist1, nobs)
# TODO: seems it has not converged for amenities    
@show "Helpman Productivity and Amenity System Converged"
@show "Check Productivity and Amenity Convergence"
@show CHaconverge, CHbconverge  
@show "Elapsed Time in Seconds"
@show CHxtic


# *****************
# **** Impacts ****
# *****************

dL = CL ./ L
ldL = log.(dL)  

dw = Cw ./ w
ldw = log.(dw)

dr = Cr ./ r
ldr = log.(dr)

dP = CP ./ P
ldP = log.(dP)

lacrwelfgain = log.(acrwelfgain)

drealw = Crealwage ./ realwage  
ldrealw = log.(drealw)

dtradesh = diag(tradesh) 
Cdtradesh = diag(Ctradesh)
ddtradesh = Cdtradesh ./ dtradesh

HdL = CHL ./ HL
lHdL = log.(HdL)

Hdw = CHw ./ Hw
lHdw = log.(Hdw) 

Hdr = CHr ./ Hr
lHdr = log.(Hdr)

HdP = CHP ./ HP
lHdP = log.(HdP)

lHacrwelfgain = log.(Hacrwelfgain)

Hdrealw = CHrealwage ./ Hrealwage
lHdrealw = log.(Hdrealw)

CHdtradesh = diag(CHtradesh)
Hddtradesh = CHdtradesh ./ dtradesh

# **** Treatment Regressions ****
# *******************************

@show "Treatment Regressions"

# Define controls
X = [ones(size(treat)) treat]   

# Population treatment
bL, bintL, rL, rintL, statsL = regress(ldL, X)  
@show "Population Treatment"
@show bL[2] 
stats = regstats(ldL, treat, "linear", ["beta", "covb"])
@show stats.beta
@show sqrt.(diag(stats.covb))

# Real wage treatment
brw, bintrw, rrw, rintrw, statsrw = regress(ldrealw, X)
@show "Real Wage Treatment"
@show brw[2]
stats = regstats(ldrealw, treat, "linear", ["beta", "covb"])
@show stats.beta
@show sqrt.(diag(stats.covb))  

# Price index treatment
bP, bintP, rP, rintP, statsP = regress(ldP, X)
@show "Price Index Treatment"  
@show bP[2]
stats = regstats(ldP, treat, "linear", ["beta", "covb"])
@show stats.beta
@show sqrt.(diag(stats.covb))

# Wage treatment
bw, bintw, rw, rintw, statsw = regress(ldw, X)  
@show "Wage Treatment"
@show bw[2]
stats = regstats(ldw, treat, "linear", ["beta", "covb"])
@show stats.beta
@show sqrt.(diag(stats.covb))

# Land price treatment
br, bintr, rr, rintr, statsr = regress(ldr, X)
@show "Land Price Treatment"
@show br[2]
stats = regstats(ldr, treat, "linear", ["beta", "covb"]) 
@show stats.beta
@show sqrt.(diag(stats.covb))

# ACR welfare
bacr, bintacr, racr, rintacr, statsacr = regress(lacrwelfgain, X) 
@show "ACR Welfare Treatment"
@show bacr[2]
stats = regstats(lacrwelfgain, treat, "linear", ["beta", "covb"])
@show stats.beta
@show sqrt.(diag(stats.covb))

# ***************************************
# **** Helpman Treatment Regressions ****
# ***************************************

@show "Helpman Treatment Regressions" 

# Define controls
X = [ones(size(treat)) treat]

# Population treatment
bL, bintL, rL, rintL, statsL = regress(lHdL, X)
@show "Population Treatment"
@show bL[2]
stats = regstats(lHdL, treat, "linear", ["beta", "covb"])
@show stats.beta 
@show sqrt.(diag(stats.covb))

# Real wage treatment
brw, bintrw, rrw, rintrw, statsrw = regress(lHdrealw, X) 
@show "Real Wage Treatment"
@show brw[2]
stats = regstats(lHdrealw, treat, "linear", ["beta", "covb"])
@show stats.beta
@show sqrt.(diag(stats.covb))

# Price index treatment
bP, bintP, rP, rintP, statsP = regress(lHdP, X)
@show "Price Index Treatment"
@show bP[2]
stats = regstats(lHdP, treat, "linear", ["beta", "covb"])
@show stats.beta
@show sqrt.(diag(stats.covb))

# Wage treatment
bw, bintw, rw, rintw, statsw = regress(lHdw, X)
@show "Wage Treatment" 
@show bw[2]
stats = regstats(lHdw, treat, "linear", ["beta", "covb"])
@show stats.beta
@show sqrt.(diag(stats.covb))

# Land price treatment
br, bintr, rr, rintr, statsr = regress(lHdr, X)  
@show "Land Price Treatment"
@show br[2]
stats = regstats(lHdr, treat, "linear", ["beta", "covb"])
@show stats.beta
@show sqrt.(diag(stats.covb))

# ACR welfare
bacr, bintacr, racr, rintacr, statsacr = regress(lHacrwelfgain, X)
@show "ACR Welfare Treatment"
@show bacr[2] 
stats = regstats(lHacrwelfgain, treat, "linear", ["beta", "covb"])
@show stats.beta
