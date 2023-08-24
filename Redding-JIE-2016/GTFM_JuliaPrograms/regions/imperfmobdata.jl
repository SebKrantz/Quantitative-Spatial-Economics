# Monte Carlo for Quantitative Spatial Model;
# Constant and increasing returns to scale model;
# Regions specification;

# SJR, November, 2015;

# *********************;
# **** Choose User ****;
# *********************;

using Random
using Statistics
using Plots
# using LaTeXStrings
using LinearAlgebra
using StatsBase: geomean
using SpecialFunctions: gamma # gamma function
# using Interpolations
using Dierckx # 2D interpolation
# using StatsModels: regress

cd("Redding-JIE-2016/GTFM_JuliaPrograms/regions")

include("../graydist.jl")

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
epsilon = 3

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

# ***********************************
# **** MEAN CHANGE IN TRADE COST ****
# ***********************************

mnrdist = mean(rdist, dims=2) 

mnrdmat = reshape(mnrdist, N, N)
XXL = range(minimum(lgd), stop=maximum(lgd), length=1000)  
YYL = range(minimum(ltd), stop=maximum(ltd), length=1000)' 

# https://stackoverflow.com/questions/41176260/julia-how-to-interpolate-non-uniformly-spaced-2-d-data-onto-a-grid
# using PyCall
# @pyimport scipy.interpolate as si
# using Dierckx
function my_spline2d(x, y, z, xx, yy)
    spl = Spline2D(collect(vec(x)), collect(vec(y)), collect(z))
    xx = collect(vec(xx))
    yy = collect(vec(yy))
    xxe = repeat(xx, inner = length(yy))
    yye = repeat(xx, outer = length(xx))
    zz = evaluate(spl, xxe, yye)
    return xx, yy, zz
end
# XXD, YYD, ZZD = griddata(lgd, ltd, mnrdmat, XXL, YYL)
XXD, YYD, ZZD = my_spline2d(lgd, ltd, mnrdmat, XXL, YYL)
p1 = contour(XXD, YYD, ZZD, levels=5, fill = true,
             title="Mean change in trade cost", 
             xlabel="Longitude", ylabel="Latitude") 
savefig(p1, "graphs/mean_rdist.pdf")

# ***********************************************
# **** Three-Dimensional Initial Equilibrium ****
# *********************************************** 

# PRODUCTIVITY
amat = reshape(a, N, N) 
XXL = range(minimum(lgd), stop=maximum(lgd), length=1000)
YYL = range(minimum(ltd), stop=maximum(ltd), length=1000)'

XXA, YYA, ZZA = my_spline2d(lgd, ltd, amat, XXL, YYL)  

# AMENITIES
bmat = reshape(b, N, N)
XXL = range(minimum(lgd), stop=maximum(lgd), length=1000)
YYL = range(minimum(ltd), stop=maximum(ltd), length=1000)' 

XXB, YYB, ZZB = my_spline2d(lgd, ltd, bmat, XXL, YYL)

# POPULATION  
Lmat = reshape(L, N, N)
XXL = range(minimum(lgd), stop=maximum(lgd), length=1000)
YYL = range(minimum(ltd), stop=maximum(ltd), length=1000)'

XXL, YYL, ZZL = my_spline2d(lgd, ltd, Lmat, XXL, YYL)

# PRICE INDEX
Pmat = reshape(P, N, N)  
XXP = range(minimum(lgd), stop=maximum(lgd), length=1000)
YYP = range(minimum(ltd), stop=maximum(ltd), length=1000)'

XXP, YYP, ZZP = my_spline2d(lgd, ltd, Pmat, XXP, YYP)

# WAGE
wmat = reshape(w, N, N)
XXw = range(minimum(lgd), stop=maximum(lgd), length=1000)
YYw = range(minimum(ltd), stop=maximum(ltd), length=1000)'

XXw, YYw, ZZw = my_spline2d(lgd, ltd, wmat, XXw, YYw)

# RELATIVE LAND PRICE
rmat = reshape(r, N, N)  
XXr = range(minimum(lgd), stop=maximum(lgd), length=1000)
YYr = range(minimum(ltd), stop=maximum(ltd), length=1000)'

XXr, YYr, ZZr = my_spline2d(lgd, ltd, rmat, XXr, YYr)

# MULTI-PANEL FIGURE
p2 = plot(
    contourf(XXA, YYA, ZZA, levels=5, title="Panel A: Productivity"),
        # xticks=(lgd[2], lgd[4], lgd[6], lgd[8], lgd[10]), 
        # yticks=(lgd[2], lgd[4], lgd[6], lgd[8], lgd[10]),
        # xticklabels=[2,4,6,8,10], yticklabels=[2,4,6,8,10]),
        
    contourf(XXB, YYB, ZZB, levels=5, title="Panel B: Amenities"),
    
    contourf(XXL, YYL, ZZL, levels=5, title="Panel C: Population"),
    
    contourf(XXP, YYP, ZZP, levels=5, title="Panel D: Price Index"),
    
    contourf(XXw, YYw, ZZw, levels=5, title="Panel E: Wages"),
    
    contourf(XXr, YYr, ZZr, levels=5, title="Panel F: Land Prices"),
    
    layout=(3, 2), size=(900, 800),
        titlefontsize = 13 # tickfontsize = 10
)
savefig(p2, "graphs/initial_equil.pdf")


# **************************************************************************
# **** Three-Dimensional Impact of Change in Trade Costs Counterfactual ****
# **************************************************************************

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

# GRID DATA

# POPULATION
dLmat = reshape(dL, N, N)
XXL = range(minimum(lgd), stop=maximum(lgd), length=1000)
YYL = range(minimum(ltd), stop=maximum(ltd), length=1000)'

XXL, YYL, ZZL = my_spline2d(lgd, ltd, dLmat, XXL, YYL)

# REAL WAGE  
drealwmat = reshape(drealw, N, N) 
XXrw = range(minimum(lgd), stop=maximum(lgd), length=1000)
YYrw = range(minimum(ltd), stop=maximum(ltd), length=1000)'

XXrw, YYrw, ZZrw = my_spline2d(lgd, ltd, drealwmat, XXrw, YYrw)

# PRICE INDEX
dPmat = reshape(dP, N, N)
XXP = range(minimum(lgd), stop=maximum(lgd), length=1000)
YYP = range(minimum(ltd), stop=maximum(ltd), length=1000)'

XXP, YYP, ZZP = my_spline2d(lgd, ltd, dPmat, XXP, YYP)

# WAGE
dwmat = reshape(dw, N, N) 
XXw = range(minimum(lgd), stop=maximum(lgd), length=1000)
YYw = range(minimum(ltd), stop=maximum(ltd), length=1000)' 

XXw, YYw, ZZw = my_spline2d(lgd, ltd, dwmat, XXw, YYw)

# RELATIVE LAND PRICE
drmat = reshape(dr, N, N)
XXr = range(minimum(lgd), stop=maximum(lgd), length=1000)  
YYr = range(minimum(ltd), stop=maximum(ltd), length=1000)'

XXr, YYr, ZZr = my_spline2d(lgd, ltd, drmat, XXr, YYr)

# ACR WELFARE
acrmat = reshape(acrwelfgain, N, N) 
XXa = range(minimum(lgd), stop=maximum(lgd), length=1000)
YYa = range(minimum(ltd), stop=maximum(ltd), length=1000)'

XXa, YYa, ZZa = my_spline2d(lgd, ltd, acrmat, XXa, YYa)

# WELFARE
welfmat = reshape(welfgain, N, N)
XXwelf = range(minimum(lgd), stop=maximum(lgd), length=1000) 
YYwelf = range(minimum(ltd), stop=maximum(ltd), length=1000)'

XXwelf, YYwelf, ZZwelf = my_spline2d(lgd, ltd, welfmat, XXwelf, YYwelf)

# MULTI-PANEL FIGURE
p3 = plot(
    contourf(XXL, YYL, ZZL, levels=5, title="Panel A: Population"),
    
    contourf(XXw, YYw, ZZw, levels=5, title="Panel B: Wages"),
    
    contourf(XXP, YYP, ZZP, levels=5, title="Panel C: Price Index"),
    
    contourf(XXr, YYr, ZZr, levels=5, title="Panel D: Land Rents"),
    
    contourf(XXrw, YYrw, ZZrw, levels=5, title="Panel E: Real Wage"),
    
    contourf(XXa, YYa, ZZa, levels=5, title="Panel F: Incorrect Immobile Welfare"),
        
    layout=(3, 2), size=(900, 800), titlefontsize = 12
)

savefig(p3, "graphs/transport_impact.pdf")


# HISTOGRAM FIGURE

# Population
p41 = histogram(ldL, bins=20, normalize=:probability, label="All")
histogram!(ldL[treat.==1], bins=20, normalize=:probability, label="Treated") 
histogram!(ldL[treat.==0], bins=20, normalize=:probability, label="Untreated")
title!("Panel A: Population")
ylabel!("Probability")

# Real wage
p42 = histogram(ldrealw, bins=20, normalize=:probability, label="All")
histogram!(ldrealw[treat.==1], bins=20, normalize=:probability, label="Treated")
histogram!(ldrealw[treat.==0], bins=20, normalize=:probability, label="Untreated") 
title!("Panel B: Real Wage")
ylabel!("Probability")

# Price index  
p43 = histogram(ldP, bins=20, normalize=:probability, label="All")
histogram!(ldP[treat.==1], bins=20, normalize=:probability, label="Treated")
histogram!(ldP[treat.==0], bins=20, normalize=:probability, label="Untreated")
title!("Panel C: Price Index") 
ylabel!("Probability")

# Wage
p44 = histogram(ldw, bins=20, normalize=:probability, label="All")  
histogram!(ldw[treat.==1], bins=20, normalize=:probability, label="Treated")
histogram!(ldw[treat.==0], bins=20, normalize=:probability, label="Untreated")
title!("Panel D: Wage")
ylabel!("Probability")

# Land price
p45 = histogram(ldr, bins=20, normalize=:probability, label="All")
histogram!(ldr[treat.==1], bins=20, normalize=:probability, label="Treated")
histogram!(ldr[treat.==0], bins=20, normalize=:probability, label="Untreated")
title!("Panel E: Land Rents")
ylabel!("Probability") 

# ACR welfare
p46 = histogram(acrwelfgain, bins=20, normalize=:probability, label="All") 
histogram!(acrwelfgain[treat.==1], bins=20, normalize=:probability, label="Treated")
histogram!(acrwelfgain[treat.==0], bins=20, normalize=:probability, label="Untreated")
title!("Panel F: Incorrect Immobile Welfare")
ylabel!("Probability")

plot(p41, p42, p43, p44, p45, p46, 
    layout=(3, 2), size=(900, 800), titlefontsize = 12) # , legend = :outertopright)
# Hard to display actual welfare gain due to bin definition;
# plot(welfgain * ones(1, 2), [0, 0.6], color="red", linestyle="-", linewidth=1.5)

savefig("graphs/transport_histogram.pdf")


# Perfect mobility welfare  
histogram(acrwelfgain, bins=20, normalize=:probability, label="ACR")
histogram!(mobwelfgain, bins=20, normalize=:probability, label="Perfect Mobility")
xlabel!("Relative Change in Welfare")
ylabel!("Probability")

savefig("graphs/perfmobil_histogram.pdf")


# HELPMAN GRID DATA

# HELPMAN POPULATION
HdLmat = reshape(HdL, N, N)  
XXL = range(minimum(lgd), stop=maximum(lgd), length=1000)
YYL = range(minimum(ltd), stop=maximum(ltd), length=1000)'

XXL, YYL, HZZL = my_spline2d(lgd, ltd, HdLmat, XXL, YYL)

# HELPMAN REAL WAGE
Hdrealwmat = reshape(Hdrealw, N, N)
XXrw = range(minimum(lgd), stop=maximum(lgd), length=1000)
YYrw = range(minimum(ltd), stop=maximum(ltd), length=1000)'

XXrw, YYrw, HZZrw = my_spline2d(lgd, ltd, Hdrealwmat, XXrw, YYrw)

# HELPMAN PRICE INDEX
HdPmat = reshape(HdP, N, N)
XXP = range(minimum(lgd), stop=maximum(lgd), length=1000)  
YYP = range(minimum(ltd), stop=maximum(ltd), length=1000)'

XXP, YYP, HZZP = my_spline2d(lgd, ltd, HdPmat, XXP, YYP)

# HELPMAN WAGE
Hdwmat = reshape(Hdw, N, N)
XXw = range(minimum(lgd), stop=maximum(lgd), length=1000)
YYw = range(minimum(ltd), stop=maximum(ltd), length=1000)'

XXw, YYw, HZZw = my_spline2d(lgd, ltd, Hdwmat, XXw, YYw)

# HELPMAN RELATIVE LAND PRICE  
Hdrmat = reshape(Hdr, N, N)
XXr = range(minimum(lgd), stop=maximum(lgd), length=1000)
YYr = range(minimum(ltd), stop=maximum(ltd), length=1000)'

XXr, YYr, HZZr = my_spline2d(lgd, ltd, Hdrmat, XXr, YYr)

# HELPMAN ACR WELFARE
Hacrmat = reshape(Hacrwelfgain, N, N)
XXa = range(minimum(lgd), stop=maximum(lgd), length=1000) 
YYa = range(minimum(ltd), stop=maximum(ltd), length=1000)'

XXa, YYa, HZZa = my_spline2d(lgd, ltd, Hacrmat, XXa, YYa)

# HELPMAN WELFARE
Hwelfmat = reshape(Hwelfgain, N, N)  
XXwelf = range(minimum(lgd), stop=maximum(lgd), length=1000)
YYwelf = range(minimum(ltd), stop=maximum(ltd), length=1000)'

XXwelf, YYwelf, HZZwelf = my_spline2d(lgd, ltd, Hwelfmat, XXwelf, YYwelf)

# HELPMAN HISTOGRAM FIGURE
p51 = histogram(lHdL, bins=20, normalize=:probability, label="All")
histogram!(lHdL[treat.==1], bins=20, normalize=:probability, label="Treated")
histogram!(lHdL[treat.==0], bins=20, normalize=:probability, label="Untreated")
title!("Panel A: Population")
ylabel!("Probability")

p52 = histogram(lHdrealw, bins=20, normalize=:probability, label="All")  
histogram!(lHdrealw[treat.==1], bins=20, normalize=:probability, label="Treated")
histogram!(lHdrealw[treat.==0], bins=20, normalize=:probability, label="Untreated")
title!("Panel B: Real Wage")
ylabel!("Probability")

p53 = histogram(lHdP, bins=20, normalize=:probability, label="All")
histogram!(lHdP[treat.==1], bins=20, normalize=:probability, label="Treated") 
histogram!(lHdP[treat.==0], bins=20, normalize=:probability, label="Untreated")
title!("Panel C: Price Index")
ylabel!("Probability")

p54 = histogram(lHdw, bins=20, normalize=:probability, label="All")
histogram!(lHdw[treat.==1], bins=20, normalize=:probability, label="Treated")
histogram!(lHdw[treat.==0], bins=20, normalize=:probability, label="Untreated")
title!("Panel D: Wage")  
ylabel!("Probability") 

p55 = histogram(lHdr, bins=20, normalize=:probability, label="All")
histogram!(lHdr[treat.==1], bins=20, normalize=:probability, label="Treated")
histogram!(lHdr[treat.==0], bins=20, normalize=:probability, label="Untreated")
title!("Panel E: Land Rents") 
ylabel!("Probability")

p56 = histogram(Hacrwelfgain, bins=20, normalize=:probability, label="All")
histogram!(Hacrwelfgain[treat.==1], bins=20, normalize=:probability, label="Treated") 
histogram!(Hacrwelfgain[treat.==0], bins=20, normalize=:probability, label="Untreated")
title!("Panel F: Incorrect Immobile Welfare")
ylabel!("Probability")  

plot(p51, p52, p53, p54, p55, p56, 
    layout=(3, 2), size=(900, 800), titlefontsize = 12) # , legend = :outertopright)
# Hard to display actual welfare gain due to bin definition;
# plot(Hwelfgain .* ones(2), [0, 0.6], color="red", linestyle="-", linewidth=1.5)

savefig("graphs/H_transport_histogram.pdf")


# HELPMAN PRODUCTIVITIES FIGURE 
p61 = scatter(log.(a), log.(Ha_i), marker=:circle, label="Productivity")
plot!(log.(a), log.(a), label="45 Degree Line")
xlabel!("Constant Returns")
ylabel!("Increasing Returns")
title!("Panel A: Log Productivity")

p62 = scatter(log.(b), log.(Hb_i), marker=:circle, label="Amenities") 
plot!(log.(b), log.(b), label="45 Degree Line")
xlabel!("Constant Returns") 
ylabel!("Increasing Returns")
title!("Panel B: Log Amenities")

plot(p61, p62, layout=(1, 2), size=(900, 300), titlefontsize = 12)

savefig("graphs/H_prod_amen.pdf")


# *******************************
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
@show sqrt.(diag(stats.covb))