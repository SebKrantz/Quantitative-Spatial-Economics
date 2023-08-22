using Random, LinearAlgebra, Statistics, Plots

# Set default random number stream
Random.seed!(1)

# Distance matrix
N = 11
NN = N * N
ltd = range(0, stop=4, length=N)
lgd = range(0, stop=4, length=N)

# Transport weights
tt0 = 7.9
tt1 = 1

tau0 = fill(tt0, N, N)
tau1 = fill(tt0, N, N)
tau1[:, 6] = tt1
tau1[6, :] = tt1

dist0 = zeros(NN)
dist1 = zeros(NN)

for z in 1:NN
    seed = falses(N, N)
    seed[z] = true
    temp = graydist(tau0, seed, "quasi-euclidean")
    dist0[z, :] = reshape(temp, 1, NN)
    temp = graydist(tau1, seed, "quasi-euclidean")
    dist1[z, :] = reshape(temp, 1, NN)
end
dist0[diagind(dist0)] .= 1
dist1[diagind(dist1)] .= 1
rdist = dist1 ./ dist0

# Define treatment
treat = zeros(N, N)
treat[:, 6] = 1
treat[6, :] = 1
treat = reshape(treat, NN)

# Define controls
X = [ones(NN) treat]

# Trade costs are a power function of effective distance
dist0 = dist0 .^ 0.33
dist1 = dist1 .^ 0.33

# Parameters
alpha = 0.75
sigma = 4

# Random shocks
a = randn(NN)
a = exp.(a)
a = a ./ geomean(a)

b = randn(NN)
b = exp.(b)
b = b ./ geomean(b)

# Other parameters
nobs = NN
H = fill(100, nobs)
LL = 153889

# Loop over parameter values
epsgrid = range(3.1, stop=5.1, length=10)
thetagrid = range(3.1, stop=5.1, length=10)
K = length(epsgrid)
KK = K * K
thetamat, epsmat = ndgrid(thetagrid, epsgrid)
epsvec = reshape(epsmat, KK)
thetavec = reshape(thetamat, KK)

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

    # Solve for region populations and wages
    w, L, tradesh, dtradesh, Lconverge, wconverge, xtic = solveLw(estparam, fund, dist0, nobs)
    convmat[p, 1] = wconverge
    convmat[p, 2] = Lconverge

    # Price index
    P = pindex(estparam, fund, w, dtradesh, nobs)

    # Land price
    r = landprice(estparam, fund, L, w, dist0, nobs)

    # Welfare
    welf = welfare(estparam, fund, L, w, tradesh, dist0, nobs)

    # Real wage
    realwage = realw(estparam, fund, L, w, tradesh, dist0, nobs)

    # Change in trade costs
    Cw, CL, Ctradesh, Cdtradesh, CLconverge, Cwconverge, Cxtic = solveLw(estparam, fund, dist1, nobs)
    convmat[p, 3] = Cwconverge
    convmat[p, 4] = CLconverge

    # Counterfactual price index
    CP = pindex(estparam, fund, Cw, Cdtradesh, nobs)

    # Counterfactual land prices
    Cr = landprice(estparam, fund, CL, Cw, dist1, nobs)

    # Counterfactual welfare
    Cwelf = welfare(estparam, fund, CL, Cw, Ctradesh, dist1, nobs)

    # Counterfactual real wage
    Crealwage = realw(estparam, fund, CL, Cw, Ctradesh, dist1, nobs)

    # Welfare gains
    welfgain = welfaregains(estparam, fund, Ctradesh, tradesh, CL, L, nobs)
    welfgain = round.(welfgain .* 10^4) ./ 10^4

    # Perfectly immobile welfare gains (ACR)
    acrwelfgain = acrwelfaregains(estparam, fund, Ctradesh, tradesh, CL, L, nobs)

    # Perfectly mobile welfare gains
    mobwelfgain = mobwelfaregains(estparam, fund, Ctradesh, tradesh, CL, L, nobs)

    # Relative changes
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

# Check all parameter combinations converged
minimum(convmat)

# Analyze how effects vary with parameter values
dLmat = reshape(dLmat, length(thetagrid), length(epsgrid))
XXL = range(minimum(epsgrid), maximum(epsgrid), length=1000)
YYL = range(minimum(thetagrid), maximum(thetagrid), length=1000)
ZZL = griddata(epsgrid, thetagrid, dLmat, XXL, YYL)

drwmat = reshape(drwmat, length(thetagrid), length(epsgrid))
XXrw = range(minimum(epsgrid), maximum(epsgrid), length=1000)
YYrw = range(minimum(thetagrid), maximum(thetagrid), length=1000)
ZZrw = griddata(epsgrid, thetagrid, drwmat, XXrw, YYrw)

dPmat = reshape(dPmat, length(thetagrid), length(epsgrid))
XXP = range(minimum(epsgrid), maximum(epsgrid), length=1000)
YYP = range(minimum(thetagrid), maximum(thetagrid), length=1000)
ZZP = griddata(epsgrid, thetagrid, dPmat, XXP, YYP)

dwmat = reshape(dwmat, length(thetagrid), length(epsgrid))
XXw = range(minimum(epsgrid), maximum(epsgrid), length=1000)
YYw = range(minimum(thetagrid), maximum(thetagrid), length=1000)
ZZw = griddata(epsgrid, thetagrid, dwmat, XXw, YYw)

drmat = reshape(drmat, length(thetagrid), length(epsgrid))
XXr = range(minimum(epsgrid), maximum(epsgrid), length=1000)
YYr = range(minimum(thetagrid), maximum(thetagrid), length=1000)
ZZr = griddata(epsgrid, thetagrid, drmat, XXr, YYr)

dacrmat = reshape(dacrmat, length(thetagrid), length(epsgrid))
XXacrb = range(minimum(epsgrid), maximum(epsgrid), length=1000)
YYacrb = range(minimum(thetagrid), maximum(thetagrid), length=1000)
ZZacrb = griddata(epsgrid, thetagrid, dacrmat, XXacrb, YYacrb)

dmobmat = reshape(dmobmat, length(thetagrid), length(epsgrid))
XXmobb = range(minimum(epsgrid), maximum(epsgrid), length=1000)
YYmobb = range(minimum(thetagrid), maximum(thetagrid), length=1000)
ZZmobb = griddata(epsgrid, thetagrid, dmobmat, XXmobb, YYmobb)

# Multi-panel figure
p1 = contourf(YYL, XXL, ZZL, levels=10, xlabel="Theta", ylabel="Epsilon", title="Panel A : Population Treatment", fontsize=8)
p2 = contourf(YYw, XXw, ZZw, levels=10, xlabel="Theta", ylabel="Epsilon", title="Panel B : Wage Treatment", fontsize=8)
p3 = contourf(YYP, XXP, ZZP, levels=10, xlabel="Theta", ylabel="Epsilon", title="Panel C : Price Treatment", fontsize=8)
p4 = contourf(YYr, XXr, ZZr, levels=10, xlabel="Theta", ylabel="Epsilon", title="Panel D : Land Rent Treatment", fontsize=8)
p5 = contourf(YYrw, XXrw, ZZrw, levels=10, xlabel="Theta", ylabel="Epsilon", title="Panel E : Real Wage Treatment", fontsize=8)
p6 = contourf(YYacrb, XXacrb, ZZacrb, levels=10, xlabel="Theta", ylabel="Epsilon", title="Panel F : Incorrect Immobile Welfare", fontsize=8)

plot(p1, p2, p3, p4, p5, p6, layout=(3, 2), size=(800, 800), legend=false)
