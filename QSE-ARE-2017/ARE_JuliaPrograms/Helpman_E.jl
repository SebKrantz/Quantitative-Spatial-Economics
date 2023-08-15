# Monte Carlo for Annual Review
# Quantitative Spatial Model
# Helpman (1998) model
# Countries and regions

# July, 2016

# *********************
# **** Load Files  **** 
# *********************

function load_dir(dir::String)
    files = readdir(dir)
    for file in files
        if endswith(file, ".jl")
            include(joinpath(dir, file))
        end
    end
end

cd("literature/Redding QSE/ARE_JuliaPrograms/")
load_dir("functions")


using Distributions
using LinearAlgebra
using Statistics
using StatsBase

# No need to declare globals in Julia (only if they are mutable)
# global alpha 
# global sigma 
# global LL 
# global LLwest 
# global LLeast 
# global F

# ************************
# **** Initialization ****
# ************************

# Set default random number stream
# s = MersenneTwister(1)

# *************************
# **** Distance matrix ****
# *************************

# Locations on a N * N latitude and longitude grid
# Implies N * N locations
# Distance matrix is N * N

N = 30
NN = N*N

# Other latitude-longitude grid
ltd = range(0,stop=4,length=N)' 
lgd = range(0,stop=4,length=N)

# Transport costs for each point on grid
tt = 1
τ = fill(tt, (N,N))

# Compute weighted distance using these transport costs

dist = Matrix{Float64}(undef, NN, NN)


# # Convert Matrix to nested tuples
# function mattup(A)
#     ntuple(i -> ntuple(j -> A[i,j], size(A,2)), size(A,1))
# end

# # Function to convert vec(A) to a tuple
# function vectotuple(v::Vector)
#     t = ntuple(i -> v[i], length(v))
#     return t
# end

#  using ImageMorphology

## Same done manually
function my_dist_transform(X)
    N = size(X,1)
    dist = Matrix{Float64}(undef, N, N)
    x, y = Tuple(findall(X .== 1)[1])
    for i in 1:N, j in 1:N
        dist[i,j] = sqrt((i-x)^2 + (j-y)^2)
    end
    return dist
end

# Distance weighting does not work yet
for z in 1:NN
    seed = falses(N,N)
    seed[z] = true
    # w = vectotuple(vec(τ))
    temp = my_dist_transform(seed) # distance_transform(feature_transform(seed)) # w # graydist(τ,seed,'quasi-euclidean')
    dist[z,:] = temp[:]
end

# Own iceberg transport costs are one
dist[diagind(dist)] .= 1 

# Trade costs are a power function of effective distance
dist = dist.^0.33

# Define east and west as two countries
Iwest = falses(N,N) 
Ieast = falses(N,N)
Iwest[:,1:Int(N/2)] .= 1
Ieast[:,Int(N/2+1):N] .= 1
Iwest = vec(Iwest)
Ieast = vec(Ieast) 

# Border friction between grid points
bord = fill(2.0, (NN,NN))  
bord[diagind(bord)] .= 1

# Border friction between countries
bordcty = ones(NN,NN)
bordcty[Iwest.==1, Ieast.==1] .= 2
bordcty[Ieast.==1, Iwest.==1] .= 2

# Counterfactual border friction between grid points
cbord = ones(NN,NN)

# Counterfactual border friction between countries  
cbordcty = ones(NN,NN)

# **************************
# **** Parameterization ****
# **************************

# Share of goods in consumption expenditure (1-housing share)
alpha = 0.75

# Elasticity of substitution
sigma = 5

# ************************************ 
# **** Random productivity shocks ****
# ************************************

# a = rand(Normal(0,1), NN, 1)
# Read a.csv
using CSV
using DataFrames
a = CSV.read("a.csv", DataFrame; header = false) # Read a.csv
a = exp.(Matrix(a))
a[Iwest.==1] .= a[Iwest.==1] ./ geomean(a[Iwest.==1])
a[Ieast.==1] .= a[Ieast.==1] ./ geomean(a[Ieast.==1]) 

println("Summary statistics productivities")
println("mean(a) std(a) max(a) min(a)")
[mean(a) std(a) maximum(a) minimum(a)]

# **************************
# **** Other Parameters ****  
# **************************:

# Land area
H = 100*ones(NN,1)   

# Aggregate labor Supply 
LL = 153889 # US civilian labor force 2010 (Statistical Abstract, millions)
LLwest = (sum(Iwest)/(sum(Iwest)+sum(Ieast)))*LL 
LLeast = (sum(Ieast)/(sum(Iwest)+sum(Ieast)))*LL

# Fixed production cost
F = 1

# ********************************
# **** Matrix of Fundamentals ****
# ********************************

fund = zeros(NN,4)
fund[:,1] = a
fund[:,2] = H
fund[:,3] = Iwest
fund[:,4] = Ieast

# ****************************************************************************
# **** Open Economy Solve for Endogenous Variables in Initial Equilibrium ****
# ****************************************************************************

println(">>>> Start Wage and Population Convergence <<<<")
w,L,tradesh,dtradesh,converge,xtic = solveHLwCtyOpen_E(fund,dist,bord,bordcty,NN)
println(">>>> Wage and Population System Converged <<<<")
println(">>>> Check Wage and Population Convergence (Should ==1) <<<<")
println(converge) 
println(">>>> Elapsed Time in Seconds <<<<")
println(xtic)

# Price index
P = Hpindex(fund,L,w,dtradesh)   

# Land prices
r = Hlandprice(fund,L,w)

# Real wage  
realwage = Hrealw(fund,L,w,tradesh)

# ************************************************************************
# ***** Counterfactual Eliminating Border Frictions Between Countries ****
# ************************************************************************

println(">>>> Start Wage and Population Convergence <<<<")
cw,cL,ctradesh,cdtradesh,cconverge,xtic = solveHLwCtyOpen_E(fund,dist,bord,cbordcty,NN)
println(">>>> Wage and Population System Converged <<<<")
println(">>>> Check Wage and Population Convergence (Should ==1) <<<<")
println(cconverge)
println(">>>> Elapsed Time in Seconds <<<<") 
println(xtic)

# Price index
cP = Hpindex(fund,cL,cw,cdtradesh)

# Land prices
cr = Hlandprice(fund,cL,cw) 

# Real wage
crealwage = Hrealw(fund,cL,cw,ctradesh)  

# Welfare gains
welfgain = Hwelfaregains(ctradesh,tradesh,cL,L)
println(">>>> Welfare Gains <<<<")
welfgain = round.(welfgain, digits=4) 
unique(welfgain)

# **************************************************************************
# ***** Counterfactual Eliminating Border Frictions Between Grid Points ****
# **************************************************************************

println(">>>> Start Wage and Population Convergence <<<<")
ccw,ccL,cctradesh,ccdtradesh,ccconverge,xtic = solveHLwCtyOpen_E(fund,dist,cbord,bordcty,NN)
println(">>>> Wage and Population System Converged <<<<")  
println(">>>> Check Wage and Population Convergence (Should ==1) <<<<")
println(ccconverge)
println(">>>> Elapsed Time in Seconds <<<<")
println(xtic)

# Price index
ccP = Hpindex(fund,ccL,ccw,ccdtradesh)

# Land prices 
ccr = Hlandprice(fund,ccL,ccw)

# Real wage
ccrealwage = Hrealw(fund,ccL,ccw,cctradesh)

# Welfare gains
welfgain = Hwelfaregains(cctradesh,tradesh,ccL,L)
println(">>>> Welfare Gains <<<<")
welfgain = round.(welfgain, digits=4)
unique(welfgain)

# ***********************************************
# **** Three-Dimensional Initial Equilibrium ****
# ***********************************************

# LOG PRODUCTIVITY  
amat = reshape(log.(a),N,N) 

# LOG POPULATION
Lmat = reshape(log.(L),N,N)   

# LOG WAGE
wmat = reshape(log.(w),N,N)

# LOG RELATIVE LAND PRICE
rmat = reshape(log.(r),N,N)  

# PRICE INDEX
Pmat = reshape(log.(P),N,N)   

using Plots

# MULTI-PANEL FIGURE

# Productivity
heatmap(amat, xlabel="Longitude", ylabel="Latitude", 
        title="Log Productivity", fontsize=12, c =:viridis) 

# MULTI-PANEL FIGURE 
p = plot(layout = (2,2))

# Population
heatmap!(p[1], Lmat, xlabel="Longitude", ylabel="Latitude",  
    title="Panel A: Log Population", fontsize=8, c =:viridis)

# Wages
heatmap!(p[2], wmat, xlabel="Longitude", ylabel="Latitude",
    title="Panel B: Log Wages", fontsize=8, c =:viridis) 

# Land prices  
heatmap!(p[3], rmat, xlabel="Longitude", ylabel="Latitude",
    title="Panel C: Log Land Prices", fontsize=8, c =:viridis)

# Price Index
heatmap!(p[4], Pmat, xlabel="Longitude", ylabel="Latitude", 
    title="Panel D: Log Price Index", fontsize=8, c =:viridis)

savefig("initial.pdf")

# **************************************************************************
# **** Three-Dimensional Eliminating Border Frictions Between Countries ****
# **************************************************************************

dL = cL./L
ldL = log.(dL)
dw = cw./w  
ldw = log.(dw)
dr = cr./r
ldr = log.(dr) 
dP = cP./P
ldP = log.(dP)

# POPULATION
dLmat = reshape(ldL,N,N)  

# PRICE INDEX
dPmat = reshape(ldP,N,N)

# WAGE
dwmat = reshape(ldw,N,N)

# RELATIVE LAND PRICE
drmat = reshape(ldr,N,N)

# MULTI-PANEL FIGURE DIFFERENCES
p = plot(layout = (2,2))

# Population
heatmap!(p[1], dLmat, xlabel="Longitude", ylabel="Latitude",  
  title="Panel A: Log Relative Population", fontsize=8, c =:viridis)

# Wage
heatmap!(p[2], dwmat, xlabel="Longitude", ylabel="Latitude",
   title="Panel B: Log Relative Wages", fontsize=8, c =:viridis)

# Land price
heatmap!(p[3], drmat, xlabel="Longitude", ylabel="Latitude", 
   title="Panel C: Log Relative Land Rents", fontsize=8, c =:viridis)
   
# Price index
heatmap!(p[4], dPmat, xlabel="Longitude", ylabel="Latitude",
   title="Panel D: Log Relative Price Index", fontsize=8, c =:viridis) 

savefig("c.pdf")

# MULTI-PANEL FIGURE LEVELS
p = plot(layout = (2,2))

# Population
heatmap!(p[1], Lmat + dLmat, xlabel="Longitude", ylabel="Latitude",
  title="Panel A: Log Population", fontsize=8, c =:viridis)

# Wage
heatmap!(p[2], wmat + dwmat, xlabel="Longitude", ylabel="Latitude", 
 title="Panel B: Log Wages", fontsize=8, c =:viridis)

# Land price
heatmap!(p[3], rmat + drmat, xlabel="Longitude", ylabel="Latitude",
 title="Panel C: Log Land Rents", fontsize=8, c =:viridis)

# Price index 
heatmap!(p[4], Pmat + dPmat, xlabel="Longitude", ylabel="Latitude", 
 title="Panel D: Log Price Index", fontsize=8, c =:viridis)

savefig("c_lev.pdf")

# ****************************************************************************
# **** Three-Dimensional Eliminating Border Frictions Between Grid Points **** 
# ****************************************************************************

ddL = ccL./L
lddL = log.(ddL)
ddw = ccw./w
lddw = log.(ddw)  
ddr = ccr./r
lddr = log.(ddr)
ddP = ccP./P
lddP = log.(ddP)

# POPULATION
ddLmat = reshape(lddL,N,N)

# PRICE INDEX 
ddPmat = reshape(lddP,N,N)

# WAGE
ddwmat = reshape(lddw,N,N)

# RELATIVE LAND PRICE
ddrmat = reshape(lddr,N,N)

# MULTI-PANEL FIGURE DIFFERENCES 
p = plot(layout = (2,2))

# Population
heatmap!(p[1], ddLmat, xlabel="Longitude", ylabel="Latitude",
  title="Panel A: Log Relative Population", fontsize=8, c =:viridis) 

# Wage
heatmap!(p[2], ddwmat, clim=(-0.05, 0.05), xlabel="Longitude", ylabel="Latitude",
 title="Panel B: Log Relative Wages (Truncated)", fontsize=8, c =:viridis)

# Land price  
heatmap!(p[3], ddrmat, xlabel="Longitude", ylabel="Latitude",
 title="Panel C: Log Relative Land Rents", fontsize=8, c =:viridis)

# Price index
heatmap!(p[4], ddPmat, clim=(-1.35, -1.15), xlabel="Longitude", ylabel="Latitude", 
 title="Panel D: Log Relative Price Index (Truncated)", fontsize=8, c =:viridis)

savefig("cc.pdf")

# MULTI-PANEL FIGURE LEVELS
p = plot(layout = (2,2))

# Population
heatmap!(p[1], Lmat + ddLmat, xlabel="Longitude", ylabel="Latitude", 
 title="Panel A: Log Population", fontsize=8, c =:viridis)

# Wage  
heatmap!(p[2], wmat + ddwmat, xlabel="Longitude", ylabel="Latitude",
 title="Panel B: Log Wages", fontsize=8, c =:viridis)

# Land price
heatmap!(p[3], rmat + ddrmat, xlabel="Longitude", ylabel="Latitude",
 title="Panel C: Log Land Rents", fontsize=8, c =:viridis)

# Price index
heatmap!(p[4], Pmat + ddPmat, xlabel="Longitude", ylabel="Latitude",
 title="Panel D: Log Price Index", fontsize=8, c =:viridis)

savefig("cc_lev.pdf")