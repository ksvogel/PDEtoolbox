#=
Kaela S. Vogel
2018-03-25
Homeworks for MAT 228

=#
include("PDEtool.jl")
using PDEtool
using PyPlot
using LatexPrint



##############################################################

#= Homework 4, problem 1
Solve the advection equation withperiodic boundry conditions using Lax-Wendroff and upwinding.
=#
# For testing
a =  1.0
h = 1/(2^5) # Grid spacin
k = 8*h/10 # Tine stepping
s = Bool(1) # Switch variable to give basic mesh
funcRHS( x, t) = 0
FDM = 1 # Lax-Wendroff
#FDM = 2 # Upwinding
#u_x0(x)= 0.5*cos.(2*pi*x) + 0.5
u_x0(x)= 0.5*cos.(2*pi*x) + 0.5
#u_x0(x) = abs(x - 0.5) < 0.25 ? 1 : 0 # initial condtion at t = 0
#v, error = PDEtool.advectionFDM(h, k, funcRHS, u_x0, a, s, FDM)


###############################################################################
# Grid refinement study Lax-Wendroof/Uwinding continuous/discontinuous I.C. table 1
# Norm Convergence table 2

a =  1.0
c = 0.9
s = Bool(1) # Switch variable to give basic mesh
u_x01(x)= 0.5*cos.(2*pi*x) + 0.5
u_x02(x) = abs(x - 0.5) < 0.25 ? 1 : 0 # initial condtion at t = 0
hstart =  4
hend = 10
Nts = [10*2^j for j in hstart:1:hend]
Nxs = round.(Int, [c*(nt) for nt in Nts ])
ks = [1/(Nt) for Nt in Nts]
hs = [1/Nx for Nx in Nxs] # I will use linspace to index from 0:1-h since we loop back around on the mesh for periodic conditions
vsLW1 = Vector{Any}(length(Nts))
vsUP1 = Vector{Any}(length(Nts))
vsLW2 = Vector{Any}(length(Nts))
vsUP2 = Vector{Any}(length(Nts))
normerrors1 = zeros(6,length(Nts))
normerrors2 = zeros(6,length(Nts))

for i in 1:length(Nts)
    # Continuous Initial data
    vsLW1[i], normerrors1[1:3, i]  = PDEtool.advectionFDM(hs[i], ks[i], Nxs[i], Nts[i], u_x01, a, s, 1)
    vsUP1[i], normerrors1[4:6, i]  = PDEtool.advectionFDM(hs[i], ks[i], Nxs[i], Nts[i], u_x01, a, s, 2)

    # Discontinuous Initial data
    vsLW2[i], normerrors2[1:3, i]  = PDEtool.advectionFDM(hs[i], ks[i], Nxs[i], Nts[i], u_x02, a, s, 1)
    vsUP2[i], normerrors2[4:6, i]  = PDEtool.advectionFDM(hs[i], ks[i], Nxs[i], Nts[i], u_x02, a, s, 2)

end

refinement_ratios1 = (normerrors1[:,1:(end-1)] ./ normerrors1[:, 2:end])'
refinement_ratios2 = (normerrors2[:,1:(end-1)] ./ normerrors2[:, 2:end])'
R = round.(refinement_ratios1; refinement_ratios2],5, 10)
R = [ [Nxs[2:end] ; Nxs[2:end] ] R ]
lap(R)
NE = signif.([normerrors1'; normerrors2'],5)
lap(NE)




###############################################################################
# Solve the advection equation with crank nicoslson

a =  1.0
c = 0.9
s = Bool(1) # Switch variable to give basic mesh
Nt = 10*2^6
Nx = Int64(c*(Nt))
k = 1/Nt
h = 1/Nx
FDM = 3 # CN
#u_x0(x)= .5*cos.(2*pi*x) + 0.5
u_x0(x) = abs(x - 0.5) < 0.25 ? 1 : 0 # initial condtion at t = 0
v, error = PDEtool.advectionFDM(h, k, Nx, Nt, u_x0, a, s, FDM)

# Plotting
fig = figure("CN_advection",figsize=(10,10))
x = collect(linspace(0, 1-h, Nx))
l = size(v[:,1])
snaps = [ 1 10 20 30 40 341 441 541 640]
subplot(339) # Create the third plot of a 4x4 group of subplots
suptitle("Crank-Nicolson for Advection, discontinuous initial conditions") # Supe title, title for all subplots combined
for sn in 1:1:9
    sp = parse(Int64, string("33", sn))
    subplot(sp)
    timestr = string("T = ", k*snaps[sn])
    title(timestr)
    ax = gca()
    setp(ax[:get_xticklabels](),visible=false) # Disable x tick labels
    setp(ax[:get_yticklabels](),visible=false) # Disable y tick labels
    y = v[:,snaps[sn]]
    #labels = [string("Time = ", kmesh[1]) string("Time = ", kmesh[50])     string("Time = ", kmesh[100])  string("Time = ", kmesh[end])  ]
    plot(x, y, linewidth=1, alpha=0.6)
    fig[:canvas][:draw]()
end # Update the figure, required when doing additional modifications


##############################################################################
# Plot the relative phase error
x = [j for j in 0:.05: pi]
y = -atan.(-0.9*sin.(x)./(1-0.9^2*sin.(x).^2./4))./(0.9.*x)
plot(x, y, linewidth=2, alpha=0.6)
title("Relative Phase Error for advection discretized with Crank-Nicolson")
ylabel(L"$\frac{\arg(g(\theta))}{-\nu \theta}$")
xlabel(L"$\theta$")
#savefig("phase_compare.png", type="png", dpi=300)

# In part A we found using Von Neumann analysis that we expect no amplitude error. The relative phase
