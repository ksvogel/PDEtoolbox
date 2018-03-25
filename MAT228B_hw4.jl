#=
Kaela S. Vogel
2018-03-25
Homeworks for MAT 228

=#
include("PDEtool.jl")
using PDEtool
using Plots
plotlyjs()


##############################################################

#= Homework 4, problem 1
Solve the advection equation withperiodic boundry conditions using Lax-Wendroff and upwinding.

# Need to use interpolation for the final timestep since linspace does not give the final time exactly at 1
=#

a =  1.0
h = 1/2^7 # Grid spacing
k = 0.9*h/a # Tine stepping
s = Bool(1) # Switch variable to give basic mesh
funcRHS( x, t) = 0
FDM = 1 # Lax-Wendroff
#FDM = 2 # Upwinding
#u_x0(x)= 0.5*cos.(2*pi*x) + 0.5
u_x0(x)= 0.5*cos.(2*pi*x) + 0.5
u_x0(x) = abs(x - 0.5) < 0.25 ? 1 : 0 # initial condtion at t = 0
v, error = advectionFDM(h, k, funcRHS, u_x0, a, s, FDM)

# Plotting
Plots.plotlyjs()
hmesh, kmesh, F, M, T = meshmaker(funcRHS, h, k, s)
l = size(v[:,1])
x = hmesh
y_vals = [v[:,1] v[:,50]-ones(l) v[:,100]-2*ones(l) v[:, end]-3*ones(l)]
labels = [string("Time = ", kmesh[1]) string("Time = ", kmesh[50])  string("Time = ", kmesh[100])  string("Time = ", kmesh[end])  ]
Plots.plot(x, y_vals, linewidth=2, alpha=0.6, labels = labels)

###############################################################################
# Grid refinement study Lax-Wendroof/Uwinding continuous/discontinuous I.C. table 1
# Norm Convergence table 2
a =  1.0
s = Bool(1) # Switch variable to give basic mesh
funcRHS( x, t) = 0 # Dummy, these problems don't have a forcing term
u_x01(x)= 0.5*cos.(2*pi*x) + 0.5
u_x02(x) = abs(x - 0.5) < 0.25 ? 1 : 0 # initial condtion at t = 0
hstart =  3
hend = 11
spacing = [1/2^j for j in hstart:1:hend]
vsLW1 = Vector{Any}(length(spacing))
vsUP1 = Vector{Any}(length(spacing))
vsLW2 = Vector{Any}(length(spacing))
vsUP2 = Vector{Any}(length(spacing))
normerrors1 = zeros(6,length(spacing))
normerrors2 = zeros(6,length(spacing))

for i in 1:length(spacing)
    # Continuous Initial data
    vsLW1[i], normerrors1[1:3, i]  = advectionFDM(spacing[i], 0.9*spacing[i]/a, funcRHS, u_x01, a, s, 1)
    vsUP1[i], normerrors1[4:6, i]  = advectionFDM(spacing[i], 0.9*spacing[i]/a, funcRHS, u_x01, a, s, 2)
    # Discontinuous Initial data
    vsLW2[i], normerrors2[1:3, i]  = advectionFDM(spacing[i], 0.9*spacing[i]/a, funcRHS, u_x02, a, s, 1)
    vsUP2[i], normerrors2[4:6, i]  = advectionFDM(spacing[i], 0.9*spacing[i]/a, funcRHS, u_x02, a, s, 2)
end

refinement_ratios1 = (normerrors1[:,1:(end-1)] ./ normerrors1[:, 2:end])'
refinement_ratios2 = (normerrors1[:,1:(end-1)] ./ normerrors1[:, 2:end])'



###############################################################################
# Solve the advection equation with crank nicoslson

a =  1.0
h = 1/2^7 # Grid spacing
k = 0.9*h/a # Tine stepping
s = Bool(1) # Switch variable to give basic mesh
funcRHS( x, t) = 0
FDM = 3 # CN
#u_x0(x)= .5*cos.(2*pi*x) + 0.5
u_x0(x) = abs(x - 0.5) < 0.25 ? 1 : 0 # initial condtion at t = 0
v, error = advectionFDM(h, k, funcRHS, u_x0, a, s, FDM)

# Plotting
heatmap(v)
hmesh, kmesh, F, M, T = meshmaker(funcRHS, h, k, s)
l = size(v[:,1])
Plots.plotlyjs()
x = hmesh
y_vals = [v[:,1] v[:,50]-ones(l) v[:,100]-2*ones(l) v[:, end]-3*ones(l)]
labels = [string("Time = ", kmesh[1]) string("Time = ", kmesh[50])  string("Time = ", kmesh[100])  string("Time = ", kmesh[end])  ]
Plots.plot(x, y_vals, linewidth=2, alpha=0.6, labels = labels)

# Plot the relative phase error

x = [j for j in 0:.05: pi]
y = -atan.(-0.9*sin.(x)./(1-0.9^2*sin.(x).^2./4))./(0.9.*x)
Plots.plot(x, y, linewidth=2, alpha=0.6)

# In part A we found using Von Neumann analysis that we expect no amplitude error. The relative phase
