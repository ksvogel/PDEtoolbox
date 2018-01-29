#=
Kaela S. Vogel
2018-01-01
Homeworks for MAT 228

=#

##############################################################

#= Homework 2, problem 2
This is the setup and the function call to solve the heat equation

u_t = \Delta u on  \Omega = (0,1) X (0,1
u(0,0,t) = 0  boundary
u(x,y,0) = exp[-100((x - 0.3)^2 + (y - 0.4)^2)]

with Crank-Nicolson 2D and Forward Euler solvers

=#
include("PDEtool.jl")
using PDEtool
using Plots
plotlyjs()

#cd("/home/kaela/Documents/GithubRepositories/PDE_solvers")

h = 1/2^4 # Grid spacing
k = 1/2^4 # Tine stepping
funcRHS( x, y, t) = 0
u_0t(t) = 0 # boundary condition at x = 0
u_1t(t) = 0 # boundary condition at x = 1
u_xy0(x, y)= exp(-100((x - 0.3)^2 + (y - 0.4)^2))
#_x0(x) = exp.(-(x-.5).^2)-exp.(-.5^2)
#u_x0(x) = x < 0.5 ? x : 1-x # initial condtion at t = 0
a =  1.0
v = PDEtool.CN_heat2D(h, k, funcRHS,  u_0t, u_1t, u_xy0, a)

##############################################################

#= Homework 1, problem 2
This is the setup and the function call to solve the heat equation

u_t = 0.01u_{xx} + 1 - \exp(-t), 0 < x < 1
u(0,t) = 0,
u(1,t) = 0
u(x,0) = 0

with Crank-Nicolson and then perform a refinement study to show it is second order accurate in space and time.

=#

h = 1/2^4 # Grid spacing
k = 1/2^4 # Tine stepping
funcRHS( x, t) = 1 - exp.( - t)
u_0t(t) = 0 # boundary condition at x = 0
u_1t(t) = 0 # boundary condition at x = 1
u_x0(x)= 0
#_x0(x) = exp.(-(x-.5).^2)-exp.(-.5^2)
#u_x0(x) = x < 0.5 ? x : 1-x # initial condtion at t = 0
a =  0.1
v = PDEtool.CN_heat(h, k, funcRHS,  u_0t, u_1t, u_x0, a)
heatmap(v)

# The grid refinement study using the infinity-norm
hstart =  3
hend = 10
spacing = [1/2^j for j in hstart:1:hend]
v_unrefined = Vector{Any}(length(spacing)) # This will be an array of matrices of different dims

for j in 1:length(spacing) # The indexing assumes h and k are equal
    v = PDEtool.CN_heat(spacing[j], spacing[j], funcRHS,  u_0t, u_1t, u_x0, a)
    v_unrefined[j] = v
end

ratios, v_restrictedsCNH = refinement_ratios(v_unrefined, spacing, spacing)
writecsv("hw_1_refinement.txt", ratios)


##### Creates a heatmap of the Crank Nicolson solution for sucessively refined grids
######### Sub plot (a)
fig1_11 = heatmap(linspace(0,1,size(v_unrefined[1],2)),1:-1/size(v_unrefined[1],1):0,v_unrefined[1], title="(a) h, k = 1/2^3")
######### Sub plot (b)
fig1_12 = heatmap(linspace(0,1,size(v_unrefined[2],2)),1:-1/size(v_unrefined[2],1):0,v_unrefined[2], title="(b) h, k = 1/2^4")
######### Sub plot (c)
fig1_21 = heatmap(linspace(0,1,size(v_unrefined[3],2)),1:-1/size(v_unrefined[3],1):0,v_unrefined[3],xlabel="t",ylabel="x", title="(c) h, k = 1/2^5")
######### Sub plot (d)
fig1_22 = heatmap(linspace(0,1,size(v_unrefined[7],2)),1:-1/size(v_unrefined[7],1):0,v_unrefined[7], title="(d) h, k = 1/2^10")

plot(fig1_11, fig1_12, fig1_21, fig1_22)

##############################################################

#= Homework 1, problem 3
This is the setup and the function call to solve the heat equation

u_t = 1u_{xx}, 0 < x < 1
u(0,t) = 1,
u(1,t) = 0
u(x,0) = 1 if x < 0.5, 0 if x >= 0.5

with Crank-Nicolson and BDF-2

=#

h = 0.02 # grid spacing
k = 0.1 # time stepping
funcRHS( x, t) = 0 #
u_0t(t) = 1 # boundary condition at x = 0
u_1t(t) = 0 # boundary condition at x = 1
u_x0(x) = x < 0.5 ? 1 : 0 # initial condtion at t = 0
a =  1.0
#v_1  = PDEtool.CN_heat(h, k, funcRHS,  u_0t, u_1t, u_x0, a)
m = PDEtool.BDF2_heat(h, k, funcRHS,  u_0t, u_1t, u_x0, a)
heatmap(linspace(0,1,size(m,2)),1:-1/size(m,1):0,m, xlabel="t",ylabel="x", title="h = 0.02, k = 0.1")

# This code block runs Crank Nicolson for a selection of grid spacings and then makes a heat map for each run combined in one figure.
v_1 = PDEtool.CN_heat(h, k, funcRHS,  u_0t, u_1t, u_x0, a)
v_2 = PDEtool.CN_heat(h, k/10, funcRHS,  u_0t, u_1t, u_x0, a)
v_3 = PDEtool.CN_heat(h, k/100, funcRHS,  u_0t, u_1t, u_x0, a)
v_4 = PDEtool.CN_heat(h/10, k, funcRHS,  u_0t, u_1t, u_x0, a)

######### Sub plot (a)
fig1_11 = heatmap(linspace(0,1,size(v_1,2)),1:-1/size(v_1,1):0,v_1, title="(a) h = 0.02, k = 0.1")
######### Sub plot (b)
fig1_12 = heatmap(linspace(0,1,size(v_2,2)),1:-1/size(v_2,1):0,v_2, title="(b) h = 0.02, k = 0.01")
######### Sub plot (c)
fig1_21 = heatmap(linspace(0,1,size(v_3,2)),1:-1/size(v_3,1):0,v_3,xlabel="t",ylabel="x", title="(c) h = 0.02, k = 0.001")
######### Sub plot (d)
fig1_22 = heatmap(linspace(0,1,size(v_4,2)),1:-1/size(v_4,1):0,v_4, title="(d) h = 0.002, k = 0.1")
gui()
plot(fig1_11, fig1_12, fig1_21, fig1_22)
#savefig(plot(fig1_11, fig1_12, fig1_21, fig1_22) , "hw1_fig2.png")
