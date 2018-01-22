#=
Kaela S. Vogel
2018-01-01
Homeworks for MAT 228

=#

##############################################################

#= Homework 1, problem 2
This is the setup and the function call to solve the heat equation with Crank-Nicolson and then perform a refinement study to show it is second order accurate in space and time.

=#
include("PDEtool.jl")
using PDEtool
using Plots
plotlyjs()

cd("/home/kaela/Documents/GithubRepositories/PDE_solvers")



h = 1/2^4
k = 1/2^4
funcRHS( x, t) = 0 # 1 - exp.( - t)
u_0t(t) = 0 # boundary condition at x = 0
u_1t(t) = 0 # boundary condition at x = 1
#u_x0(x)= 0
u_x0(x) = exp.(-(x-.5).^2)-exp.(-.5^2)
#u_x0(x) = x < 0.5 ? x : 1-x # initial condtion at t = 0
a =  1.0
v = PDEtool.CN_heat(h, k, funcRHS,  u_0t, u_1t, u_x0, a)
heatmap(v)

# The revision  study
hstart =  3
hend = 10
spacing = [1/2^j for j in hstart:1:hend]
v_unrefined = Vector{Any}(length(spacing)) # This will be an array of matrices of different dims

for j in 1:length(spacing) # The indexing assumes h and k are equal
    v = PDEtool.CN_heat(spacing[j], spacing[j], funcRHS,  u_0t, u_1t, u_x0, a)
    v_unrefined[j] = v
end

ratios, v_restrictedsCNH = refinement_ratios(v_unrefined)


##################################################################

h = 1/2^1 # Space length
k = 1/2^1 # Time length
funcRHS( x, t) = 0 #1 - exp.( - t)
u_0t(t) = 0 # boundary condition at x = 0
u_1t(t) = 0 # boundary condition at x = 1
#u_x0(x)= 0
u_x0(x) = exp.(-(x-.5).^2)-exp.(-.5^2)
#u_x0(x) = x < 0.5 ? 1 : 0 # initial condtion at t = 0
#u_x0(x) = x < 0.5 ? x : 1-x # initial condtion at t = 0
#a =  0.1
#u_x0(x) = x < 0.5 ? 1 : 0 # initial condtion at t = 0
a =  1.0
v = PDEtool.BDF2_heat(h, k, funcRHS,  u_0t, u_1t, u_x0, a)
heatmap(v)

# The revision  study
hstart =  3
hend = 10
spacing = [1/2^j for j in hstart:1:hend]
v_unrefined = Vector{Any}(length(spacing)) # This will be an array of matrices of different dims

for j in 1:length(spacing) # The indexing assumes h and k are equal
    v = PDEtool.BDF2_heat(spacing[j], spacing[j], funcRHS,  u_0t, u_1t, u_x0, a)
    v_unrefined[j] = v
end

ratios, v_restricteds = refinement_ratios(v_unrefined)
v_unrefined[1]
v_unrefined[2]
#################################################################

#HW1 Problem 3


h = 0.02
k = 0.1
funcRHS( x, t) = 0 #
u_0t(t) = 1 # boundary condition at x = 0
u_1t(t) = 0 # boundary condition at x = 1
u_x0(x) = x < 0.5 ? 1 : 0 # initial condtion at t = 0
a =  1.0
#v_1  = PDEtool.CN_heat(h, k, funcRHS,  u_0t, u_1t, u_x0, a)
m = PDEtool.CN_heat(h/10, k/10, funcRHS,  u_0t, u_1t, u_x0, a)
heatmap(m)

# This code block runs Crank Nicolson for a selection of grid spacings and then makes a heat map for each run combined in one figure.
v_1 = PDEtool.CN_heat(h, k, funcRHS,  u_0t, u_1t, u_x0, a)
v_2 = PDEtool.CN_heat(h, k/10, funcRHS,  u_0t, u_1t, u_x0, a)
v_3 = PDEtool.CN_heat(h, k/100, funcRHS,  u_0t, u_1t, u_x0, a)
v_4 = PDEtool.CN_heat(h/10, k, funcRHS,  u_0t, u_1t, u_x0, a)

hmesh, kmesh, F, M, T = meshmaker(funcRHS, h, k)

######### Sub plot (a)
anames1 = ["" for x=kmesh]
anames2 = ["" for x=hmesh]
anames1[1:2:T] = ["$(x)" for x in kmesh[1:2:T]]
anames2[1:8:M] = ["$(x)" for x in hmesh[M:-8:1]]
fig1_11 = heatmap(linspace(0,1,size(v_1,2)),1:-1/size(v_1,1):0,v_1,xlabel="foo",ylabel="bar",xticks= (kmesh, anames1), yticks= (hmesh, bnames2))

######### Sub plot (b)
fig1_12 = heatmap(v_2)

######### Sub plot (c)
fig1_21 = heatmap(v_3)

######### Sub plot (d)
fig1_22 = heatmap(v_4)

plot(fig1_11, fig1_12, fig1_21, fig1_22)
