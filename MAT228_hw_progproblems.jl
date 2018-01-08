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

cd("/home/kaela/Documents/GithubRepositories/PDE_solvers")

h = 1/2^4
k = 1/2^4
funcRHS( x, t) = 1 - exp.( - t)
u_0t(t) = 0 # boundary condition at x = 0
u_1t(t) = 0 # boundary condition at x = 1
u_x0(x) = 0 # initial condtion at t = 0
a = 0.01
v = PDEtool.CN_heat(h, k, funcRHS,  u_0t, u_1t, u_x0, a)
