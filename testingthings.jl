# testing ideas script

B = ones(10,10)
M = 10
rb = 2:M-1
red = 1
for j in 2:M-1
red = filter(k -> iseven(k+j), 2:M-1)
#println(red)
end

#include("PDEtool.jl")
using PDEtool
h = 2.0^(-3)

####### Set up our R.H.S. function
funcRHS = (x,y) -> -exp(-(x - 0.25)^2 - (y - 0.6)^2)
mesh, F = PDEtool.h_space(funcRHS, h)
# Set up  initial guess and boundary data, this is homogenous Dirichlet
u = zeros(size(F))
println(u)
turtles = 3
s1 = 3
s2 = 3
tol = 10.0^(-8)

uSOR, iters = SOR(h, 1000, tol, F)
Ua, residual = MG_vcycle(h, F, uSOR, turtles, s1, s2, tol)
uout, res, maxiter = gauss_sidel(u, h, 1000, tol, 1, F)
iterdiff = vecnorm(uout-Ua[turtles],1)

uSOR, iters = SOR(h, 1000, tol, F)
iterdiff = vecnorm(uSOR-Ua[turtles],1)
iterdiff = vecnorm(uSOR-uout,1)

c = foomp(zeros(size(F)))
#=
h = 2.0^(-2)
funcRHS = (x,y) -> -exp(-(x - 0.25)^2 - (y - 0.6)^2)
mesh, F = PDEtool.h_space(funcRHS, h)
u = zeros(size(F))
u, res, maxiter = gauss_sidel(u, h, 50, .000001, 1, F)
=#
