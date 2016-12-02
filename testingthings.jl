# testing ideas script

B = ones(10,10)
M = 10
rb = 2:M-1
for j in 2:M-1
red = filter(k -> iseven(k+j), 2:M-1)
println(red)
end

#include("PDEtool.jl")
using PDEtool
h = 2.0^(-3)

# Set up  initial guess and boundary data, this is homogenous Dirichlet
u = zeros(size(F))

####### Set up our R.H.S. function
funcRHS = (x,y) -> -exp(-(x - 0.25)^2 - (y - 0.6)^2)
mesh, F = PDEtool.h_space(funcRHS, h)

v = PDEtool.squish(F)
u = PDEtool.foomp(v)
