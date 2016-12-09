# testing ideas script
using PDEtool

turtles = 3
h = 2.0^(-turtles)
funcRHS = (x,y) -> -exp(-(x - 0.25)^2 - (y - 0.6)^2)
#funcRHS = (x,y) -> 2*(x^ + y^2 - x - y)
mesh, F = PDEtool.h_space(funcRHS, h)
# Set up  initial guess and boundary data, this is homogenous Dirichlet
u0 = zeros(size(F))
tol = 10.0^(-10)
v = [1 3]
precon = "np"
maxiter = 100
u, r, k = PDEtool.PCG(u0, F, h, tol, precon)
uturtle, residual, maxiter = PDEtool.multigrid(u0, F, maxiter, tol, turtles, v)

M = size(F,2)
r0 = ones(M,M)
z0 = eye(M)
z0[1,1] = 0
z0[M,M] = 0
w = 2.0/(1+sin(h*pi))

for k in M-1:-1:2, j in M-1:-1:2
    # Backward sweep of SOR
    z0[j,k] = w*( 0.25 * (z0[j-1,k] + z0[j+1,k] + z0[j,k-1] + z0[j,k+1] - h^(2.0) * r0[j,k])) + (1-w)*z0[j,k]
  end


for j in 2:M-1, k in 2:M-1
      # Forward Sweep SOR
      z0[j,k] = w*( 0.25 * (z0[j-1,k] + z0[j+1,k] + z0[j,k-1] + z0[j,k+1] - h^(2.0) * r0[j,k])) + (1-w)*z0[j,k]
  end



#=
println(k)
println(abs(u - uturtle))
println(abs(r - residual))

b = abs(rand(9,9))
B = .5(b' + b) + 100*eye(9)
eig(B)
=#


#= testing vcycle
B = ones(10,10)
M = 10
rb = 2:M-1
red = 1
for j in 2:M-1
red = filter(k -> iseven(k+j), 2:M-1)
#println(red)
end

include("PDEtool.jl")
using PDEtool
####### Set up our R.H.S. function
turtles = 10
h = 2.0^(-turtles)
funcRHS = (x,y) -> -exp(-(x - 0.25)^2 - (y - 0.6)^2)
mesh, F = PDEtool.h_space(funcRHS, h)
# Set up  initial guess and boundary data, this is homogenous Dirichlet
u0 = zeros(size(F))
v = [1 3]
tol = 10.0^(-4)
maxiter = 100
#uturtle = PDEtool.vcycle(u0, F, turtles, v)
uturtle, residual, iter = @time PDEtool.multigrid(u0, F, maxiter, tol, turtles, v)
uSOR, iterst = PDEtool.SOR(h, 30, tol, F)
iterdiff = vecnorm(uSOR-uturtle)


uSOR, iterst = SOR(h, 30, tol, F)
Ua, residual = MG_vcycle(h, F, u, turtles, s1, s2, tol)
uout, res, maxiter = gauss_sidel(u, h, 1000, tol, 1, F)
iterdiff = vecnorm(uout-Ua[turtles],1)

uSOR, iterst = SOR(h, 1000, tol, F)
iterdiff = vecnorm(uSOR-Ua[turtles],1)
iterdiff = vecnorm(uSOR-uout,1)

c = squish(ones(3,3))
d = squish(c)
q = foomp(ones(size(F)))

for turtles in [5 6 7]
  h = 2.0^(-turtles)
  funcRHS = (x,y) -> -exp(-(x - 0.25)^2 - (y - 0.6)^2)
  mesh, F = PDEtool.h_space(funcRHS, h)
  # Set up  initial guess and boundary data, this is homogenous Dirichlet
  u = zeros(size(F))
  s1 = 1
  s2 = 3
  tol = 10.0^(-4)
  Ua, residual = MG_vcycle(h, F, u, turtles, s1, s2, tol)
end

h = 2.0^(-2)
funcRHS = (x,y) -> -exp(-(x - 0.25)^2 - (y - 0.6)^2)
mesh, F = PDEtool.h_space(funcRHS, h)
u = zeros(size(F))
u, res, maxiter = gauss_sidel(u, h, 50, .000001, 1, F)
=#
