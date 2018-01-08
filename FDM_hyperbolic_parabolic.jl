#=
Kaela S. Vogel
December 2017

This is a collection of PDE solvers using finite difference methods

=#


#=
meshmaker is a basic function to set up a mesh for the finite difference methods below

=#

function meshmaker(funcRHS::Function, h::Float64, k::Float64)
hmesh = [j for j in 0:h:1]
kmesh = [j for j in 0:k:1]
M = length(hmesh)
T = length(kmesh)
F = zeros(M,T)

for k in 1:T, j in 1:M,
  F[j,k] = funcRHS(hmesh[j], kmesh[k])
end

return hmesh, kmesh, F, M, T

end #function]


#= CRANK-NICOLSON SOLVER FOR THE HEAT EQUATION
This is a solver for the heat equation  for 0 < x < 1, 0 <= t <= 1
u_t(x,t) = a * u_xx(x,t) + F(x,y)
u(0,t) = u_0t
u(1,t) = u_1t
u(x,0) = u_x0 initial condtions
v - output solution matrix where rows are x values at column time t
h - space step
k - time step
funcRHS - function of x, t, and constants
spaceDim - amount of spatial dimensions

Crank-Nicolson is an implicit methods and requires solving a matrix at each time step
=#


function CN_heat(h::Float64, k::Float64, funcRHS::Function,  u_0t::Function, u_1t::Function, u_x0::Function, a::Float64)

##############################################################

#################### Set up the mesh and th e solution vector

hmesh, kmesh, F, M, T = meshmaker(funcRHS, h, k)
v = zeros(M,T) # time stepping vector
rhs = zeros(M,1) # temp right hand side array
v[1:M, 1] = u_x0(hmesh) # initial condtions

v[1, 1:T] = u_0t(kmesh) # boundary condition 1
v[M, 1:T] = u_1t(kmesh) # boundary condition 2

#################### Set up the matrix to solve at each time step

r = (a*k)/(2*h^2)
du = -r * ones(M-1)
d  = (1+2*r) * ones(M)
A = Tridiagonal(du, d, du)

################### Run the solver for time steps 2 to M

for i = 1: (T-1) # time stepping loop

    # boundary condition t = 0 setup for RHS vector
    rhs[1] = r*(u_0t(kmesh[i]) + u_0t(kmesh[i+1])) + (1-2*r)*v[1,i] + r*v[2,i]

    # form the right hand side
    for j = 2: M-1
        rhs[j] = r*v[j-1, i] + (1 - 2*r)*v[j, i] + r*v[j+1, i] + (k/2)*(F[j, i]+F[j, i+1])
    end # end rhs setup loop

    # boundary condition t = 1 setup for RHS vector
    rhs[M] = r*(v[M-1,i]) + (1-2r)*v[M,i] + r*(u_1t(kmesh[i]) + u_1t(kmesh[i+1]))

    # Solve the system at each timestep
    v[1:M, i+1] = A\rhs

end # end time step loop


    return v
end
