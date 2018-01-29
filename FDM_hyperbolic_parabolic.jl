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

    # This evaluates F at all the inner points and the last time boundary point
    for k in 2:T, j in 2:M-1
        F[j,k] = funcRHS.(hmesh[j], kmesh[k])
    end

return (hmesh, kmesh, F, M, T)

end #function]


#= CRANK-NICOLSON SOLVERS FOR THE HEAT EQUATION
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
using PDEtool

####################### Crank Nicolson, 2 spatial dimensions
# So rather than save a million matrices I will only save the one for the final timestep

function CN_heat2D(h::Float64, k::Float64, funcRHS::Function,  u_0t::Function, u_1t::Function, u_xy0::Function, a::Float64)


    #################### Set up the mesh and th e solution vector
    # The solution vector v has increasing y-space going down the colums and increasing x-space going right across columns

    hmesh, kmesh, F, M, T = meshmaker(funcRHS, h, k)
    v = zeros(M,M) # numeric solution, only latest timestep
    rhs = zeros(M-2,M-2) # temp right hand side matrix
    lhs = zeros(1,(M-2)*(M-2)) # temp left hand side hand side vector

    ################### Set up inital and boundary conditions in solution matrix

    # Initial conditions
    for j in 1:M, k in 1:M
        v[j,k] = uxy0(hmesh[j], hmesh[k])
    end

    # Boundary conditions, this function assume boundary is 0
    v[1:M, 1] = 0
    v[1, 1:M] = 0
    v[M, 1:M] = 0
    v[1:M, M] = 0

    #################### Set up the matrix to solve at each time step
    # Note: We only use this matrix to solve for the inner points, so it
    # has dimensions (M-2)(M-2) X M-2)(M-2).

    # This is for the LHS, it requires you to form the input spatial matrix into a single long vector.
    r = (a*k)/(h^2)
    du = 1.0 * ones(M-3)
    d  = - 2.0 * ones(M-2)
    L = Tridiagonal(du, d, du) # Form the 1D laplacian matrix
    L2D = kron(eye(M-2), L) + kron(L, eye(M-2)) # Form 2D Laplacian
    A = eye((M-2)*(M-2)) - (r/2)*L2D
    A = factorize(A)

    # This is for the %HS, it requires you to form the input spatial matrix into a single long vector.
    B = eye((M-2)*(M-2)) + (r/2)*L2D
    B = factorize(B)


    ################### Run the solver for time steps 2 to T
    # This solver uses the 2D laplacian which requires the (M-2)X (M-2) be temporarily transformed to a 1 X(M-2)(M-2) vector in order to do matrix multiplication.

    for i = 2: (T-1) # time stepping loop
        # I'm just ignoring boundar condition setup because the boundary is 0.  This is the part where we form (1-r/2 L)v^n.
        rhs = copy[v[2:M-1,2:M-1]]
        rhs = reshape(matrix,(M-2)*(M-2), 1)
        rhs = B*rhs # This will be slow, I don't care

        # I'm just ignoring boundar condition setup because the boundary is 0.  This is the part where we solve (1+r/2 L)v^n+1 = rhs
        lhs = A\rhs
        v[2:M-1, 2:M-1] = reshape(lhs, M-2, M-2)

    end


    return v

end #function




####################### Crank Nicolson, 1 for i = 2: (T-1) # time stepping loopspatial dimension

function CN_heat(h::Float64, k::Float64, funcRHS::Function,  u_0t::Function, u_1t::Function, u_x0::Function, a::Float64)


    #################### Set up the mesh and th e solution vector
    # The solution vector v has increasing space as going down rows and increasing time as going right in columns

    hmesh, kmesh, F, M, T = meshmaker(funcRHS, h, k)
    v = zeros(M,T) # numeric solution
    rhs = zeros(M-2,1) # temp right hand side array
    v[1:M, 1] = u_x0.(hmesh) # initial condtions

    v[1, 1:T] = u_0t.(kmesh) # boundary condition 1
    v[M, 1:T] = u_1t.(kmesh) # boundary condition 2

    #################### Set up the matrix to solve at each time step
    # Note: We only use this matrix to solve for the inner points, so it
    # has dimensions (M-2)x(M-2)
    r = (a*k)/(2*h^2)
    du = -r * ones(M-3)
    d  = (1+2*r) * ones(M-2)
    A = Tridiagonal(du, d, du)
    A = factorize(A)

    ################### Run the solver for time steps 2 to M
    # Be careful with indexing of v because of avoiding boundary, julia starts indexing at 1, so the boundary is located at index 1 and M, the M-2 inner points are located in index 2:M-1.  We are still indexing w.r.t. the original grid when we build rhs, so relative to the original grid, rhs is a M-2x1 vector
    # Also be careful that F is a forward temporal average and is not present at time step 1
    # The loop indexing is off by one for what actual time step index, but while F is MxT it is only evaluated at the inner points
    for i = 1: (T-1) # time stepping loop

        # This is the first equation for the inner point
        rhs[1] = r*(u_0t(kmesh[i]) + u_0t(kmesh[i+1])) + (1-2*r)*v[2,i] + r*v[3,i] + (k/2)*(F[2, i]+ F[2, i+1]) #k*F[2, i] #

        # form the right hand side
        for j = 2: M-3
            rhs[j] = r*v[j, i] + (1 - 2*r)*v[j+1, i] + r*v[j+2, i] +  +(k/2)*(F[j+1, i]+F[j+1, i+1]) #k*F[j+1, i] #
        end # end rhs setup loop

        # equation for the last inner point
        rhs[M-2] = r*(v[M-2,i]) + (1-2r)*v[M-1,i] + r*(u_1t(kmesh[i]) + u_1t(kmesh[i+1])) +(k/2)*(F[M-1, i]+F[M-1, i+1]) # k*F[M-1,i] #

        # Solve the system at each timestep
        v[2:M-1, i+1] = A\rhs

    end # end time step loop


    return v
end #function

###########################################################################

#=
This function performs BDF-2 on the heat equation
which has the form
u^{n+1} = B_1 u^{n} + B_2 u^{n-1} + b
=#
function BDF2_heat(h::Float64, k::Float64, funcRHS::Function,  u_0t::Function, u_1t::Function, u_x0::Function, a::Float64)


    #################### Set up the mesh and th e solution vector

    hmesh, kmesh, F, M, T = meshmaker(funcRHS, h, k)
    v = zeros(M,T) # numeric solution
    rhs = zeros(M-2,1) # temp right hand side array
    BEstep = zeros(M-2,1) # temp right hand side array for BE
    v[1:M, 1] = u_x0.(hmesh) # initial condtions

    v[1, 1:T] = u_0t.(kmesh) # boundary condition 1
    v[M, 1:T] = u_1t.(kmesh) # boundary condition 2

    #################### Set up the laplacian maatrix and the transform it to match BDF2 form. # Note: We only use this matrix to solve for the inner points, so it has dimensions (M-2)x(M-2)

    r = (2*a*k)*(1/h^2)
    du = 1.0 * ones(M-3)
    d  = - 2.0 * ones(M-2)
    L = Tridiagonal(du, d, du)
    A = 3.0 * eye(M-2) - r*L
    A = factorize(A)

    ################### Run the solver for time steps 3 to T

    # BDF2 requires two past points, so we need to run backwaards Eular for the second time step

    # Set up Backwards Euler matrix
    du =  1.0 * ones(M-3)
    d  = - 2.0 * ones(M-2)
    L = Tridiagonal(du, d, du)
    B = eye(M-2) - r*L
    B = factorize(B)

    # Deal with boundary conditions
    BEstep[1,1] = r*v[1, 2] + v[2,1]
    BEstep[M-2,1] = r*v[M, 2] + v[M-1,1]

    # Inner, inner points for the rhs
    BEstep[2:M-3] = v[3:M-2, 1]

    # Run Backward Euler for 1 step to get step 2
    v[2:M-1, 2] = B\BEstep



    for i = 2: (T-1) # time stepping loop

        # Deal with boundary points, since g_0(t+1), g_1(t+1) is known
        rhs[1, 1] = r*v[1, i+1]+ 4*v[2, i] - v[2, i-1]
        rhs[M-2, 1] = r*v[M, i+1]+ 4*v[M-1, i] - v[M-1, i-1]

        # Put together the right hand side vector
        rhs[2:M-3,1] = 4*v[3:M-2, i] - v[3:M-2, i-1]

        # Solve the system at each timestep
        v[2:M-1, i+1] = A\rhs

    end # end time step loop


    return v
end #function



############################################################################

#=
This function is to set up a mesh refinement study for a numerically solved PDE.  It accomplishes this by restricting the finer mesh solution to a coarser grid by taking every 2i+1 index .  Then the ratio of discrete l1 norm of the difference at sucessive timesteps is taken.

=#


function refinement_ratios(v_unrefined::Vector{Any}, hs::Array{Float64,1}, ks::Array{Float64,1})
# Typically we only compare one time column
# Set up the grid restriction operator
# WARNING: This if for grid spacing using powers of 2, e.g. h = 1/2^6
# Restriction operator takes every 2i + 1 index of the finer grid

    v_restricteds = Vector{Any}(length(v_unrefined)-1) # This will store the newly restricted solutions v

    # Restrict every MxM grid to a (M-1)/2 X (M-1)/2 grid
    for j in 2: length(v_unrefined)
        squishing = copy(v_unrefined[j])
        newgridlength = Int64((size(squishing,1)-1)/2)
        newgridwidth = Int64((size(squishing,2)-1)/2)
        restricted = squishing[[2*i+1 for i in  0:1:newgridlength], [2*i+1 for i in  0:1:newgridwidth]]
        v_restricteds[j-1] = copy(restricted)
    end # restricting


    # Calculate the ratios
    refinedratio_length = (length(v_restricteds)-1)
    refinedratio = Vector{Float64}(refinedratio_length)

    for j in 1: (refinedratio_length)
        # The vecnorm is the appropriate norm here because we want to compare entries, we are actual doing the norm of a discrete integration

        refinedratio[j] = ( hs[j]*ks[j]* vecnorm( v_unrefined[j] - v_restricteds[j], 1) )/( hs[j+1]*ks[j+1]*vecnorm( v_unrefined[j+1] - v_restricteds[j+1],1))

        # refinedratio[j] = ( hs[j]*ks[j]* vecnorm( v_unrefined[j] - v_restricteds[j], 1) )/( hs[j+1]*ks[j+1]*vecnorm( v_unrefined[j+1] - v_restricteds[j+1],1))

    end

    return (refinedratio, v_restricteds)

end #function


# When doing a norm to see order of accuracy, don't use a regular norm, treat it as a discrete integral in time (or space) and don't forget to add the weights

# "waterfall plot", or look at a few time snapshots in space, heat maps are hard to process.
#  Discountinuous = lots of high frequency spatial components
# BDF2 is centered at the the n+1 time point (so evaluate F at n+1)
# Direct solve may be expensive for the large 2d spatial systems
# 2nd homework use a sparse solver
# 228A look at direct solve for
# ADI solve is a bunch of 1-d direct solves.
