#=
Kaela S. Vogel
December 2017

This is a collection of PDE solvers using finite difference methods

=#


#=
meshmaker is a basic function to set up a mesh for the finite difference methods below

=#

function meshmaker(funcRHS::Function, h::Float64, k::Float64, s::Bool)

    # s is a boolean, True gives normal mesh, false give cell-centered
    if s
        hmesh = [j for j in 0:h:1]
        kmesh = [j for j in 0:k:1]


    else
        hmesh = [j for j in 0:h:(1-h)] .+ 0.5*h
        kmesh = [j for j in 0:k:1]

    end
    M = length(hmesh)
    T = length(kmesh)
    F = zeros(M,T)

    # This evaluates F at all the inner points and the last time boundary point
    for i in 2:T, j in 2:M-1
        F[j,i] = funcRHS.(hmesh[j], kmesh[i])
    end

return (hmesh, kmesh, F, M, T)

end #function]

###############################################################################
# Begin Hyperbolic Solvers
###############################################################################

#= advectionMAT makes the stepping matrices for advectionFDM, it takes values for the diagonal, off diagonal, and corners and builds a N x N matrix where N is the number of mesh points

c1                  Value for upper right corner
c2                  Value for lower left corner
dc                   Value repeated along diagonal
duc                  Value repeated along superdiagonal
dlc                  Value repeated along subiagonal
N                   Dimension of the matrix
=#

function advectionMAT(c1::Float64, c2::Float64, dlc::Float64, dc::Float64, duc::Float64, N::Int64)

    du = duc * ones(N-1)
    d  = dc* ones(N)
    dl = dlc*ones(N-1)
    Amat = eye(N) *Tridiagonal(dl, d, du)
    Amat[1,end] = c1
    Amat[end,1] = c2

    return sparse(Amat)

end # function


#= advectionFDM can call several finite difference methods for solving the advection equation. Options

FDM = 1              Lax-Wendroff
FDM = 2              Upwinding
FDM = 3              Crank-Nicolson

=#

function advectionFDM(h::Float64, k::Float64, Nx::Int64, Nt::Int64, u_x0::Function, a::Float64, s::Bool, FDM::Int64)

###################### Initial setup of initial conditions and mesh

# I am assuming periodic conditions which is taken care of in the upper right and lower left corners of the time stepping matrix and so no boundary conditions are needed

    hmesh = collect(linspace(0, 1-h, Nx))
    v = zeros(Nx, Nt+1) # numeric solution
    v[:, 1] = u_x0.(hmesh) # initial condtions
    nu = a*k/h
    error = zeros(3,1) # Vector to hold the error in the 1, 2, and infinity norms

###################### Switch between stepping methods

    if FDM == 1 # Step with Lax-Wendroff
        # v_{j}^{n+1} = v_j^n - \frac{ak}{2h}(v_{j+1}^n - v_{j-1}^n) + \frac{a^2k^2}{2h^2}(v_{j+1}^n - 2v_j^n + v_{j-1}^n).
        # Generate MxM stepping matrix
        A = advectionMAT(nu/2 + nu^2/2, -nu/2 + nu^2/2, nu/2 + nu^2/2, 1 - nu^2, -nu/2+ nu^2/2, Nx)

        # Timestepping
        for i = 1: (Nt) # time stepping loop
            # Advance one timestep
            v[:, i+1] = A*v[:, i]
        end # end time step loop

        # This is a fix for the fact the the interval does not end at exactly one
        #corr =( kmesh[end]-1)/.9
        #v0 = u_x0.(hmesh.-a*kmesh[end]) # .- kmesh[end]

        error[1] = h*vecnorm((v[:,end] - v[:,1]), 1)
        error[2] = h^(1/2)*vecnorm((v[:,end] - v[:,1]), 2)
        error[3] = vecnorm((v[:,end] - v[:,1]), Inf)

        return (v, error)
#####################


    elseif FDM == 2 # Step with Upwinding
        # Generate MxM stepping matrix, this is for a >= 0

        A = advectionMAT(nu, 0.0, nu, 1 - nu, 0.0, Nx)

        # Timestepping
        for i = 1: (Nt) # time stepping loop
            # Advance one timestep
            v[:, i+1] = A*v[:, i]
        end # end time step loop

        error[1] = h*vecnorm((v[:,end] - v[:,1]), 1)
        error[2] = h^(1/2)*vecnorm((v[:,end] - v[:,1]), 2)
        error[3] = vecnorm((v[:,end] - v[:,1]), Inf)

        return (v, error)

#####################

    elseif FDM == 3 # Step with Crank-Nicolson
        # Generate MxM stepping matrix, this is for the LHS
        A = advectionMAT(-nu/4, nu/4, 1.0, nu/4, -nu/4, Nx)
        # Generate MxM stepping matrix, this is for the RHS
        B = A'

        # Timestepping
        for i = 1: (Nt) # time stepping loop
            # Advance one timestep
            v[:, i+1] = A\(B*v[:, i])
        end # end time step loop

        error[1] = h*vecnorm((v[:,end] - v[:,1]), 1)
        error[2] = h^(1/2)*vecnorm((v[:,end] - v[:,1]), 2)
        error[3] = vecnorm((v[:,end] - v[:,1]), Inf)

        return (v, error)

#####################

    else
        error("No method type selected")
    end



end #function
###############################################################################
# Begin Parabolic Solvers
###############################################################################

#= FitzHugh-Nagumo equations
This is a solver for the FitzHugh-Nagumo equations. It uses a fractional step method with Peaceman-Rachford ADI for the diffusion solve and Backward-Euler for the reaction solve. This basically works out to making a Backward-Euler sandwhich with the two steps of PR-ADI done over a half time step as bread.


Spatial dimensions on unit square 0 < x < 1, 0 < y < 1
u_t(x,t) = a * u_xx(x,t) + F(x,y)
u(0,t) = u_0t
u(1,t) = u_1t
u(x, y, 0) = u_x0 initial condtions
v - output solution matrix where rows are x values at column time t
h - space step
k - time step
funcRHS - function of x, t, and constants
spaceDim - amount of spatial dimensions

=#

function FHN2D(h::Float64, k::Float64, v_xy0::Function, w_xy0::Function, Tend::Int64, s::Bool)

    #################### Parameters specific to The FitzHugh-Nagumo equations
    a = 0.1
    gam = 2.0
    eps = 0.005
    I = 0.0
    D = 0.00001
    fv_vw(v,w) =  (a - v).*(v - 1).*v - w + I # Used for F.E. solve
    fw_vw(v,w) = eps.*(v - gam.*w) # Used for F.E. solve

    #################### Set up the mesh and the solution vector
    hmesh = [j for j in 0:h:(1-h)] .+ 0.5*h
    kmesh = [j for j in 0:k:Tend]
    M = length(hmesh)
    T = length(kmesh)
    v = zeros(M,M) # numeric solution
    w = zeros(M,M) # numeric solution
    stor = 34 # after how many timesteps to store a solution
    numsol = Int64(ceil(T/stor)+1)
    sol_cnt = 2 # counter for how many solutions have been stored
    vstored = Vector{Any}(numsol)

    # Intermediate variables pre-allocation
    v_star = zeros(M,M) # numeric solution after diffusing in x-direction
    RHS = zeros(M,M)
    LHS = zeros(M,M)
    v_half = zeros(M,M) # numeric solution after PR-ADI timestep k/2
    v_FE = zeros(M,M) # numeric solution after FE step
    w_FE = zeros(M,M) # numeric solution after FE step

    ################### Set up inital conditions in solution matrix

    # Initial conditions
    for i in 1:M, j in 1:M
        v[i,j] = v_xy0(hmesh[i], hmesh[j])
        w[i,j] = w_xy0(hmesh[i], hmesh[j])
    end
    vstored[1] = v # Store the initial voltage

    #################### Set up the left and right matrices to solve at each time step.

    r = (D*k)/(2*2*h^2) # The extra factor of 2 is due to the half-time step
    du = 1.0 * ones(M-1)
    d  = - 2.0 * ones(M)
    d[1]  = -1.0
    d[end] = -1.0
    L =  Tridiagonal(du, d, du) # Form the 1D laplacian matrix

    # This is for the LHS for the implicit stepping of the form (I - ak/2 L)
    A = speye(M) - r*L

    # This is for the RHS for the explicit stepping of the form (I + ak/2 L)
    B = speye(M) + r*L

    for i in 2:T-1

        ###################################### PR-ADI solve k/2 timestep
        # Diffuse in the x-direction (I - ka/2 Lx) u* = (I + ka/2 Ly) u^n
        RHS = v*B
        v_star = A\RHS
        # Diffuse in the y-direction (I - ka/2 Ly) u^(n+1) = (I + ka/2 Lx) u*
        RHS = (B*v_star)'
        LHS = A\RHS
        # Assign v^(*)
        v_half = LHS'

        ###################################### Forward Euler solve
        v_FE = v_half + k*fv_vw(v_half,w)
        w = w + k*fw_vw(v_half,w)

        ###################################### PR-ADI solve k/2 timestep
        # Diffuse in the x-direction (I - ka/2 Lx) u* = (I + ka/2 Ly) u^n
        RHS = v_FE*B
        v_star = A\RHS
        # Diffuse in the y-dvstored[sol_cnt] = virection (I - ka/2 Ly) u^(n+1) = (I + ka/2 Lx) u*
        RHS = (B*v_star)'
        LHS = A\RHS
        v = LHS'

        if mod(i, stor) == 0 #This stores a solution ever
            vstored[sol_cnt] = v
            sol_cnt += 1
        end

        #vstored[i] = v
    end # time step loop
    vstored[end] = v

    return vstored

end # function

#=

=#
#= PEACEMAN-RACHFORD ADI FOR THE 2D HEAT EQUATION
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

=#

function prADI_heat2D(h::Float64, k::Float64, funcRHS::Function,  u_0t::Function, u_1t::Function, u_xy0::Function, a::Float64, s::Bool)

    #################### Set up the mesh and th e solution vector

    hmesh, kmesh, F, M, T = meshmaker(funcRHS, h, k, s)
    v = zeros(M,M) # numeric solution, only latest timestep

    # Intermediate variables pre-allocation
    v_star = zeros(M,M) # numeric solution after diffusing in x-direction
    RHS = zeros(M,M)
    LHS = zeros(M,M) # numeric solution, only latest timestep


    ################### Set up inital conditions in solution matrix

    # Initial conditions
    for i in 1:M, j in 1:M
        v[i,j] = u_xy0(hmesh[i], hmesh[j])
    end

    #################### Set up the left and right matrices to solve at each time step.

    r = (a*k)/(2*h^2)
    du = 1.0 * ones(M-1)
    d  = - 2.0 * ones(M)
    d[1]  = -1.0
    d[end] = -1.0
    L =  Tridiagonal(du, d, du) # Form the 1D laplacian matrix

    # This is for the LHS for the implicit stepping of the form (I - ak/2 L)
    A = speye(M) - r*L

    # This is for the RHS for the explicit stepping of the form (I + ak/2 L)
    B = speye(M) + r*L

    for i in 2:T-1
        # In the P-R ADI method we diffuse in the x and y directions seperately and     therefore end up with 2 1D systems

        # Diffuse in the x-direction (I - ka/2 Lx) u* = (I + ka/2 Ly) u^n
        RHS = v*B
        v_star = A\RHS

        # Diffuse in the y-direction (I - ka/2 Ly) u^(n+1) = (I + ka/2 Lx) u*
        RHS = (B*v_star)'
        LHS = A\RHS

        # Assign v^(n+1)
        v = LHS'

    end # time step loop

    return v

end # function

####################### Forward Euler, 2 spatial dimensions
# So rather than save a million matrices I will only save the one for the final timestep.  Requires k <= h^2/2a

function FE_heat2D(h::Float64, k::Float64, funcRHS::Function,  u_0t::Function, u_1t::Function, u_xy0::Function, a::Float64, s::Bool)


    #################### Set up the mesh and th e solution vector

    hmesh, kmesh, F, M, T = meshmaker(funcRHS, h, k, s)
    v = zeros(M,M) # numeric solution, only latest timestepA = kron(Iy, Lx) + kron(Ly,Ix)

    ################### Set up inital and boundary conditions in solution matrix

    # Initial conditionsA = kron(Iy, Lx) + krA = kron(Iy, Lx) + kron(Ly,Ix)on(Ly,Ix)
    for j in 1:M, i in 1:M
        v[j,i] = u_xy0(hmesh[j], hmesh[i])
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
    L =  Tridiagonal(du, d, du) # Form the 1D laplacian matrix
    # Am I doing kroneker delta in right place
    L2D = kron(speye(M-2), sparse(L)) + kron(sparse(L), speye(M-2))
    A = speye((M-2)*(M-2)) + r*L2D

    ################### Run the solver for time steps 2 to T
    # This solver uses the 2D laplacian which requires the (M-2)X (M-2) be temporarily transformed to a 1 X(M-2)(M-2) vector in order to do matrix multiplication.

    for i = 2: (T-1) # time stepping loop
        # I'm just ignoring boundar condition setup because the boundary is 0.
        rhs = v[2:M-1,2:M-1]
        rhs = reshape(rhs,(M-2)*(M-2), 1) # julia reshapes columnwise
        lhs = A*rhs # This will be slow, I don't care

        v[2:M-1, 2:M-1] = reshape(lhs, M-2, M-2)

    end


    return v


end #function

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

####################### Crank Nicolson, 2 spatial dimensions
# So rather than save a million matrices I will only save the one for the final timestep

function CN_heat2D(h::Float64, k::Float64, funcRHS::Function,  u_0t::Function, u_1t::Function, u_xy0::Function, a::Float64, s::Bool)


    #################### Set up the mesh and th e solution vector


    hmesh, kmesh, F, M, T = meshmaker(funcRHS, h, k, s)
    v = zeros(M,M) # numeric solution, only latest timestepA = kron(Iy, Lx) + kron(Ly,Ix)
    #rhs = zeros(M-2,M-2) # temp right hand side matrix
    lhs = zeros(1,(M-2)*(M-2)) # temp left hand side hand side vector

    ################### Set up inital and boundary conditions in solution matrix

    # Initial conditionsA = kron(Iy, Lx) + kron(Ly,Ix)
    for j in 1:M, i in 1:M
        v[j,i] = u_xy0(hmesh[j], hmesh[i])
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
    L =  Tridiagonal(du, d, du) # Form the 1D laplacian matrix
    # Am I doing kroneker delta in right place
    L2D = kron(speye(M-2), sparse(L)) + kron(sparse(L), speye(M-2))
    A = speye((M-2)*(M-2)) - (0.5*r)*L2D
    A = factorize(A)

    # This is for the %HS, it requires you to form the input spatial matrix into a single long vector.
    B = speye((M-2)*(M-2)) + (0.5*r)*L2D




    ################### Run the solver for time steps 2 to T
    # This solver uses the 2D laplacian which requires the (M-2)X (M-2) be temporarily transformed to a 1 X(M-2)(M-2) vector in order to do matrix multiplication.

    for i = 2: (T-1) # time stepping loop
        # I'm just ignoring boundar condition setup because the boundary is 0.  This is the part where we form (1-r/2 L)v^n.
        lhs = zeros(1,(M-2)*(M-2)) # temp left hand side hand side vector
        #rhs = v[2:M-1,2:M-1]
        #rhs = reshape(rhs,(M-2)*(M-2), 1) # julia reshapes columnwise
        #rhs = B*rhs # This will be slow, I don't care
        rhs = zeros(M,M) # temp right hand side matrix
        # form the right hand side

        for j = 2: M-1, i = 2:M-1
            rhs[j,i] = (1 - 2*r)*v[j, i] + 0.5*r*(v[j+1, i] +v[j, i+1] + v[j-1, i] + v[j, i-1])
        end # end rhs setup loop
        moop = rhs[2:M-1,2:M-1]
        moop = reshape(moop,(M-2)*(M-2), 1) # julia reshapes columnwise

        # I'm just ignoring boundar condition setup because the boundary is 0.  This is the part where we solve (1+r/2 L)v^n+1 = rhs
        lhs = A\moop

        #lhs = gauss_sidelCN(v, 1000, 10^(-7.0), r, rhs)


        v[2:M-1,2:M-1] = reshape(lhs, M-2, M-2) #

    end


    return v


end #function




####################### Crank Nicolson, 1 for i = 2: (T-1) # time stepping loopspatial dimension

function CN_heat(h::Float64, k::Float64, funcRHS::Function,  u_0t::Function, u_1t::Function, u_x0::Function, a::Float64, s::Bool)


    #################### Set up the mesh and th e solution vector
    # The solution vector v has increasing space as going down rows and increasing time as going right in columns

    hmesh, kmesh, F, M, T = meshmaker(funcRHS, h, k, s)
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
function BDF2_heat(h::Float64, k::Float64, funcRHS::Function,  u_0t::Function, u_1t::Function, u_x0::Function, a::Float64, s::Bool)


    #################### Set up the mesh and th e solution vector

    hmesh, kmesh, F, M, T = meshmaker(funcRHS, h, k, s)
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


function refinement_ratios(v_unrefined::Vector{Any}, hs::Array{Float64,1}, ks::Array{Float64,1}, s::Bool)
# Typically we only compare one time column
# Set up the grid restriction operator

# s is a Bool that switches based on whether a regular or cell centered mesh is used, s = Bool(0) gives cell centered

    v_restricteds = Vector{Any}(length(v_unrefined)-1) # This will store the newly restricted solutions v

    if s # Restricts a generic mesh
        # WARNING: This if for grid spacing using powers of 2, e.g. h = 1/2^6
        # Restriction operator takes every 2i + 1 index of the finer grid

        # Restrict every MxM grid to a (M-1)/2 X (M-1)/2 grid
        for j in 2: length(v_unrefined)
            squishing = copy(v_unrefined[j])
            newgridlength = Int64((size(squishing,1)-1)/2)
            newgridwidth = Int64((size(squishing,2)-1)/2)
# Typically we only compare one time column
            restricted = squishing[[2*i+1 for i in  0:1:newgridlength],[2*i+1 for  i in  0:1:newgridwidth]]
            v_restricteds[j-1] = copy(restricted)
        end # restricting
# Typically we only compare one time column

    else # Restricts a cell-centered mesh
        # WARNING: This if for grid spacing using powers of 3, e.g. h = 1/3^6
        # Restriction operator takes every 2i + 1 index of the finer grid

        # Restrict every MxM grid to a (M)/3 X (M)/3 grid
        for j in 2: length(v_unrefined)
            squishing = copy(v_unrefined[j])
            newgridlength = Int64(size(squishing,1)/3)
            newgridwidth = Int64(size(squishing,2)/3)
            restricted = squishing[[3*i+2 for i in  0:1:newgridwidth-1],[3*i+2 for  i in  0:1:newgridwidth-1]]
            v_restricteds[j-1] = copy(restricted)
        end # restricting

    end #conditional switch

    # Calculate the ratios
    refinedratio_length = (length(v_restricteds)-1)
    refinedratio = Vector{Float64}(refinedratio_length)

    for j in 1: (refinedratio_length)
        # The vecnorm is the appropriate norm here because we want to compare entries, we are actual doing the norm of a discrete integration

        refinedratio[j] = ( hs[j]^2*ks[j]* vecnorm( v_unrefined[j] - v_restricteds[j], 1) )/( hs[j+1]^2*ks[j+1]*vecnorm( v_unrefined[j+1] - v_restricteds[j+1],1))

# hs[j]*ks[j]*
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

##################################################################

####### Gauss-Siedel iteration specifically for Crank-Nicolson system

function gauss_sidelCN( v::Array{Float64,2}, maxiter::Int64, tol::Float64, r::Float64, F::Array{Float64,2})
    u = copy(F)
    M = size(u, 2)
    # Run the iterations
    V = zeros(size(u))
    s = 1/(1+2r)


    for iter in 0:maxiter


    for j in 2:M-1, k in 2:M-1
      u[j,k] =  s* (0.5*r*u[j-1,k] + 0.5*r*u[j+1,k] + 0.5*r*u[j,k-1] + 0.5*r*u[j,k+1] + F[j,k])
    end


    iterdiff = vecnorm(u-V,1)


    #if vecnorm(residual,1) < tol*vecnorm(F,1)

    if iterdiff < tol*vecnorm(u,1)
    #println("GS Tolerance reached after $iter iterations")

    return u
    end

    V = copy(u)
    end # iteration loop
    println("Tolerance not reached")

    return u

    #####
end # function
