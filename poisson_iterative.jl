# Kaela S. Vogel
# This is a module of some of the different iterative methods to solve the Poisson Equation
# To call it you need "using(poisson_iterative.jl)"
# h - mesh spacing
#



####### Set up our R.H.S. function
funcRHS(x,y) = -exp(-(x - 0.25)^2 - (y - 0.6)^2)

####### Function to set up mesh and RHS for a given spacing h
function h_space(h)
  mesh = [j for j in 0:h:1]
  M = length(mesh)
  F = zeros(M,M)

  for j in 1:M, k in 1:M
      F[j,k] = funcRHS(mesh[j],mesh[k])
  end

  return mesh,F

end #function



####### Jacobi Iteration
# use differance between successive iterations for tolerance
function jacobi_iter(h, maxiter, tol)
  # Set up mesh and F (right hand side)
  mesh, F = h_space(h)

  # Set up  initial guess and b oundary data, this is homogenous Dirichlet
  u = zeros(size(F))
  u[1,:] = zeros(1,M)
  u[M,:] = zeros(1,M)

  for iter in 0:maxiter
    
    for j in 2:M-1, k in 2:M-1
        V[j,k] = 0.25 * (u[j-1,k] + u[j+1,k] + u[j,k-1] + u[j,k+1] - h^(2.0) * F[j,k])
    end

    iterdiff = vecnorm(u - V,1)
    if iterdiff < tol*vecnorm(V)
        #println("Tolerance reac hed after $iter iterations")
      return V, iter
    end
    u = copy(V)

  end
  #println("Tolerance not reached")
  return u, maxiter

  #####
end # function

##################################################################

####### Gauss-Siedel iteration

function gauss_sidel(h::Float64, maxiter::Int64, tol::Float64, rb::Bool)
  # Set up mesh and F (right hand side)
  mesh = [j for j in 0:h:1]
  M = length(mesh)
  F = zeros(M,M)

  for j in 1:M, k in 1:M
      F[j,k] = funcRHS(mesh[j],mesh[k])
  end


  # Set up  initial guess and boundary data, this is homogenous Dirichlet
  u = zeros(size(F))

  # Run the iterations
  V = zeros(size(u))

  # Indices for red/black version
  red = filter(k -> iseven(k+j), 2:M-1)
  black = filter(k -> isodd(k+j), 2:M-1)

  for iter in 0:maxiter

    if rb
      for j in 2:M-1
        red = filter(k -> iseven(k+j), 2:M-1)
        black = filter(k -> isodd(k+j), 2:M-1)
          for k in red
            u[j,k] = 0.25 * (u[j-1:(M-1)/21,k] + u[j+1,k] + u[j,k-1] + u[j,k+1] -  h^(2.0) * F[j,k])
          end # red loop

          for k in black
            u[j,k] = 0.25 * (u[j-1,k] + u[j+1,k] + u[j,k-1] + u[j,k+1] -  h^(2.0) * F[j,k])
          end # black loop
        end #outer

      else
        for j in 2:M-1, k in 2:M-1
          u[j,k] = 0.25 * (u[j-1,k] + u[j+1,k] + u[j,k-1] + u[j,k+1] - h^(2.0)  * F[j,k])
        end
      end # if/else

      iterdiff = vecnorm(u - V,1)
      if iterdiff < tol*vecnorm(u)
        #println("Tolerance reached after $iter iterations")
        return u, iter
      end
      V = copy(u)
  end # iteration loop
  #println("Tolerance not reached")

  A = laplace_mat2D(h)

  res = F - A*u

  return u, maxiter

  #####
end # function

##################################################################

####### Successive Over Relaxation Iteration

function SOR(h::Float64, maxiter::Int64, tol::Float64)
  # Calculate optimal F
  w = 2.0/(1+sin(h*pi))
  # Set up mesh and F (right hand side)
  mesh = [j for j in 0:h:1]
  M = length(mesh)
  F = zeros(M,M)

  for j in 1:M, k in 1:M
      F[j,k] = funcRHS(mesh[j],mesh[k])
  end

  # Set up  initial guess and boundary data, this is homogenous Dirichlet
  u = zeros(M,M)
  V = zeros(M,M)
  # Run the iterations
  for iter in 0:maxiter

    for j in 2:M-1
      for k in 2:M-1
        u[j,k] = w*( 0.25 * (u[j-1,k] + u[j+1,k] + u[j,k-1] + u[j,k+1] - h^(2.0) * F[j,k])) + (1-w)*u[j,k]
      end
    end

    iterdiff = vecnorm(u-V,1)
    if iterdiff < tol*vecnorm(u)
      #println("Tolerance reached after $iter iterations")
      return u, iter
    end

    V = copy(u)

  end
  #println("Tolerance not reached")
  return u, maxiter


  #####
end # function

###################################################
# Outputs n x n sparse matrix for discrete Laplacian in 1D
function laplace_mat1D(h)
  n = 1/h - 1
  n = convert(Int,n)# spacing
  ############ Set up A
  dl = ones(n-1)
  d = -2*ones(n)
  du = ones(n-1)
  A = full((1/h^2)*Tridiagonal(dl, d, du))
  return A
end

# Outputs n x n sparse matrix for discrete Laplacian in 2D
function laplace_mat2D(h)
  n = 1/h - 1
  n = convert(Int,n)# spacing
  ############ Get 1D laplace
  Lx = laplace_mat1D(h)
  Ly = laplace_mat1D(h)
  ############ 1D identities
  Ix = eye(n)
  Iy = eye(n)
  ############ Set up A
  A = kron(Iy, Lx) + kron(Ly,Ix)
  return A
end

function laplaceDir2d(h)
  A = laplace_mat2D(h)
  n = 1/h - 1
  n = convert(Int,n)# spacing
  grid = collect(linspace(0+1/(n+1),1-1/(n+1),n))
  M = length(grid)
  F = zeros(M,M)
  for j in 1:M, k in 1:M
      F[j,k] = funcRHS(grid[j],grid[k])
  end
  f = copy(vec(F))
  U = pinv(A)*f
end


function laplace_mat3d(h)
  n = 1/h - 1
  n = convert(Int,n)# spacing
  Lx = laplace_mat1D(h)
  Ly = laplace_mat1D(h)
  Lz = laplace_mat1D(h)

  Ix = eye(n)
  Iy = eye(n)
  Iz = eye(n)

  A = kron(Lx, kron(Iy,Iz)) + kron(Ix, kron(Ly,Iz)) + kron(Ix, kron(Iy, Lz))
  return A
end

function laplaceDir3d(h)
  A = laplace_mat2D(h)
  n = 1/h - 1
  n = convert(Int,n)# spacing
  grid = collect(linspace(0+1/(n+1),1-1/(n+1),n))
  M = length(grid)
  F = zeros(M,M,M)
  for j in 1:M, k in 1:M, l in 1:M
      F[j,k] = -exp(-(grid[j] - 0.25)^2 - (grid[k] - 0.6)^2 - (grid[l] - 0.5)^2)
  end
  f = copy(vec(F))
  U = pinv(A)*f
end
##############################################################################
