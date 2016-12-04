# Kaela S. Vogel
# This is a module of some of the different iterative methods to solve the Poisson Equation
# To call it you need "using(poisson_iterative.jl)"
# h - mesh spacing
#





####### Function to set up mesh and RHS for a given spacing h
#=
funcRHS anonymous function, e.g. (x,y) -> x^2+y

=#

function h_space(funcRHS::Function, h::Float64)
  mesh = [j for j in 0:h:1]
  M = length(mesh)
  F = zeros(M,M)

  for j in 1:M, k in 1:M
      F[j,k] = funcRHS(mesh[j],mesh[k])
  end

  return mesh,F

end #function

###### Function to calculate the residual
function rescalc(u, h, F)
  M = size(u,1)
  lap2D = zeros(M, M)

  for j in 2:M-1, k in 2:M-1
    lap2D[j,k] = h^(-2) * ( u[j-1,k] + u[j+1,k] - 4* u[j,k] + u[j,k-1] + u[j,k+1] )
  end

  return (F-lap2D)

end



####### Jacobi Iteration
# use differance between successive iterations for tolerance
function jacobi_iter(h, maxiter, tol, F)
  # Set up mesh and F (right hand side)
M = size(F,1)
  # Set up  initial guess and b oundary data, this is homogenous Dirichlet
  u = zeros(size(F))
  u[1,:] = zeros(1,M)
  u[M,:] = zeros(1,M)
  V = zeros(size(u))
  for iter in 0:maxiter

    for j in 2:M-1, k in 2:M-1
        V[j,k] = 0.25 * (u[j-1,k] + u[j+1,k] + u[j,k-1] + u[j,k+1] - h^(2.0) * F[j,k])
    end

    iterdiff = vecnorm(u - V,1)
    if iterdiff < tol*vecnorm(V)
        println("Tolerance reac hed after $iter iterations")
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

function gauss_sidel(u::Array{Float64,2}, h::Float64, maxiter::Int64, tol::Float64, rb, F::Array{Float64,2})

  local M = size(u, 2)
  # Run the iterations
  V = zeros(size(u))

  for iter in 0:maxiter


      for j in 2:M-1
        rs = filter(k -> iseven(k+j), 2:M-1)
        red = convert(Array{Int64,1}, rs)
        bs = filter(k -> isodd(k+j), 2:M-1)
        black = convert(Array{Int64,1}, bs)
          for k in red
            u[j,k] = 0.25 * (u[j-1,k] + u[j+1,k] + u[j,k-1] + u[j,k+1] -  h^(2.0) * F[j,k])
          end # red loop

          for k in black
            u[j,k] = 0.25 * (u[j-1,k] + u[j+1,k] + u[j,k-1] + u[j,k+1] -  h^(2.0) * F[j,k])
          end # black loop
        end #outer
        #=
        if rb == 'rb'
      else
        for j in 2:M-1, k in 2:M-1
          u[j,k] = 0.25 * (u[j-1,k] + u[j+1,k] + u[j,k-1] + u[j,k+1] - h^(2.0)  * F[j,k])
        end
      end # if/else
=#
      iterdiff = vecnorm(u-V,1)
      residual = rescalc(u, h, F)

      #if vecnorm(residual,1) < tol*vecnorm(F,1)
      if iterdiff < tol*vecnorm(u,1)
        println("GS RB Tolerance reached after $iter iterations")
        #residual = rescalc(u)
        return u, residual, iter
      end

      V = copy(u)
  end # iteration loop
  #println("Tolerance not reached")
  residual = rescalc(u, h, F)
  return u, residual, maxiter

  #####
end # function

##################################################################

####### Successive Over Relaxation Iteration

function SOR(h::Float64, maxiter::Int64, tol::Float64, F::Array{Float64,2})
  # Calculate optimal F
  w = 2.0/(1+sin(h*pi))

  M = size(F,2)


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

    iterdiff = vecnorm(u-V)
    residual = rescalc(u, h, F)
    if iterdiff < tol*vecnorm(u)
      println("SOR Tolerance reached after $iter iterations")
      #residual = rescalc(u)
      return u, residual, iter
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
