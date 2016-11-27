# Kaela S. Vogel
# This is a module of some of the different iterative methods to solve the Poisson Equation
# To call it you need "using(poisson_iterative.jl)"
# h - mesh spacing
#

####### Set up our R.H.S. function
funcRHS(x,y) = -exp(-(x - 0.25)^2 - (y - 0.6)^2)



####### Jacobi Iteration
# use differance between successive iterations for tolerance
function jacobi_iter(h, maxiter, tol)
  # Set up mesh and F (right hand side)
  mesh = [j for j in 0:h:1]
  M = length(mesh)
  F = zeros(M,M)

  for j in 1:M, k in 1:M
      F[j,k] = funcRHS(mesh[j],mesh[k])
  end


  # Set up  initial guess and boundary data, this is homogenous Dirichlet
  u = zeros(size(F))
  u[1,:] = zeros(1,M)
  u[M,:] = zeros(1,M)
  u[:,1] = zeros(M,1)
  u[:,M] = zeros(M,1)
  V = zeros(size(u))
  # Run the iterations
  for iter in 0:maxiter
    for j in 2:M-1, k in 2:M-1
        V[j,k] = 0.25 * (u[j-1,k] + u[j+1,k] + u[j,k-1] + u[j,k+1] - h^(2.0) * F[j,k])
    end

    iterdiff = vecnorm(u - V,1)
    if iterdiff < tol*vecnorm(V)
        #println("Tolerance reached after $iter iterations")
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

function gauss_sidel(h, maxiter, tol)
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
  for iter in 0:maxiter

    for j in 2:M-1, k in 2:M-1
        u[j,k] = 0.25 * (u[j-1,k] + u[j+1,k] + u[j,k-1] + u[j,k+1] - h^(2.0) * F[j,k])
    end
    iterdiff = vecnorm(u - V,1)
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

##################################################################

####### Successive Over Relaxation Iteration

function SOR(h, maxiter, tol)
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
    iterdiff = vecnorm(u - V,1)
    if iterdiff < tol*vecnorm(u)
      #println("Tolerance reached after $iter iterations")
      return u, iter
    end
    V = copy(u)

  end
  #println("Tolerance not reached")
  return u, maxiter, w


  #####
end # function
#=
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
#####################################################################################
# code to run for problem 1 HW 3

hs = [2.0^(-5) 2.0^(-6) 2.0^(-7)]
maxiter = 1000
tol = 10.0^(-2)



for j in 1:3
  h = hs[j]
  u1, iter1 = jacobi_iter(hs[j], maxiter, tol)
  u2, iter2 = gauss_sidel(hs[j], maxiter, tol)
  u3, iter3, w = SOR(hs[j], maxiter, tol)
  println("$h $iter1 $iter2 $iter3 $w")
end
=#
#########################################################################################
# 2D direct solve vs SOR
maxiter = 1000
tol = 10.0^(-2)
hs = [2.0^(-7) 2.0^(-8) 2.0^(-9) 2.0^(-10)]
hs = [2.0^(-5) 2.0^(-6) 2.0^(-7)]

# Set up F for Direct solve matrices

for j in 1:3
@time u3, iter3, w = SOR(hs[j], maxiter, tol)
#@time laplaceDir2d(hs[j])
println("Iteration $j Done, $iter3")
end

# 3D
hs = [2.0^(-5) 2.0^(-6) 2.0^(-7)]

#=
for j in 1:3
#@time u3, iter3, w = SOR(hs[j], maxiter, tol)
@time laplaceDir3d(hs[j])
println("Iteration $j Done")
end
  # @time FUNCTION println("Done $x")
=#
