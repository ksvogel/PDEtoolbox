#=

Kaela S. Vogel 2016
Functions for performing the preconditioned conjugate gradient methods

Preconditioner options, pass a string literal

"np"      No preconditioning
"SSOR"    SSOR preconditioning
"MG"      Multi-grid preconditioning

WHAT IS INITIAL GUESS?????????????
=#

function PCG(u0::Array{Float64,2}, F::Array{Float64,2}, h::Float64, tol::Float64, precon)

  # Inititialize things
  r = rescalc(u0, h, F)
  p = r
  u = u0


  # Precondition
  z = preconditioner(r, r, h, precon)
  ztrnext = dot(vec(z), vec(r))

  for k in 1:length(z)

    w = applylap(p, h, F)
    ztr = copy(ztrnext)
    alpha = ztr/dot(vec(p),vec(w))
    u = u + alpha * p
    r = r - alpha * w

    # Stopping condition
    if vecnorm(r) < tol*vecnorm(F,1)
      println("PCG tolerance reached after $k iterations")
      return u, r, k
    end

    z = preconditioner(z, r, h, precon)
    ztrnext = dot(vec(z), vec(r))
    beta = ztrnext/ztr

    if precon == "np"
      p = r
    else
      p = z + beta * p
    end


  end #main loop
  k = length(z)
  return u, r, k

end #function

#################################################################################

############################## function Preconditioner

#=

Preconditioner options, pass a string literal

"np"      No preconditioning
"SSOR"    SSOR preconditioning
"MG"      Multi-grid preconditioning

=#

function preconditioner(z0, r0, h, precon)
  # Preconditioner options
  M = size(z0,2)

  if precon == "np"
    z0 = r0

  elseif precon == "SSOR"
    # Calculate optimal omega
    w = 2.0/(1+sin(h*pi))

    # 1 Sweep of SSOR
    for j in 2:M-1
      for k in 2:M-1
        # Forward Sweep SOR
        z0[j,k] = w*( 0.25 * (z0[j-1,k] + z0[j+1,k] + z0[j,k-1] + z0[j,k+1] - h^(2.0) * r0[j,k])) + (1-w)*z0[j,k]
        # Backward Sweep SOR
        z0[j,k] = w*( 0.25 * (z0[j-1,k] + z0[j+1,k] + z0[j,k-1] + z0[j,k+1] - h^(2.0) *r0[j,k])) + (1-w)*z0[j,k]
      end
    end

  elseif precon == "MG"
    turtles = Int(- log(h)/log(2)) # will need to convert to int
    v = [2 2]
    z0, residual, maxiter = multigrid(zeros(size(r0)), r0, 1, 10.0^(-4), turtles, v)

  else
    println("No precondition option specified, program exited")
  end # precondition conditional

  return z0

end
