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
  M = size(u0,2)
  meep = F[2:end-1,2:end-1]
  F = zeros(M,M)
  F[2:end-1,2:end-1] = meep
  r = rescalc(u0, h, F)
  u = copy(u0)
  # Precondition
  z = preconditioner(u0, r, h, precon)
  p = copy(z)

  for k in 1:length(z)

    w = applylap(p, h)
    alpha = dot(vec(z), vec(r))/dot(vec(p),vec(w))
    r0 = copy(r)
    r = r0 - alpha * w
    u = u + alpha * p
    # Stopping condition
    if vecnorm(r) < tol
      println("PCG tolerance reached after $k iterations")
      return u, r, k
    end
    z0 = copy(z)
    z = preconditioner(u0, r, h, precon)
    beta = dot(vec(z), vec(r))/dot(vec(z0), vec(r0))
    p = z + beta * p


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
    return z0

  elseif precon == "SSOR"
    z0 = zeros(size(r0))

    # Calculate optimal omega
    w = 2.0/(1+sin(h*pi))


    for j in 2:M-1, k in 2:M-1
        # Forward Sweep SOR
        z0[j,k] = (w/4) *(z0[j-1,k] + z0[j+1,k] +  z0[j,k-1] + z0[j,k+1] - h^(2.0) * r0[j,k]) + (1-w)*z0[j,k]
    end

    for j in M-1:-1:2, k in M-1:-1:2
        # Backward sweep of SOR
        z0[j,k] = (w/4) *(z0[j-1,k] + z0[j+1,k] +  z0[j,k-1] + z0[j,k+1] - h^(2.0) * r0[j,k]) + (1-w)*z0[j,k]
    end

    return z0

  elseif precon == "MG"
      z0 = zeros(size(r0))
      turtles = Int(- log(h)/log(2)) # will need to convert to int
      v = [1 1]
      z0 = vcycle(z0, r0, turtles, v)
      return z0

  else
    println("No precondition option specifiedprogram exited")
  end # precondition conditional
  println("I got here")
end #function
