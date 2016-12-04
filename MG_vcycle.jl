# Kaela S. Vogel
#=
This is an implementation of multigrid methods to solve PDEs, specific type is V-cyle.
=#

# load PDEtool, which contains all of these functions
# remember to use module name with module.function notation

######################## main function to use multigrid solver
#=
turtles       keeps track of which grid level we are on
lagturtle     estimate for our function from previous level
F             current R.H.S.
v             v[1] is our pre-smoothing iterations, v[2] is post-smoothing
=#

function multigrid(u0::Array{Float64,2}, F::Array{Float64,2}, maxiter::Int64, tol::Float64, turtles::Int64, v::Array{Int64,2})

  # Set up
  uturtle = copy(u0)
  uprevious = copy(u0) # Varible to keep track of prev sol
  h = 2.0^(-turtles) # Get the grid spacing

  # Run the larger loop
  for iter in 1:maxiter

    # Run the V-cycle
    uturtle = vcycle(uturtle, F, turtles, v)

    # Check for convergence
    if vecnorm(uturtle-uprevious) < tol*vecnorm(uturtle)
      println("V-cycle - tolerance reached after $iter iterations")
      residual = rescalc(uturtle, h, F)
      res = vecnorm(residual,2)
      return uturtle, residual, iter
    end

    uprevious = copy(uturtle)
  end #main for loop

  residual = rescalc(uturtle, h, F)

  return uturtle, residual, maxiter

end #function

#################################################################################

######################## recursive function to perform V-cycle
#=
turtles       keeps track of which grid level we are on
lagturtle     estimate for our function from previous level
F             current R.H.S.
v             v[1] is our pre-smoothing iterations, v[2] is post-smoothing
=#

function vcycle(lagturtle::Array{Float64,2}, F::Array{Float64,2}, turtles::Int64, v::Array{Int64,2})

  h = 2.0^(-turtles) # Get the grid spacing

  if turtles > 1

    # Relax v[1] times
    uturtle, res, iterz = gauss_sidel(lagturtle, h, v[1], 10.0^(-20), 3, F)

    # Restrict residual to coarser grid
    rescoarse = squish(res)

    # Initial guess for the error we will solve for at each level
    errturtle = zeros(size(rescoarse))

    # SEND THE ERROR TURTLES DOWN, i.e., recursively call vcycle
    # This solves for the error A e = -r
    errturtle = vcycle(errturtle, - rescoarse, turtles-1, v)

    # Interpolate error from previous level to a finer grid
    errturtle = foomp(errturtle)

    # Correct current level
    uturtle = uturtle - errturtle

    # Relax v[2] times
    uturtle, res, iterz = gauss_sidel(uturtle, h, v[2], 10.0^(-20), 3, F)

    return uturtle# SEND THE TURTLES BACK UP

  else # When we reach the 3x3 grid with only 1 interior point
    # Just copy the turtle from the previous level
    uturtle = lagturtle

    # Direct solve for the single point
    uturtle[2,2] = -h^2/4 * F[2,2]

    # This ends the turtles going all the way down
    # function will exit and send the turtles all the way up
    return uturtle

  end #conditional statement


end #function
########################################################################

######################## function to apply restriction transformation
#=
 Maps a grid with spacing h to a grid with spacing 2h
In: M x M 2D array of function values
Out: (M+1)/2 x (M+1)/2 array of function values
=#
function squish(vf::Array{Float64,2})

  M = size(vf,1) # This includes the boundary
  mm = convert(Int,(M+1)/2) # restricted array size including boundary
  n=mm
  vS = zeros(mm, mm)

  for j in 4:2:M-1, k in 4:2:M-1
    l = convert(Int,j/2)
    p = convert(Int, k/2)
    vS[l,p] = 1/16 * (vf[j-1,k-1] + vf[j-1,k+1] + vf[j+1,k-1] + vf[j+1,k+1] + 2*(vf[j,k+1] + vf[j-1,k] + vf[j, k-1] + vf[j+1,k]) + 4*vf[j,k])
  end


  for j in 2:mm-1, k in 2:mm-1
          p = convert(Int, 2*j - 1);
          q = convert(Int, 2*k - 1);
          vS[i,j] = (1/16)*( vf[p-1,q-1] + 2*vf[p-1,q] + vf[p-1,q+1] + 2*vf[p,q-1] + 4*vf[p,q] + 2*vf[p,q+1] + vf[p+1,q-1] + 2*vf[p+1,q] + vf[p+1,q+1] )
  end

  return vS

end #function

######################## function to apply interpolation transformation
#=
 Maps a grid with spacing h to a grid with spacing 2h
In: (M+1)/2 x (M+1)/2 array of function values
Out: M x M 2D array of function values
=#
function foomp(vs::Array{Float64,2})
  mm = size(vs,1)
  M = convert(Int,2*mm-1) # doing inner points and then adding boundary
  vF = zeros(M, M)

# Map each point of smaller matrix to a 9x9 matrix in the larger ones
  Ires = 0.25*[1 2 1; 2 4 2; 1 2 1]

  for j in 2:mm - 1, k in 2:mm-1
    p = 2*j-2
    q = 2*k-2
    vF[p:p+2,q:q+2] = vF[p:p+2,q:q+2] + vs[j,k]*Ires;
  end

  return vF

end #function

##################### Function to move recursively down the grid levels

function turtlesallthewaydown(Ua, turtles, Fs, Rs,  hs, s1, s2)

  turtles -= 1
  #println("We have $turtles turtles to go, h is $(hs[turtles])")
  local f_prev = copy(Rs[turtles+1])
  local f_2hd = squish(f_prev) # restrict the residual to coarser grid
  Fs[turtles] = copy(f_2hd)

  # Relax s1 times
  e_2hd, rz, iterz = gauss_sidel(zeros(size(f_2hd)), hs[turtles], s1, 10.0^(-20), 1, -Fs[turtles])
  Ua[turtles] = copy(e_2hd)
  Rs[turtles] = copy(rz)
  #println(vecnorm(rz,2))
  if turtles > 2 # Call this function again, move to next grid
    turtlesallthewaydown(Ua, turtles, Fs, Rs,  hs, s1, s2)
  else # Direct solve rather than relax on coarsest grid
    #CURRENTLY STILL ITERATIVE
    turtles -= 1
    #println("This is the final turtle")
    f_prevfinal = copy(Rs[turtles+1])
    f_4h = squish(f_prevfinal) # restrict the residual to coarser grid
    E = zeros(3,3)
    hfinal = copy(hs[turtles])
    E[2,2] = -hfinal/4 * f_4h[2,2]
    Ua[turtles] = copy(E)

    return Ua, Fs, Rs
  end #conditional

end #function

##################### Function to move recursively up the grid levels

function turtlesallthewayup(Ua, Rs, turtles, Fs, hs, s1, s2)

  turtles += 1
  #println("We are at turtle level $turtles")
  e_prev = copy(Ua[turtles-1])
  corr = foomp(e_prev) # interpolate error from previous grid
  Uacorrect = copy(Ua[turtles]) - corr
  cnorm = vecnorm(corr, 2)
  #println("corr is $cnorm")
  # Relax s2 times
  e_hup, r_hup, iterz = gauss_sidel(Uacorrect, hs[turtles], s2, 10.0^(-50), 1, -Fs[turtles])
  Ua[turtles] = copy(e_hup)
  Rs[turtles] = copy(r_hup)
  #println(vecnorm(r_hup,2))


  if turtles <  size(hs,1)# Call this function again, move up to next grid
    turtlesallthewayup(Ua, Rs, turtles, Fs, hs, s1, s2)
  else
    return Ua, Fs, Rs
  end #conditional

end #function

##################### Function to run the multigrid V-cycle
#=
s1   number of tiems to smooth on restricion steps
s2    number of times to smooth on prolongation steps
u   initial guess for solution
h   fine grid spacing
turtles  keeps track of which grid level we are on
F   RHS of original PDE descritized and evaluated on fine grid

  =#

function MG_vcycle(h, F, u0, turtles, s1, s2, tol)

  hs = [(2.0^j)*h for j in turtles-1:-1:0] # Vector to store all the grid spacing
  Rs = [u0 for k in 1:turtles] # Vector to store residuals at each level
  turtletrack = copy(turtles)

  M = size(F,2)
  # Array of arrays that will hold the approx sol. at each level
  Ua = [u0 for k in 1:turtles]

  # Array of arrays that will hold the RHS at each level
  Fs = [F for k in 1:turtles]

  maxiter = 1000
  V = copy(Ua[turtles])
  for iter in 1:maxiter
    # Reset level counter
    turtles = copy(turtletrack)

    # Relax s1 times
    uoutz, r_h, iterz = gauss_sidel(Ua[turtles], hs[end], s1, 10.0^(-20), 3, F)
    #println(vecnorm(r_h,2))

    # Store the new estimation of u, uout
    Ua[turtles] = copy(uoutz)
    Rs[turtles] = copy(r_h)

    # Start the process of restriction and working on coarser grids
    Ua, Fs, Rs = turtlesallthewaydown(Ua, turtles, Fs, Rs,  hs, s1, s2)

    # AND NOW TO START CORRECTING INITIAL GUESSES
    Ua, Fs, Rs = turtlesallthewayup(Ua, Rs, 1, Fs, hs, s1, s2)

    # Check if relative residual error is less than tolerance
    v = copy(Ua[turtles])
    residual = rescalc(v, hs[end], F)
    res = vecnorm(residual,2)
    #println("residual after cycle $res")

      if vecnorm(Ua[turtles]-V,1) < tol*vecnorm(Ua[turtles],1)
        println("V-cycle - tolerance reached after $iter iterations")
        return Ua, residual, iter
      end
          V = copy(Ua[turtles])
    end
    # Check if relative residual error is less than tolerance

    residual = rescalc(V, hs[end], F)
  println("beep")
  return Ua, residual, iter

end #function
