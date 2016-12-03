# Kaela S. Vogel
#=
This is an implementation of multigrid methods to solve PDEs, specific type is V-cyle.
=#

# load PDEtool, which contains all of these functions
# remember to use module name with module.function notation


######################## function to apply restriction transformation
#=
 Maps a grid with spacing h to a grid with spacing 2h
In: M x M 2D array of function values
Out: M/2 + 1 x M/2 + 1 array of function values
*     M IS AN ODD NUMBER IN GENERAL
ADD A WAY TO CATCH UNEVEN M-1
=#
function squish(vf)

  local M = size(vf,1) # This includes the boundary
  local mm = convert(Int,(M+1)/2) # restricted array size including boundary

  vS = zeros(mm, mm)

  for j in 4:2:M-1, k in 4:2:M-1
    l = convert(Int,j/2)
    p = convert(Int, k/2)
    vS[l,p] = 1/16 * (vf[j-1,k-1] + vf[j-1,k+1] + vf[j+1,k-1] + vf[j+1,k+1] + 2*(vf[j,k+1] + vf[j-1,k] + vf[j, k-1] + vf[j+1,k]) + 4*vf[j,k])
  end

  return vS

end #function

######################## function to apply interpolation transformation
#=
 Maps a grid with spacing h to a grid with spacing 2h
In: M/2 + 1 x M/2 + 1 array of function values
Out: M x M 2D array of function values

ADD A WAY TO CATCH UNEVEN M-1
=#
function foomp(vs)
  local mm = size(vs,1)
  local M = convert(Int,2*mm-1) # doing inner points and then adding boundary
  vF = zeros(M, M)

# Map each point of smaller matrix to a 9x9 matrix in the larger ones
  Ires = 0.25*[1 2 1; 2 4 2; 1 2 1]

  for j in 2:mm - 1, k in 2:mm-1
    p = 2*j - 2
    q = 2*k - 2
    vF[p:p+2,q:q+2] = vF[p:p+2,q:q+2] + vs[j,k]*Ires;
  end
#=
# Map the points that line up in the course grid and fine grid
  for j in 0:mm-1, k in 0:mm-1
    vF[2*j+1,2*k+1] = v[j+1,k+1]
  end


# Interpolate the points directly right and below each corse grid point on fine grid
  for j in 1:M
    black = filter(k -> isodd(k+j), 1:M)
    for k in black
      if isodd(j)
        vF[j, k] = 1/2*(vF[j,k-1] + vF[j,k+1])
      else
        vF[j, k] = 1/2*(vF[j-1,k] + vF[j+1,k])
      end # conditional
    end # red loop
  end #j loop

# Interpolate the points on the oposing corner of each 4x4 unit of fine grid
  for j in 2:2:M-1, k in 2:2:M-1
    vF[j,k] = 1/4*(vF[j-1,k-1] + vF[j+1, k-1] + vF[j-1,k+1] + vF[j+1,k+1])
  end
=#
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
  e_2hd, rz, iterz = gauss_sidel(zeros(size(f_2hd)), hs[turtles], s1, 10.0^(-20), 1, Fs[turtles])
  Ua[turtles] = copy(e_2hd)
  Rs[turtles] = copy(rz)
  println(vecnorm(rz,2))
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
  println("We are at turtle level $turtles")
  e_prev = copy(Ua[turtles-1])
  corr = foomp(e_prev) # interpolate error from previous grid
  Uacorrect = copy(Ua[turtles]) - corr
  cnorm = vecnorm(corr, 2)
  println("corr is $cnorm")
  # Relax s2 times
  e_hup, r_hup, iterz = gauss_sidel(Uacorrect, hs[turtles], s2, 10.0^(-50), 1, Fs[turtles])
  Ua[turtles] = copy(e_hup)
  Rs[turtles] = copy(r_hup)
  println(vecnorm(r_hup,2))


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

  maxiter = 100

  for iter in 1:maxiter
    # Reset level counter
    turtles = copy(turtletrack)

    # Relax s1 times
    uoutz, r_h, iterz = gauss_sidel(Ua[turtles], hs[end], s1, 10.0^(-20), 3, F)
    println(vecnorm(r_h,2))

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
    println("residual after cycle $res")

      if vecnorm(residual,1) < tol*vecnorm(F,1)
        println("V-cycle - tolerance reached after $iter iterations")
        return Ua, residual
      end

    end
    # Check if relative residual error is less than tolerance
    v = copy(Ua[turtles])
    residual = rescalc(v, hs[end], F)
  println("beep")
  return Ua, residual

end #function
