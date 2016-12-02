# Kaela S. Vogel
#=
This is an implementation of multigrid methods to solve PDEs, specific type is V cyles\
=#

#load PDEtool, remember to use module name with module.function notation


######################## function to apply restriction transformation
#=
 Maps a grid with spacing h to a grid with spacing 2h
In: M x M 2D array of function values
Out: M/2 + 1 x M/2 + 1 array of function values
*     M IS AN ODD NUMBER IN GENERAL
ADD A WAY TO CATCH UNEVEN M-1
=#
function squish(v)

  M = size(v,1) # This includes the boundary
  mm = convert(Int,(M+1)/2) # restricted array size including boundary

  vS = zeros(mm, mm)

  for j in 4:2:M-1, k in 4:2:M-1
    l = convert(Int,j/2)
    p = convert(Int, k/2)
    vS[l,p] = 1/16 * (v[j-1,k-1] + v[j-1,k+1] + v[j+1,k-1] + v[j+1,k+1] + 2*(v[j,k+1] + v[j-1,k] + v[j, k-1] + v[j+1,k]) + 4*v[j,k])
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
function foomp(v)
  mm = size(v,1)
  M = convert(Int,2*mm-1) # doing inner points and then adding boundary
  vF = ones(M, M)

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

  return vF

end #function

##################### Function to move recursively down the grid levels

function turtlesallthewaydown(Ua, turtles, Fs, r_h,  h, s1, s2)

  turtles -= 1
  println("We have $turtles turtles to go")
  hr = 2*h # coarse the grid
  f_2h = squish(r_h) # restrict the residual to coarser grid

  # Relax s1 times
  e_2h, r_2h, maxiter = gauss_sidel(zeros(size(f_2h)), hr, s1, 1.0^(-6), 1, f_2h)
  Ua[turtles] = e_2h
  Fs[turtles] = f_2h

  if turtles > 2 # Call this function again, move to next grid
    turtlesallthewaydown(Ua, turtles, Fs, r_2h,  hr, s1, s2)
  else # Direct solve rather than relax on coarsest grid
    #CURRENTLY STILL ITERATIVE
    turtles -= 1
    println("This is the final turtle")
    he = 2*hr # coarse the grid
    f_4h = squish(r_2h) # restrict the residual to coarser grid
    E, R, maxiter = gauss_sidel(zeros(size(f_4h)), he, 1000, 1.0^(-6), 1, f_4h)
    Ua[turtles] = E
    Fs[turtles] = f_2h
    turtlesallthewayup(Ua, turtles, Fs, r_h,  h, s1, s2)
  end #conditional

end #function

##################### Function to move recursively up the grid levels

function turtlesallthewayup(Ua, turtles, Fs, r_h,  h, s1, s2)

  turtles += 1
  println("We are at turtle level $turtles")
  hr = 2*h # fine the grid, it should pay
  f_2h = squish(r_h) # restrict the residual to coarser grid

  # Relax s1 times
  e_2h, r_2h, maxiter = gauss_sidel(zeros(size(f_2h)), hr, s1, 1.0^(-6), 1, f_2h)
  Ua[turtles] = e_2h
  Fs[turtles] = f_2h

  if turtles > 2 # Call this function again, move to next grid
    turtlesallthewaydown(Ua, turtles, Fs, r_2h,  hr, s1, s2)
  else # Actuall solve rather than relax on coarsest grid
    turtles -= 1
    println("This is the final turtle")
    he = 2*hr # coarse the grid
    f_4h = squish(r_2h) # restrict the residual to coarser grid
    E, R, maxiter = gauss_sidel(zeros(size(f_4h)), he, 1000, 1.0^(-6), 1, f_4h)
    Ua[turtles] = E
    Fs[turtles] = f_2h
    turtlesallthewayup(v, turtles, residual, h)
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

function MG_vcycle(h, F, u, turtles, s1, s2)

  # Array of arrays that will hold the approx sol. at each level
  Ua = [u for k in 1:turtles]

  # Array of arrays that will hold the RHS at each level
  Fs = [F for k in 1:turtles]

  # Relax s1 times
  uout, r_h, maxiter = gauss_sidel(Ua[turtles], h, s1, 1.0^(-6), 1, F)

  # Store the new estimation of u, uout
  Ua[turtles] = uout

  # Start the process of restriction and working on coarser grids
  turtlesallthewaydown(Ua, turtles, Fs, r_h,  h, s1, s2)

end #function
