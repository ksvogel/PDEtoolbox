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
