# Kaela S. Vogel
#=
This is an implementation of multigrid methods to solve PDEs, specific type is V cyles\
=#

######################## function to apply restriction transformation
#=
 Maps a grid with spacing h to a grid with spacing 2h
In: M x M 2D array of function values
Out: M/2 + 1 x M/2 + 1 array of function values

ADD A WAY TO CATC UNEVEN M-1
=#
function squish(v::Array{Int64})

  vS = zeros(M/2 -1, M/2 - 1)

  for j in 2:2:M-1, k in 2:2:M-1
    l = convert(Int,j/2)
    p = convert(Int, k/2)
    vS[l,p] = 1/16 * (v[j-1,k-1] + v[j-1,k+1] + v[j+1,k-1] + v[j+1,k+1] + 2*(v[j,k+1] + v[j-1,k] + v[j, k-1] + v[j+1,k]) + 4*v[j,k])
  end

  vSbdry = zeros(M/2 + 1, M/2 + 1)
  vSbdry[2:M/2-1, 2:M/2-1] = vS

  return vSbdry

end #function
