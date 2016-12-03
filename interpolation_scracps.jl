# This is part of the working code for some other bilinear interpolation scheme

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
