# Counterintuitive behavior in julia

# WHAT THE HELL, VARIABLE SCOPE

function parentalunit(x::Int64)
    println("we passed x = $x")
    y = x + 7
    z = evilchild(y)
    println("Wait, x is now $x")
    println("Wait, y is now $y")
end

function evilchild(y::Int64)
    println("we passed y = $y")
    x = 2
    println("we now have x = $x")
    z = x + y
    y = 2
    return z
end



function recursioncheck(x)
  x -= 1
  y = x
  println(y)
  if x != 0
    recursioncheck(x)
  else
    return x
  end

end


b = [1 2; 3 4]
d = Vector{Int64}(3)

for i in 1:3
  a = copy(b)
  println(a)
  c = a[1,1]
  d[i] = c

end

A = zeros(9,9)

for i in 1:9, j in 1:9
  A[i,j] = i+j

end
