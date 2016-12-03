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
  y = x
  x -= 1
  println(y-x)
  if x > 0
    recursioncheck(x)
  else
    return y
  end

end
