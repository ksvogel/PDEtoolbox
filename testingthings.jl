# testing ideas script

B = ones(10,10)
M = 10
rb = 2:M-1
for j in 2:M-1
red = filter(k -> iseven(k+j), 2:M-1)
println(red)
end
