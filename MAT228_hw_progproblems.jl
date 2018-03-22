#=
Kaela S. Vogel
2018-01-01
Homeworks for MAT 228

=#
include("PDEtool.jl")
using PDEtool
using Plots
plotlyjs()


##############################################################

#= Homework 4, problem 1
Solve the advection equation withperiodic boundry conditions using Lax-Wendroff and upwinding.
=#

a =  1.0
h = 1/2^4 # Grid spacing
k = 0.9*a*h # Tine stepping
s = Bool(1) # Switch variable to give basic mesh
funcRHS( x, t) = 0

u_xy0(x, y)= exp.(-10((x-.3).^2 + (y-.4).^2))

##############################################################

#= Homework 3, problem 1
Solve the heat equation with Neumann boundry conditions using the Peaceman-Rachford ADI scheme on a cell-centered grid.
=#


a =  0.1
h = 1/3^4 # Grid spacing
k = 1/3^4# Tine stepping
s = Bool(0) # Switch variable to give cell centered grid
funcRHS( x, t) = 0
u_0t(t) = 0 # boundary condition at x = 0
u_1t(t) = 0 # boundary condition at x = 1
u_xy0(x, y)= exp.(-10((x-.3).^2 + (y-.4).^2))
#v = prADI_heat2D(h, k, funcRHS,  u_0t, u_1t, u_xy0, a, s)


# The grid refinement study using the infinity-norm
hstart =  2
hend = 7
spacing = [1/3^j for j in hstart:1:hend]
v_unrefined = Vector{Any}(length(spacing)) # This will be an array of matrices of different dimsfuncRHS

for j in 1:length(spacing) # The indexing assumes h and k are equal
    v = PDEtool.prADI_heat2D(spacing[j], spacing[j], funcRHS,  u_0t, u_1t, u_xy0, a, s)
    v_unrefined[j] = v
end

ratios, v_restricteds = refinement_ratios(v_unrefined, spacing, spacing, s)
writecsv("hw_3_refinement.txt", ratios)


# Problem 1d, showing discrete conservation law
mass_beginning = zeros(length(v_unrefined))
mass_end = zeros(length(v_unrefined))
mass_diff = zeros(length(v_unrefined))
funcRHS
for l in 1:length(spacing)funcRHS
    hmesh, kmesh, F, M, T = meshmaker(funcRHS, spacing[l], k, s)
    m = zeros(M,M)
    for j in 1:M, i in 1:M
        m[j,i] = u_xy0(hmesh[j], hmesh[i])
    endfuncRHS
    mass_beginning[l] = sum(m)
end


for j in 1:length(v_unrefined)funcRHS
    mass[j] = sum(v_unrefined[j])
    mass_diff[j] = mass[j] - mass_beginning[j]
end
conserved = [mass_beginning mass mass_diff]
writecsv("hw_3_conserved.txt", conserved)


#= Homework 3, problem 2
Solve the The FitzHugh-Nagumo equations with Neumann boundry conditions using the Peaceman-Rachford ADI scheme on a cell-centered grid and Forward Euler
=#


h = 1/2^8  # Grid spacing
k = 1.0 # Time stepping
T = 600
s = Bool(0) # Switch variable to give cell centered grid
v_xy0(x, y) = exp.(-100(x.^2 + y.^2))
w_xy0(x, y) = 0.0
vsol1 = FHN2D(h, k, v_xy0, w_xy0, T, s)


# Part c
h = 1/2^8 # Grid spacing
k = 1.0 # Time stepping
T = 600
s = Bool(0) # Switch variable to give cell centered grid
v_xy0(x, y) = 1 - 2.*x
w_xy0(x, y) = 0.05.*y
vsol2 = FHN2D(h, k, v_xy0, w_xy0, T, s)
vsol2 = vsol


fig1 = surface( vsol1[1] )
fig2 = surface( vsol1[2] )
fig3 = surface( vsol1[3] )
fig4 = surface( vsol1[4] )
fig5 = surface( vsol1[5] )
fig6 = surface( vsol1[6] )
fig7 = surface( vsol1[7] )
fig8 = surface( vsol1[8] )
fig9 = surface( vsol1[10])
#data = [ fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9 ]
plot(fig1, fig2, fig1, fig4, fig5, fig6, fig7, fig8, fig9)

fig1 = surface( vsol1[1])
fig2 = surface( vsol1[2])
fig3 = surface( vsol1[3])
fig4 = surface( vsol1[4])
fig5 = surface( vsol1[5])
fig6 = surface( vsol1[6])
fig7 = surface( vsol1[7])
fig8 = surface( vsol1[8])
fig9 = surface( vsol1[9])
fig1 = surface( vsol1[10])
fig1 = surface( vsol1[11])
fig2 = surface( vsol1[12])
fig3 = surface( vsol1[13])
fig4 = surface( vsol1[14])
fig5 = surface( vsol1[15])
fig6 = surface( vsol1[16])
fig7 = surface( vsol1[17])
fig8 = surface( vsol1[18])
fig9 = surface( vsol1[19])
#data = [ fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9 ]
plot(fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9)

##############################################################

#= Homework 2, problem 2
This is the setup and the function call to solve the heat equation
size(v_unrefined[3],2)),1:-1/size(v_unrefined[3]
u_t = \Delta u on  \Omega = (0,1) X (0,1
u(0,0,t) = 0  boundary
u(x,y,0) = exp[-100((x - 0.3)^2 + (y - 0.4)^2)]

with Crank-Nicolson 2D and Forwara =  0.1
h = 1/2^5 # Grid spacing
k = h^2/(4*a)# Tine stepping
funcRHS( x, t) = 0
u_0t(t) = 0 # boundary condition at x = 0
u_1t(t) = 0 # boundary condition at x = 1
u_xy0(x, y)= exp.(-100((x-.5).^2 + (y-.5).^2))d Euler solvers. Both are timed for varying grid sizes.

=#


#cd("/home/kaela/Documents/GithubRepositories/PDE_solvers")

a =  0.1
h = 1/2^5 # Grid spacing
k = h^2/(4*a)# Tine stepping
funcRHS( x, t) = 0
u_0t(t) = 0 # boundary condition at x = 0
u_1t(t) = 0 # boundary condition at x = 1
u_xy0(x, y)= exp.(-100((x-.5).^2 + (y-.5).^2))
v = PDEtool.FE_heat2D(h, k, funcRHS,  u_0t, u_1t, u_xy0, a)

# The indexing assumes h and k are equal

a = 1.0
h = 1/2^6 # Grid spacing
k = 0.01
funcRHS( x, t) = 0
u_0t(t) = 0 # boundary condition at x = 0
u_1t(t) = 0 # boundary condition at x = 1
u_xy0(x, y)= exp.(-100((x-.5).^2 + (y-.5).^2)) #-exp.(-.5^2) #exp(-100((x -
v = PDEtool.CN_heat2D(h, k, funcRHS,  u_0t, u_1t, u_xy0, a)
surface(v)

hstart =  3
hend = 8
FE_timeeelapsed = zeros(1,hend-hstart+1)
CN_timeeelapsed = zeros(1,hend-hstart+1)
spacing = [1/2^j for j in hstart:1:hend]


for j in 1:length(spacing)

    FE_timeeelapsed[j] = @elapsed PDEtool.FE_heat2D(spacing[j], spacing[j]^2/(4*a) , funcRHS,  u_0t, u_1t, u_xy0, a)

    #CN_timeeelapsed[j] = @elapsed PDEtool.CN_heat2D(spacing[j], 0.01, funcRHS,  u_0t, u_1t, u_xy0, a)

end
ratio_CN = zeros(5,1)
ratio_FE = zeros(5,1)
for i in 1:5
    #ratio_CN[i] = CN_timeeelapsed[i+1]/CN_timeeelapsed[i]
    ratio_FE[i] = FE_timeeelapsed[i+1]/FE_timeeelapsed[i]

end
##############################################################

#= Homework 1, problem 2
This is the setup and the function call to solve the heat equation

u_t = 0.01u_{xx} + 1 - \exp(-t), 0 < x < 1
u(0,t) = 0,
u(1,t) = 0
u(x,0) = 0

with Crank-Nicolson and then perform a refinement study to show it is second order accurate in space and time.

=#

h = 1/2^5 # Grid spacing
k = .1 # Tine stepping
funcRHS( x, t) = 1 - exp.( - t)
u_0t(t) = 0 # boundary condition at x = 0
u_1t(t) = 0 # boundary condition at x = 1
u_x0(x)= 0
#_x0(x) = exp.(-(x-.5).^2)-exp.(-.5^2)
#u_x0(x) = x < 0.5 ? x : 1-x # initial condtion at t = 0
a =  0.1
v = PDEtool.CN_heat(h, k, funcRHS,  u_0t, u_1t, u_x0, a)size(v_unrefined[3],2)),1:-1/size(v_unrefined[3]
heatmap(v)

# The grid refinement study using the infinity-norm
hstart =  3
hend = 10
spacing = [1/2^j for j in hstart:1:hend]
v_unrefined = Vector{Any}(length(spacing)) # This will be an array of matrices of different dims

for j in 1:length(spacing) # The indexing assumes h and k are equal
    v = PDEtool.CN_heat(spacing[j], spacing[j], funcRHS,  u_0t, u_1t, u_x0, a)
    v_unrefined[j] = v
end

ratios, v_restrictedsCNH = refinement_ratios(v_unrefined, spacing, spacing)
writecsv("hw_1_refinement.txt", ratios)


##### Creates a heatmap of the Crank Nicolson solution for sucessively refined grids
######### Sub plot (a)
fig1_11 = heatmap(linspace(0,1,size(v_unrefined[1],2)),1:-1/size(v_unrefined[1],1):0,v_unrefined[1], title="(a) h, k = 1/2^3")
######### Sub plot (b)
fig1_12 = heatmap(linspace(0,1,size(v_unrefined[2],2)),1:-1/size(v_unrefined[2],1):0,v_unrefined[2], title="(b) h, k = 1/2^4")
######### Sub plot (c)
fig1_21 = heatmap(linspace(0,1,size(v_unrefined[3],2)),1:-1/size(v_unrefined[3],1):0,v_unrefined[3],xlabel="t",ylabel="x", title="(c) h, k = 1/2^5")
######### Sub plot (d)
fig1_22 = heatmap(linspace(0,1,size(v_unrefined[7],2)),1:-1/size(v_unrefined[7],1):0,v_unrefined[7], title="(d) h, k = 1/2^10")

plot(fig1_11, fig1_12, fig1_21, fig1_22)

##############################################################

#= Homework 1, problem 3
This is the setup and the function call to solve the heat equation

u_t = 1u_{xx}, 0 < x < 1
u(0,t) = 1,
u(1,t) = 0
u(x,0) = 1 if x < 0.5, 0 if x >= 0.5

with Crank-Nicolson and BDF-2

=#

h = 0.02 # grid spacing
k = 0.1 # time stepping
funcRHS( x, t) = 0 #
u_0t(t) = 1 # boundary condition at x = 0
u_1t(t) = 0 # boundary condition at x = 1
u_x0(x) = x < 0.5 ? 1 : 0 # initial condtion at t = 0
a =  1.0
#v_1  = PDEtool.CN_heat(h, k, funcRHS,  u_0t, u_1t, u_x0, a)
m = PDEtool.BDF2_heat(h, k, funcRHS,  u_0t, u_1t, u_x0, a)
heatmap(linspace(0,1,size(m,2)),1:-1/size(m,1):0,m, xlabel="t",ylabel="x", title="h = 0.02, k = 0.1")

# This code block runs Crank Nicolson for a selection of grid spacings and then makes a heat map for each run combined in one figure.
v_1 = PDEtool.CN_heat(h, k, funcRHS,  u_0t, u_1t, u_x0, a)
v_2 = PDEtool.CN_heat(h, k/10, funcRHS,  u_0t, u_1t, u_x0, a)
v_3 = PDEtool.CN_heat(h, k/100, funcRHS,  u_0t, u_1t, u_x0, a)
v_4 = PDEtool.CN_heat(h/10, k, funcRHS,  u_0t, u_1t, u_x0, a)

######### Sub plot (a)
fig1_11 = heatmap(linspace(0,1,size(v_1,2)),1:-1/size(v_1,1):0,v_1, title="(a) h = 0.02, k = 0.1")
######### Sub plot (b)
fig1_12 = heatmap(linspace(0,1,size(v_2,2)),1:-1/size(v_2,1):0,v_2, title="(b) h = 0.02, k = 0.01")
######### Sub plot (c)
fig1_21 = heatmap(linspace(0,1,size(v_3,2)),1:-1/size(v_3,1):0,v_3,xlabel="t",ylabel="x", title="(c) h = 0.02, k = 0.001")
######### Sub plot (d)
fig1_22 = heatmap(linspace(0,1,size(v_4,2)),1:-1/size(v_4,1):0,v_4, title="(d) h = 0.002, k = 0.1")
plot(fig1_11, fig1_12, fig1_21, fig1_22)
#savefig(plot(fig1_11, fig1_12, fig1_21, fig1_22) , "hw1_fig2.png")
