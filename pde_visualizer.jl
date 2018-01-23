#=
Kaela S. Vogel
2018-01-07
Homeworks for MAT 228

Visualizations of solved PDEs
=#
#using LaTeXStrings
#PGFPlots.pushPGFPlotsPreamble("\\usepackage{amssymb}")

# This code block runs Crank Nicolson for a selection of grid spacings and then makes a heat map for each run combined in one figure.
#=
v_1 = PDEtool.CN_heat(h, k, funcRHS,  u_0t, u_1t, u_x0, a)
v_2 = PDEtool.CN_heat(h, k/10, funcRHS,  u_0t, u_1t, u_x0, a)
v_3 = PDEtool.CN_heat(h, k/100, funcRHS,  u_0t, u_1t, u_x0, a)
v_4 = PDEtool.CN_heat(h/10, k, funcRHS,  u_0t, u_1t, u_x0, a)

hmesh, kmesh, F, M, T = meshmaker(funcRHS, h, k)

######### Sub plot (a)
anames1 = ["" for x=kmesh]
anames2 = ["" for x=hmesh]
anames1[1:2:T] = ["$(x)" for x in kmesh[1:2:T]]
anames2[1:8:M] = ["$(x)" for x in hmesh[M:-8:1]]
fig1_11 = heatmap(linspace(0,1,size(v_1,2)),1:-1/size(v_1,1):0,v_1,xlabel="foo",ylabel="bar",xticks= (kmesh, anames1), yticks= (hmesh, anames2), title="Meep")

######### Sub plot (b)
bnames1 = ["" for x=kmesh]
bnames2 = ["" for x=hmesh]
bnames1[1:2:T] = ["$(x)" for x in kmesh[1:2:T]]
bnames2[1:8:M] = ["$(x)" for x in hmesh[M:-8:1]]
fig1_12 = heatmap(linspace(0,1,size(v_2,2)),1:-1/size(v_2,1):0,v_2,xlabel="foo",ylabel="bar",xticks= (kmesh, bnames1), yticks= (hmesh, bnames2), title="Meep")

######### Sub plot (c)
bnames1 = ["" for x=kmesh]
bnames2 = ["" for x=hmesh]
bnames1[1:2:T] = ["$(x)" for x in kmesh[1:2:T]]
bnames2[1:8:M] = ["$(x)" for x in hmesh[M:-8:1]]
fig1_21 = heatmap(v_3)

######### Sub plot (d)
fig1_22 = heatmap(v_4)

plot(fig1_11, fig1_12, fig1_21, fig1_22)

#=


















#=
    # pyplot_surfaceplot.jl
    #
    #	Surface Plot demonstration
    #
    # Daniel Høegh (https://gist.github.com/dhoegh)
    # Julia 0.6.0
    # 09.12.2014
    # Last Edit: 20.07.17

    # Reference: https://groups.google.com/d/msg/julia-users/eVtZdp3htTM/TJOt3exCxKgJ

    using PyPlot
    using Distributions

    ###################
    ##  Create Data  ##
    ###################
    n = 100
    x = linspace(-3, 3, n)
    y = linspace(-3,3,n)

    xgrid = repmat(x',n,1)
    ygrid = repmat(y,1,n)

    z = zeros(n,n)

    for i in 1:n
    for j in 1:n
    z[i:i,j:j] = pdf(MvNormal(eye(2)),[x[i];y[j]])
    end
    end

    ############
    ##  Plot  ##
    ############
    fig = figure("pyplot_surfaceplot",figsize=(10,10))
    ax = fig[:add_subplot](2,1,1, projection = "3d")
    ax[:plot_surface](xgrid, ygridu
    =#
