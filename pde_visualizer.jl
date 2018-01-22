#=
Kaela S. Vogel
2018-01-07
Homeworks for MAT 228

Visualizations of solved PDEs
=#

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
