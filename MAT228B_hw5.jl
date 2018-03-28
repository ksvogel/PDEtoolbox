#=
Kaela S. Vogel
2018-03-25
Homeworks for MAT 228

=#
include("PDEtool.jl")
using PDEtool
using PyPlot
using LatexPrint



##############################################################

#= Homework 5, problem 1
Solve the 1D acoustic equations using Lax-Wendroff.  We were given to solve it on (0,1), but no parameter values.  I chose the parameter values so that the wave speed was less than one.
=#
C = 0.9 # Courant number
h = 1/2^8 # grid spacing
k = C*h # time step
K = 0.5 # Bulk Modulus
r = 1.0 #
Tend = 3.0
u_x0(x)=cos.(16*pi*x).*exp.(-100*(x - 0.5)^2)
p_x0(x)=cos.(16*pi*x).*exp.(-100*(x - 0.5)^2)
u, p = acoustics1D(h, k, p_x0, u_x0, K, r, Tend)

# Plotting
fig = figure("LW_acoutstic",figsize=(10,10))
x = collect(0:h:(1-h)) + 0.5*h
snaps = [ 1 10 20 30 40 100 150 200 250]
subplot(339) # Create the third plot of a 4x4 group of subplots
suptitle("Linearized acoustic equations") # Supe title, title for all subplots combined
for sn in 1:1:9
    sp = parse(Int64, string("33", sn))
    subplot(sp)
    titlestr = string("t = ", k*snaps[sn])
    title(titlestr)
    ax = gca()
    setp(ax[:get_xticklabels](),visible=false) # Disable x tick labels
    setp(ax[:get_yticklabels](),visible=false) # Disable y tick labels
    ylim((-1,1))
    y1 = p[2:end-1,snaps[sn]]
    y2 = u[2:end-1,snaps[sn]]
    y = [y1 y2]

    #labels = [string("Time = ", kmesh[1]) string("Time = ", kmesh[50])     string("Time = ", kmesh[100])  string("Time = ", kmesh[end])  ]
    if sn == 1
        plot(x, y1, linewidth=1, alpha=0.6, label = "p(x,t)")
        plot(x, y2, linewidth=1, alpha=0.6, label = "u(x,t)")
        legend(loc="upper left",fancybox="true")
    else
        plot(x, y, linewidth=1, alpha=0.6)
    end
    fig[:canvas][:draw]()
end # Update the figure, required when doing additional modifications

############### Watch the waves travel and an "animation"
for sn in 1:1:600
    #fig = figure("LW_acoutstic",figsize=(10,10))
    x = collect(0:h:(1-h)) + 0.5*h
    y1 = p[2:end-1, sn]
    y2 = u[2:end-1, sn]
    y = [y1 y2]
        ylim((-1.5,1.5))
    #labels = [string("Time = ", kmesh[1]) string("Time = ", kmesh[50])     string("Time = ", kmesh[100])  string("Time = ", kmesh[end])  ]
    plot(x, y, linewidth=1, alpha=0.6)
    pause(0.005)
    clf()
    #fig[:canvas][:draw]()
end


##############################################################
###################################################################

#= Homework 5, problem 2
Solve the 1D acoustic equations using Lax-Wendroff.  We were given to solve it on (0,1), but no parameter values.  I chose the parameter values so that the wave speed was less than one.
=#
C = 0.9 # Courant number
a = 1.0
h = 1/2^8 # grid spacing
k = C*h # time step
Tend = 5.0
u_x0(x)=cos.(16*pi*x).*exp.(-50*(x - 0.5)^2)
u_x0(x)=sin.(2*pi*x)*sin.(4*pi*x)
u_x0(x) = abs(x - 0.5) < 0.25 ? 1 : 0
vs = Vector{Any}(7)
#v = advectionFV(h, k, a, u_x0, FL, Tend)
######## Loop through all the limiter function options
for j in 1:7
    FL = j
    vs[j] = advectionFV(h, k, a, u_x0, FL, Tend)
end



###############################################################################
# Plotting each limiter function per each initial condition

FLtype = ["upwinding" "Lax-Wendroff" "Beam-Warming" "minmod" "superbee" "MC" "van Leer"]
fig = figure("limitercompare",figsize=(10,10))
x = collect(0:h:(1-h)) + 0.5*h        ##### Annotate infty norm difference
subplot(337) # Create the third plot of a 4x4 group of subplots
#suptitle(L"Initial condition: $u(x,0) = \cos(16 \pi x)\exp(-50(x-0.5)^2)$")
#suptitle(L"Initial condition: $u(x,0) = \sin(2 \pi x) \sin(4 \pi x)$")
suptitle(L"Initial condition: $u(x,0) = \mathcal{X}_{[\frac{ 1}{4 }, \frac{ 3}{ 4}]}$")

for sn in 1:1:7
    v = vs[sn]
    sp = parse(Int64, string("33", sn))
    subplot(sp)
    titlestr = string("Limiter function: ", FLtype[sn])
    title(titlestr)
    ax = gca()
    setp(ax[:get_xticklabels](),visible=false) # Disable x tick labels
    setp(ax[:get_yticklabels](),visible=false) # Disable y tick labels
    ylim((-1.5,1.5))
    y1 = v[3:end-2, 1]
    y2 = v[3:end-2, end]
    y = [y1 y2]

    ##### Annotate infty norm difference
    diff = signif(vecnorm(v[3:end-2, 1] - v[3:end-2, end], Inf),5)
    strdiff = string(diff)

    annotate(L"$|| v^0 - v^N ||_{\infty} = $",
	xy=[0;0],
	xycoords="axes fraction",
	xytext=[10,10],
	textcoords="offset points",
	fontsize=12.0,
	ha="left",
	va="bottom")

    annotate(strdiff,
	xy=[1;0],
	xycoords="axes fraction",
	xytext=[-10,10],
	textcoords="offset points",
	fontsize=12.0,
	ha="right",
	va="bottom")


    if sn == 1
        plot(x, y1, linewidth=1, alpha=0.6, label = L"$v^0$")
        plot(x, y2, linewidth=1, alpha=0.6, label = L"$v^N$")
        legend(loc="upper left",fancybox="true")
    else
        plot(x, y, linewidth=1, alpha=0.6)
    end
    fig[:canvas][:draw]()
end # Update the figure, required when doing additional modifications

############### Watch the waves travel and an "animation"
for sn in 1:10:1420
    v = vs[6]
    #fig = figure("LW_acoutstic",figsize=(10,10))
    x = collect(0:h:(1-h)) + 0.5*h
    y = v[3:end-2, sn]
    ylim((-1.5,1.5))
    #labels = [string("Time = ", kmesh[1]) string("Time = ", kmesh[50])     string("Time = ", kmesh[100])  string("Time = ", kmesh[end])  ]
    plot(x, y, linewidth=1, alpha=0.6)
    pause(0.0005)
    clf()
    #fig[:canvas][:draw]()
end
