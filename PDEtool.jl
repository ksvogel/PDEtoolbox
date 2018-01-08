#= PDEtool module
# put
include("PDEtool.jl")
using PDEtool
at top of file to avoid having to use the module.function notation

=#
push!(LOAD_PATH, "/home/kaela/Documents/GithubRepositories/PDE_solvers")

module PDEtool

include("poisson_iterative.jl")
export h_cycle
export gauss_sidel
export SOR
export jacobi_iter
export laplace_mat2D

include("MG_vcycle.jl")
export squish
export foomp
export vcycle

include("badbehavior.jl")
export recursioncheck

include("krylov_methods.jl")
export PCG
export preconditioner


include("finitedifferencemethods_hyperbolic_parabolic.jl")
export meshmaker
export CN_heat

end
