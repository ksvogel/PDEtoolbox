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

include("MG_vcycle.jl")
export squish
export foomp
export vcycle

include("badbehavior.jl")
export parentalunit
export evilchild
export recursioncheck
end
