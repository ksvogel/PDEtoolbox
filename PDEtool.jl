# PDEtool module

push!(LOAD_PATH, "/home/kaela/Documents/GithubRepositories/PDE_solvers")

module PDEtool

include("poisson_iterative.jl")
export h_cycle
export gauss_sidel

include("MG_vcycle.jl")
export squish
export foomp
end
