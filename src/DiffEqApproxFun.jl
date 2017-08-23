__precompile__()

module DiffEqApproxFun

using Reexport
@reexport using ApproxFun
@reexport using DiffEqBase

include("problem_solution_types.jl")
include("boundary_handling.jl")

export ApproxFunProblem, BoundaryCallback

end # module
