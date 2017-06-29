__precompile__()

module DiffEqApproxFun

using ApproxFun, DiffEqBase

include("problem_solution_types.jl")
include("boundary_handling.jl")

export ApproxFunProblem, BoundaryCallback

end # module
