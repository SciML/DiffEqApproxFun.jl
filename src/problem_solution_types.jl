immutable ApproxFunProblem{P,S} <: DEProblem
    prob::P
    space::S
end

function ApproxFunProblem(prob)
  sp = space(prob.u0)
  n = ncoefficients(prob.u0)
  _prob = ODEProblem((t,u)->pad!(prob.f(t,Fun(sp,u)).coefficients,n),prob.u0.coefficients,prob.tspan)
  return ApproxFunProblem(_prob,sp)
end

immutable ApproxFunSolution{SOL,S}
    sol::SOL
    space::S
end

(sol::ApproxFunSolution)(t) = Fun(sol.space,sol.sol(t))
(sol::ApproxFunSolution)(t,x) = sol(t)(x)

@inline Base.getindex(sol::ApproxFunSolution, I::Int...) = Fun(sol.space,sol.sol[I[end]])[Base.front(I)...]
@inline Base.getindex(sol::ApproxFunSolution, I::Int) = Fun(sol.space,sol.sol[I])
@inline Base.getindex(sol::ApproxFunSolution, I::Colon) = Fun.(sol.space,sol.sol[I])
@inline Base.getindex(sol::ApproxFunSolution, I::AbstractArray{Int}) = Fun.(sol.space,sol.sol[I])

@inline Base.length(sol::ApproxFunSolution) = length(sol.sol)
@inline Base.size(sol::ApproxFunSolution) = size(sol.sol)
@inline Base.eachindex(sol::ApproxFunSolution) = Base.OneTo(length(sol.u))

function DiffEqBase.init(prob::ApproxFunProblem,alg,args...;kwargs...)
    init(prob.prob,alg,args...;kwargs...)
end
function DiffEqBase.solve(prob::ApproxFunProblem,alg,args...;kwargs...)
    sol = solve(prob.prob,alg,args...;kwargs...)
    ApproxFunSolution(sol,prob.space)
end
