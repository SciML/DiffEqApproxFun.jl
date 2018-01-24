using DiffEqBase, OrdinaryDiffEq, Sundials, ApproxFun, DiffEqApproxFun
using Base.Test

S=Fourier()
u0=Fun(θ->cos(cos(θ-0.1))-cos(cos(0-0.1)),S)
c=Fun(cos,S)


ode_prob = ODEProblem((u,p,t)->u''+(c+1)*u',u0,(0.,1.))
prob = ApproxFunProblem(ode_prob)
@time sol=solve(prob,Euler(),dt=1/1000)
@time sol=solve(prob,Tsit5())
sol(0.5)
sol(0.5,0.2)

@time sol=solve(prob,CVODE_BDF())

sol[1]
sol[5](0.2)

function bc(t,u)
  C=eye(S)[3:end,:]
  tmp = [Evaluation(0);
         Evaluation(2π);
         C]\[0.;0.;u]
end

u0=Fun(θ->cos(cos(θ)) - cos(cos(0)),S)
u0 = bc(0.0,u0)
c=Fun(cos,S)
ode_prob = ODEProblem((u,p,t)->u''+(c+1)*u',u0,(0.,1.))

prob = ApproxFunProblem(ode_prob)

condition(t,u,integrator) = true
function affect!(integrator)
  S=Fourier()
  u = Fun(S,integrator.u)
  B=Dirichlet()
  C=eye(S)[3:end,:]
  tmp = [Evaluation(0);
         Evaluation(2π);
         C]\[0.;0.;u]
  integrator.u=pad!(tmp.coefficients,length(integrator.u))
end
boundary_projection = DiscreteCallback(condition,affect!)

@time sol=solve(prob,Tsit5())
@test sol(0.0,0.0) ≈ 0.0 atol=1e-12
@test sol(0.0,2π) ≈ 0.0 atol=1e-12
@test sol(0.5,0.0) > 0.3
@test sol(0.5,2π) > 0.3
@test sol(1.0,0.0) > 0.3
@test sol(1.0,2π) > 0.3

@time sol=solve(prob,Tsit5(),callback=boundary_projection)
@test sol(0.0,0.0) ≈ 0.0 atol=1e-12
@test sol(0.0,2π) ≈ 0.0 atol=1e-12
@test sol(0.5,0.0) ≈ 0.0 atol=1e-3
@test sol(0.5,2π) ≈ 0.0 atol=1e-3
@test sol(1.0,0.0) ≈ 0.0 atol=1e-3
@test sol(1.0,2π) ≈ 0.0 atol=1e-3

boundary_projection2 = BoundaryCallback(bc)

@time sol=solve(prob,Tsit5(),callback=boundary_projection2)
@test sol(0.0,0.0) ≈ 0.0 atol=1e-12
@test sol(0.0,2π) ≈ 0.0 atol=1e-12
@test sol(0.5,0.0) ≈ 0.0 atol=1e-3
@test sol(0.5,2π) ≈ 0.0 atol=1e-3
@test sol(1.0,0.0) ≈ 0.0 atol=1e-3
@test sol(1.0,2π) ≈ 0.0 atol=1e-3
