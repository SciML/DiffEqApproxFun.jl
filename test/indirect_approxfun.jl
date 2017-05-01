using DiffEqBase, OrdinaryDiffEq, Sundials, ApproxFun, DiffEqApproxFun
using ODE, ODEInterfaceDiffEq
using Base.Test

S=Fourier()
u0=Fun(θ->cos(cos(θ-0.1)),S)
c=Fun(cos,S)


ode_prob = ODEProblem((t,u)->u''+(c+1)*u',u0,(0.,1.))
prob = ApproxFunProblem(ode_prob)
@time sol=solve(prob,Tsit5())
sol(0.5)
sol(0.5,0.2)


@time sol=solve(prob,CVODE_BDF())
@time sol=solve(prob,ode45())
@time sol=solve(prob,radau())

sol[1]
sol[2,1]
