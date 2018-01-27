using DiffEqBase, OrdinaryDiffEq, ApproxFun, DiffEqApproxFun
using Base.Test

S=Fourier()
u0=Fun(θ->cos(cos(θ-0.1))-cos(cos(0-0.1)),S)
c=Fun(cos,S)

prob = ODEProblem((u,p,t)->u''+(c+1)*u',u0,(0.,1.))

@time sol=solve(prob,Euler(),dt=1/1000)

tstops =  [0.0,0.00630546,0.0196241,0.0365809,0.0580536,0.0840865,0.1151,
           0.151529,0.195094,0.24658,0.292971,0.316851,0.338594,0.36075,
           0.384077,0.408548,0.433682,0.458922,0.483926,0.508619,0.533097,
           0.557498,0.581922,0.606414,0.630967,0.655551,0.68014,0.704717,
           0.729278,0.753828,0.778375,0.802923,0.827475,0.85203,0.876586,
           0.901142,0.925698,0.950253,0.974808,0.999362,1.0]

@time sol=solve(prob,Tsit5(),dt=1/40, adaptive=false, tstops=tstops)

#@time sol=solve(prob,Tsit5())
#@time sol=solve(prob,Rosenbrock23(),dt=1/40, adaptive=false, tstops=tstops)
