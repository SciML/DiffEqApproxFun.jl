using Base.Test, DiffEqApproxFun

tic()
@time @testset "Direct ApproxFun Test" begin include("direct_approxfun.jl") end
@time @testset "Indirect ApproxFun Test" begin include("indirect_approxfun.jl") end
toc()
