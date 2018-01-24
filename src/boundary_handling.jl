function BoundaryCallback(bc)
  affect! = function (integrator)
    S=Fourier()
    u = Fun(S,integrator.u)
    tmp = bc(integrator.t,u)
    integrator.u=pad!(tmp.coefficients,length(integrator.u))
  end
  DiscreteCallback((u,t,integrator)->true,affect!)
end
