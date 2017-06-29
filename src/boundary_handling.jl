function BoundaryCallback(bc)
  affect! = function (integrator)
    S=Fourier()
    u = Fun(S,integrator.u)
    tmp = bc(integrator.t,u)
    integrator.u=pad!(tmp.coefficients,length(integrator.u))
  end
  DiscreteCallback((t,u,integrator)->true,affect!)
end
