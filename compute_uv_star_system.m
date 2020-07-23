function [Auu, Auv, Avu, Avv, bu,bv] = compute_uv_star_system(dt,N,Dx,Dy,L,u,v,p,density,kinematic_viscosity)
%
I = eye(N,N);
Auu = 1/dt*I+diag(Dx*u)-kinematic_viscosity*L;
Auv = diag(Dy*u);
Avu = diag(Dx*v);
Avv = 1/dt*I+diag(Dy*v)-kinematic_viscosity*L;
bu = 1/dt * u - 1/density*Dx*p;
bv = 1/dt * v - 1/density*Dy*p;
end

