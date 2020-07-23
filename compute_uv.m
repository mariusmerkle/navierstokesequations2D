function [u, v] =  compute_uv(u_star,v_star,dt,density,Dx,Dy,pc)
%
u = u_star - dt/density*Dx*pc;
v = v_star - dt/density*Dy*pc;

end

