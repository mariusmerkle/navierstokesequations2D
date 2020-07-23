function pc = compute_pressure_correction(density,dt,Dx,Dy,L,u_star,v_star,P_bc_matrix,P_bc_rhs)

    rhs = density/dt*(Dx*u_star+Dy*v_star)+P_bc_rhs;
    pc = (L+P_bc_matrix)\rhs;
end

