function [u_star , v_star] = compute_uv_star(Auu, Auv, Avu, Avv, bu,bv,U_bc_matrix,U_bc_rhs,V_bc_matrix,V_bc_rhs,N)
%

A = [ (Auu+U_bc_matrix) , Auv ; Avu , (Avv+V_bc_matrix)];
b = [ bu + U_bc_rhs ; bv + V_bc_rhs];
star = A\b;
u_star = star(1:N);
v_star = star(N+1:end);
end

