function [Dx,Dy,L] = build_discrete_operators(dx,dy,inner_nodes, N, N1, N1_y,N_y)
%Builds Discrete Operator Matrices for use in 2D NS solution
% Dx
    e1 = [inner_nodes(1:(N1-N1_y)); zeros(N-(N1-N1_y),1)];
    e2 = [zeros(N1-N1_y, 1); (inner_nodes((N1-N1_y+1):N1)); zeros(N-N1, 1)];
    e3 = [zeros(N1,1);(inner_nodes(N1+1:N))];
    
    Dx = 1/(2*dx)*spdiags([ -circshift(e3,-N_y) -circshift(e1+e2,-N1_y) circshift(e1,N1_y) circshift(e2+e3,N_y) ], [-N_y, -N1_y, N1_y, N_y], N,N);
    Dy = 1/(2*dy)*spdiags([ 1*circshift(e1+e2+e3,-1) -1*circshift(e1+e2+e3,1)], [-1, 1], N,N);
    
    DxDx = 1/dx^2*spdiags([ circshift(e3,-N_y) circshift(e1+e2,-N1_y) -2*(e1+e2+e3) circshift(e1,N1_y) circshift(e2+e3,N_y)], [-N_y,-N1_y, 0 , N1_y,N_y], N,N);
    DyDy = 1/dy^2*spdiags([ circshift(e1+e2+e3,1) -2*(e1+e2+e3) circshift(e1+e2+e3,-1)], [1, 0 , -1], N,N);

    L = DxDx + DyDy;
end

