function [U_bc_matrix,V_bc_matrix,U_bc_rhs,V_bc_rhs] = build_uv_boundary_system(dx,index_backstep,index_bottom,index_inlet, index_outlet,index_top, u_backstep,u_bottom,u_inlet,u_top,du_dx_outlet,v_backstep,v_bottom,v_inlet,v_top,dv_dx_outlet,N1_y,N_y,N)
%Compute Matrix & Rhs of System for boundary nodes
    U_bc_rhs = zeros(N,1);
    V_bc_rhs = zeros(N,1);
% Defining Boundary conditions
    % Initialize Diagonal Vectors for spdiags
        diagonal = zeros(N,1);
        %west
        west1 = zeros(N,1);
        west3 = zeros(N,1);
 
        westwest1 = zeros(N,1);
        westwest3 = zeros(N,1);
 
        %east
        east1 = zeros(N,1);
        east3 = zeros(N,1);
 
        easteast1 = zeros(N,1);
        easteast3 = zeros(N,1);
 
        %north
        north = zeros(N,1);
        northnorth = zeros(N,1);
        %south
        south = zeros(N,1);
        southsouth = zeros(N,1);
 
 
    %inlet
    % case Dirichlet
    diagonal(index_inlet) = ones(length(index_inlet),1);
    U_bc_rhs(index_inlet) = u_inlet;
    V_bc_rhs(index_inlet) = v_inlet;
    
    %top
    % case Dirichlet
    diagonal(index_top) = ones(length(index_top),1);
    U_bc_rhs(index_top) = u_top;
    V_bc_rhs(index_top) = v_top;
 
    %bottom
    % case Dirichlet
    diagonal(index_bottom) = ones(length(index_bottom),1);
    U_bc_rhs(index_bottom) = u_bottom;
    V_bc_rhs(index_bottom) = v_bottom;

 
    %backstep
    % case Dirichlet
    diagonal(index_backstep) = ones(length(index_backstep),1);
    U_bc_rhs(index_backstep) = u_backstep;
    V_bc_rhs(index_backstep) = v_backstep;

 
    %outlet
    %case Neumann
            %for neumann at the outlet (downwind scheme):
            %set diagonal to -3/(2*dx)
            diagonal(index_outlet) = ones(length(index_outlet),1)*(3/(2*dx));
            %set west3 to -4/(2*dx)
            west3(index_outlet) = ones(length(index_outlet),1)*(-4/(2*dx));
            %set westwest3  to +1/(2*dx)
            westwest3(index_outlet) = ones(length(index_outlet),1)*(1/(2*dx));
            %set right hand side to du_dx_outlet
            U_bc_rhs(index_outlet) = du_dx_outlet;
            V_bc_rhs(index_outlet) = dv_dx_outlet;
 
    %------------------------------------------------------------------------
    % Constructing Matrix U_bc_matrix & V_bc_matrix;
 
    %off diagonal vectors for spdiags ... circshift is needed for spdiags to
    %work, since A is mxn with m=n but we need the behavior of spdiags for m<n
    %more info in the matlab documentation on spdiags
    %west
        west1 = circshift(west1,-N1_y);
        westwest1 = circshift(westwest1,-2*N1_y);
        west3 = circshift(west3,-N_y);
        westwest3 = circshift(westwest3,-2*N_y);
 
    %east
    east1 = circshift(east1, N1_y);
    easteast1 = circshift(easteast1,2*N1_y);
    east3 = circshift(east3,N_y);
    easteast3 = circshift(easteast3,2*N_y);
 
   %north    
    north = circshift(north,-1);
    northnorth = circshift(northnorth,-2);
 
    %south
    south = circshift(south,1);
    southsouth = circshift(southsouth,2);
 
 
   %define matrix with vectors for the diagonals specified in d
    B = [ westwest3 westwest1 west3 west1 northnorth north diagonal south southsouth east1 east3 easteast1 easteast3];
    d = [-2*N_y, -2*N1_y, -N_y, -N1_y, -2, -1, 0, 1, 2, N1_y, N_y, 2*N1_y, 2*N_y];
    U_bc_matrix= spdiags(B, d, N, N);
    
    V_bc_matrix = U_bc_matrix;
end

