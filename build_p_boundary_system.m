function [P_bc_matrix,P_bc_rhs]= build_p_boundary_system(dx,dy,index_backstep,index_bottom,index_inlet, index_outlet,index_top,dp_dx_backstep,dp_dy_bottom,dp_dx_inlet,dp_dy_top,p_outlet,N1_y,N_y,N)
% Compute Matrix & Rhs of System for boundary nodes
P_bc_rhs = zeros(N,1);
 
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
   % Boundary conditions
    % inlet
 
        % case 'Neumann
            %for neumann at the inlet (upwind scheme):
            %set diagonal to -3/(2*dx)
            diagonal(index_inlet) = ones(length(index_inlet),1)*(-3/(2*dx));
            %set east1 to 4/(2*dx)
            east1(index_inlet) = ones(length(index_inlet),1)*(4/(2*dx));
            %set easteast1  to -1/(2*dx)
            easteast1(index_inlet) = ones(length(index_inlet),1)*(-1/(2*dx));
            %set right hand side to dp_dx_inlet
            P_bc_rhs(index_inlet) = dp_dx_inlet;
 
 
    % top
 
    % case 'Neumann
            %for neumann at the top (upwind scheme):
            %set diagonal to -3/(2*dy)
            diagonal(index_top) = ones(length(index_top),1)*(-3/(2*dy));
            %set south to 4/(2*dy)
            south(index_top) = ones(length(index_top),1)*(4/(2*dy));
            %set southsouth  to -1/(2*dy)
            southsouth(index_top) = ones(length(index_top),1)*(-1/(2*dy));
            %set right hand side to dp_dy_top
            P_bc_rhs(index_top) = dp_dy_top;
 
    % bottom
 
        %case 'Neumann'
            %for neumann at the southern border (downwind scheme):
                %set diagonal to 3/(2*dy)
            diagonal(index_bottom) = ones(length(index_bottom),1)*(3/(2*dy));
            %set north to -4/(2*dy)
            north(index_bottom) = ones(length(index_bottom),1)*(-4/(2*dy));
            %set northnorth to 1/(2*dy)
            northnorth(index_bottom) = ones(length(index_bottom),1)*(1/(2*dy));
            %set right hand side to dp_dy_bottom
            P_bc_rhs(index_bottom) = dp_dy_bottom;
%backstep
%case Neumann
    %for neumann at the backstep (upwind scheme):
            %set diagonal to -3/(2*dx)
            diagonal(index_backstep) = ones(length(index_backstep),1)*(-3/(2*dx));
            %set east3 to 4/(2*dx)
            east3(index_backstep) = ones(length(index_backstep),1)*(4/(2*dx));
            %set easteast3  to -1/(2*dx)
            easteast3(index_backstep) = ones(length(index_backstep),1)*(-1/(2*dx));
            %set right hand side to dp_dx_backstep
            P_bc_rhs(index_backstep) = dp_dx_backstep;
 
    
    %outlet
 
        %case 'Dirichlet'
            %for dirichlet: Set diagonal entry corresponding to point P to 1 and right hand side to p_outlet
            diagonal(index_outlet) = ones(length(index_outlet), 1);
            P_bc_rhs(index_outlet) = ones(length(index_outlet), 1)*p_outlet;
 
 
 
    %------------------------------------------------------------------------
    % Constructing Matrix P_bc_matrix
    
    %diagonal vector for spdiags
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
    P_bc_matrix = spdiags(B, d, N, N);
end

