function [X1,Y1,X2,Y2, index_top, index_bottom,index_outlet, index_backstep,index_inlet, u,v,p,E_kin,divergence,divergence_norm,Re,dt,a,profile,density,dynamic_viscosity] = backstep_incremental_pressure_correction(a, Re, dt, T,profile)
%% Main script for Navier Stokes
% Backstep flow solver in 2D
% Standard Incremental Pressure Correction Scheme 
% Reference: Discretization of the Navier-Stokes equations on a co-located grid (Cesar Paez, 2019) 

%% Input Parameters

density = 1260; %Glycerol %997; %Water % [kg/m^3] 
dynamic_viscosity = 1; %Glycerol %8.90e-4; %Water % [Pa*s] 
kinematic_viscosity = dynamic_viscosity/density; %[m^2/s]


%{
                                    l_x
-----------------------------------------------------------------------
 
h1_y
        l1_x
-------------------------                                               h_y
                        |
                        |h2_y
                        |               l2_x
                        -----------------------------------------------
%}
%% Geometric Parameters: fully determined by h1_y
% a defined in function head gives height of step compared to inlet height h1_y
c = 1; % gives length of l1_x compared to inlet height h1_y
d = 6; % gives length of l2_x compared to outlet channel height h_y should be around 10 for fully developed flow
h1_y = 0.05; % [m] inlet channel height
h2_y = a*h1_y; % [m]
h_y = h1_y + h2_y; % [m]
l1_x = c*h1_y; %  [m]
l2_x = d*h_y; % [m]
l_x = l1_x + l2_x; % [m]
 
 
%% Discretization Parameters: fully determined by N1_x and N1_y
N1_x = 15; % has to be odd
N1_y = 15; % has to be odd
N2_y = a*(N1_y-1);
N2_x = d/c*(a+1)*(N1_x - 1);
N_y = N1_y + N2_y;
N_x = N1_x + N2_x;
N1 = N1_y*N1_x;
N2 = N_y*N2_x;
N = N1 + N2_y + N2;
if(mod(N2_y,1) ~= 0 || mod(N2_x,1) ~= 0)
    error('N2_y and N2_x have to be integers, set a,c,d accordingly');
end
 
%% Index Vectors
%index_inlet = 2:N1_y-1;
index_inlet = 1:N1_y;

%index_bottom_1 = N1_y:N1_y:(N1);
index_bottom_1 = 2*N1_y:N1_y:(N1);
index_bottom_2 = (N1 + N2_y):N_y:N;
index_bottom = [index_bottom_1, index_bottom_2];
%index_top_1 = 1:N1_y:(N1- (N1_y-1));
index_top_1 = N1_y+1:N1_y:(N1- (N1_y-1));
index_top_2 = (N1 + N2_y + 1):N_y:(N - (N_y-1));
index_top = [index_top_1, index_top_2];
index_backstep = (N1+1):1:(N1+N2_y-1);
index_outlet = (N-N_y+2):1:N-1;
 
 
 
%% Convergence Parameters
minEr = 1e-6;
 
%% Domain
[X1,Y1] = meshgrid(linspace(0,l1_x,N1_x),linspace(h_y,h2_y, N1_y));
[X2,Y2] = meshgrid(linspace(l1_x,l_x,N2_x+1),linspace(h_y,0, N_y));

dx = X1(1, 2) - X1(1, 1);
dy = Y1(1, 1) - Y1(2, 1);
  


%% Initialisation
u = zeros(N,1);
v = zeros(N,1);
p = zeros(N,1); 
 
 
%% Boundary Conditions
% for Poiseuille Flow
% Re = u_avg*h_y/kinematic_viscosity
u_avg = Re*kinematic_viscosity/(h1_y);
u_max = 3/2*u_avg;


y_inlet = Y1(1:end,1);
% X-direction
switch profile
    case 'flat'
        u_inlet = u_avg;
    case 'parabolic'
        %dp_dx_poiseuille = -8*dynamic_viscosity*u_max/(h1_y^2)
        u_inlet = -h1_y^2/(2*dynamic_viscosity)*(-8*dynamic_viscosity*u_max/(h1_y^2)).*(y_inlet-h2_y)/h1_y.*(1 - (y_inlet-h2_y)/h1_y); % [m/s]
    case 'triangular'
        
        
end
        %

u_bottom = 0;
u_top = 0;
u_backstep = 0;
du_dx_outlet = 0;
 
% Y-direction
v_inlet = 0;
v_bottom = 0;
v_top = 0;
v_backstep = 0;
dv_dx_outlet = 0;
 
% Pressure
dp_dx_inlet = 0;
dp_dy_top = 0;
dp_dx_backstep = 0;
dp_dy_bottom = 0;
p_outlet = 0;
disp(['inlet max velocity = ', num2str(u_max),' m/s with Reynolds number of Re = ',num2str(Re)]);
 
u(index_inlet) = u_inlet;

%% Stability 
delta = kinematic_viscosity*dt/(min(dx,dy)^2);
%dt_delta = delta*(min(dx,dy)^2)/kinematic_viscosity;
courant = u_avg*(dt/min(dx,dy));
%dt_courant = courant*min(dx,dy)/max(u_inlet);
%dt = min(dt_courant,dt_delta);
disp(['convergence parameters: diffusion number = ',num2str(delta),' and courant number = ',num2str(courant)]);
disp(['time step dt = ',num2str(dt)]);

%% Solving
inner_nodes = ones(N,1);
inner_nodes(index_backstep) = 0;
inner_nodes(index_bottom) = 0;
inner_nodes(index_top) = 0;
inner_nodes(index_outlet) = 0;
inner_nodes(index_inlet) = 0;

[Dx,Dy,L] = build_discrete_operators(dx,dy,inner_nodes, N, N1, N1_y,N_y);
[U_bc_matrix,V_bc_matrix,U_bc_rhs,V_bc_rhs] = build_uv_boundary_system(dx,index_backstep,index_bottom,index_inlet, index_outlet,index_top, u_backstep,u_bottom,u_inlet,u_top,du_dx_outlet,v_backstep,v_bottom,v_inlet,v_top,dv_dx_outlet,N1_y,N_y,N);

[P_bc_matrix,P_bc_rhs]= build_p_boundary_system(dx,dy,index_backstep,index_bottom,index_inlet, index_outlet,index_top,dp_dx_backstep,dp_dy_bottom,dp_dx_inlet,dp_dy_top,p_outlet,N1_y,N_y,N);

maxIt = round(T/dt);

iter = 0;
E_kin = zeros(maxIt,1);
divergence_norm = zeros(maxIt,1);

while(iter < maxIt)
    [Auu, Auv, Avu, Avv, bu,bv] = compute_uv_star_system(dt,N,Dx,Dy,L,u,v,p,density,kinematic_viscosity);
    [u_star , v_star] = compute_uv_star(Auu, Auv, Avu, Avv, bu,bv,U_bc_matrix,U_bc_rhs,V_bc_matrix,V_bc_rhs,N);
    pc = compute_pressure_correction(density,dt,Dx,Dy,L,u_star,v_star,P_bc_matrix,P_bc_rhs);
    [u, v] =  compute_uv(u_star,v_star,dt,density,Dx,Dy,pc);
    p = p + 1*pc;
    divergence = Dx*u + Dy*v;
    iter = iter +1;
    divergence_norm(iter) = norm(divergence);
    E_kin(iter) = 1/2*density*[u;v]'*[u;v];
    if(isnan(E_kin(iter))) 
        warning('kinetic energy diverged!! Solution sequence stopped')
        break;
    end
    if (((iter > 1) && (divergence_norm(iter) < minEr) ))
        break;
    end
end
E_kin = E_kin(1:iter);
disp(['iterations needed: ', num2str(iter)]);
disp(['with a relative kinetic energy error of: ', num2str((E_kin(iter)-E_kin(iter-1))/E_kin(iter))]);
 



end

