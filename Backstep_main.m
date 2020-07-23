%% run Backstep solver
clear all
close all
clc
% profile 'flat' 'parabolic' 'triangular'
tic
[X1,Y1,X2,Y2, index_top, index_bottom,index_outlet, index_backstep,index_inlet, u,v,p,E_kin,divergence,divergence_norm,Re,dt,a,profile,density,dynamic_viscosity] = backstep_incremental_pressure_correction(1,50,0.001,5,'parabolic');
visualization(X1,Y1,X2,Y2, index_top, index_bottom,index_outlet, index_backstep,index_inlet, u,v,p,E_kin,divergence,divergence_norm,Re,dt,a,profile,density,dynamic_viscosity);
save([profile,'_Re_',num2str(Re),'_dt_',num2str(dt),'_a_',num2str(a),'.mat']);
toc