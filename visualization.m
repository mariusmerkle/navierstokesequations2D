function []= visualization(X1,Y1,X2,Y2, index_top, index_bottom,index_outlet, index_backstep,index_inlet, u,v,p,E_kin,divergence,divergence_norm,Re,dt,a,profile,density,dynamic_viscosity)
    params = [profile,'_Re_',num2str(Re),'_dt_',num2str(dt),'_a_',num2str(a),'_'];
    % Geometric Parameters
    h_y = Y1(1,1);
    h2_y = Y1(end,1);
    h1_y = h_y - h2_y;
    l1_x = X1(1,end);
    l_x = X2(1,end);
    % Mesh Parameters
    N1_x = size(X1,2);
    N1_y = size(X1,1);
    N1 = N1_y*N1_x;
    
    N2_x = size(X2,2)-1;
    N_y = size(X2,1);
    N2_y = N_y - N1_y;
    N2 = N_y*N2_x;
    
    N_x = N1_x + N2_x;
    N = N1 + N2_y + N2;
    
    boundary_bottom = [0, h2_y ; l1_x , h2_y ; l1_x , 0; l_x , 0];
    boundary_top = [0, h_y ; l1_x , h_y ; l_x , h_y];
    
    
    %% Hagen-Poiseuille Flow
    p_in = mean(p(index_inlet));
    p_out = mean(p(index_outlet));
    poiseuille =@(y) -h_y^2/(2*dynamic_viscosity)*((p_out-p_in)/l_x).*y/h_y.*(1 - y/h_y);
    %% Visualization
    x1 = X1(:);
    y1 = Y1(:);
    x2 = X2(:);
    y2 = Y2(:);
    x = [x1(1:N1-N1_y);x2];
    y = [y1(1:N1-N1_y);y2];
 
    index_bottom_1 = N1_y:N1_y:(N1);
    index_bottom_2 = (N1 + N2_y):N_y:N;
    index_top_1 = 1:N1_y:(N1- (N1_y-1));
    index_top_2 = (N1 + N2_y + 1):N_y:(N - (N_y-1));
    [X,Y] = meshgrid(linspace(0,l_x,N_x),linspace(h_y,0,N_y));
%     p_mat = zeros(size(X));
%     p_mat(1:N1_y,1:N1_x) = reshape(p(1:N1),size(X1));
%     p_mat(1:N_y,N1_x:end)= reshape(p(N1-N1_y+1:N),size(X2));
    u_mat = zeros(size(X));
    u_mat(1:N1_y,1:N1_x) = reshape(u(1:N1),size(X1));
    u_mat(1:N_y,N1_x:end)= reshape(u(N1-N1_y+1:N),size(X2));
    v_mat = zeros(size(Y));
    v_mat(1:N1_y,1:N1_x) = reshape(v(1:N1),size(Y1));
    v_mat(1:N_y,N1_x:end)= reshape(v(N1-N1_y+1:N),size(Y2));
    
    %Domain and Boundaries
    
    figure();
    pcolor(X1,Y1,ones(size(X1)));
    hold on
    pcolor(X2,Y2,ones(size(X2)));
    colorbar;
    axis([0 l_x -l_x/2 l_x/2]);
    plot(x(index_bottom_1),y(index_bottom_1),'r*');
    plot(x(index_top_1),y(index_top_1),'b*');
    plot(x(index_inlet),y(index_inlet),'y*');
    plot(x(index_top_2),y(index_top_2),'g*');
    plot(x(index_bottom_2),y(index_bottom_2),'m*');
    plot(x(index_outlet),y(index_outlet),'c*');
    plot(x(index_backstep),y(index_backstep),'k*');
    title('Mesh-Plot')
    xlabel('x');
    ylabel('y');
    frame = getframe(gcf);
    im = frame2im(frame);
    imwrite(im,[params,'Mesh','.jpeg']);

    %Pressure Field
    figure();
    contour(X1,Y1,reshape(p(1:N1),size(X1)),N1_x);
    hold on
    contour(X2,Y2,reshape(p(N1-N1_y+1:N),size(X2)),N2_x);
    colorbar;
    axis([0 l_x 0 h_y]);
    plot(boundary_bottom(:,1),boundary_bottom(:,2),'-k');
    plot(boundary_top(:,1),boundary_top(:,2),'-k');
    title('Contour Plot of Pressure Field')
    xlabel('x');
    ylabel('y');
    frame = getframe(gcf);
    im = frame2im(frame);
    imwrite(im,[params,'Contour_pressure','.jpeg']);
    

    figure()
    plot(E_kin);
    title('total kinetic Energy for each iteration');
    xlabel('iterations')
    ylabel('E_kin');
    frame = getframe(gcf);
    im = frame2im(frame);
    imwrite(im,[params,'Kinetic_Energy_per_Iteration','.jpeg']);
    

    figure();
    
    subplot(2, 1, 1);
    plot(u([index_top(end),index_outlet,index_bottom(end)]),'b-');
    hold on;
    plot(poiseuille(Y2(:,end)),'r-');
    title('Velocity Profile at the outlet & Poiseuille Profile');
    ylabel('m/s');
    xlabel('y');
    subplot(2, 1, 2);
    surf(X1, Y1, reshape(u(1:N1),size(X1)))
    hold on;
    surf(X2, Y2,reshape(u(N1-N1_y+1:N),size(X2)));
    view([-10,20]);
    frame = getframe(gcf);
    im = frame2im(frame);
    imwrite(im,[params,'Hagen-Poiseuille','.jpeg']);
    title('Velocity Profile in x-Direction')
    xlabel('x');
    ylabel('y');
    zlabel('m/s');

    figure()
    
    
    %s = stream2(X1,Y1,reshape(ux(1:N1),size(X1)),reshape(uy(1:N1),size(Y1)),X1(:,19),Y1(:,1));
%     streamline(s);
%     hold on
%     streamline(X2,Y2,reshape(ux(N1-N1_y+1:N),size(X2)),reshape(uy(N1-N1_y+1:N),size(Y2)),X2(:,1)+0.01,Y2(:,1));
    streamline(X,Y,u_mat,v_mat,X1(:,1),Y1(:,1));
    hold on
    plot(boundary_bottom(:,1),boundary_bottom(:,2),'-k');
    plot(boundary_top(:,1),boundary_top(:,2),'-k');
    axis([0 l_x 0 h_y]);
    title('Streamline of Velocity-Field')
    xlabel('x');
    ylabel('y');
    frame = getframe(gcf);
    im = frame2im(frame);
    imwrite(im,[params,'Streamline_Plot','.jpeg']);
    
    
    figure();
    %quiver(X(1:3:end,1:3:end),Y(1:3:end,1:3:end),u_mat(1:3:end,1:3:end),v_mat(1:3:end,1:3:end));
    quiver(X,Y,u_mat,v_mat);
    hold on
    plot(boundary_bottom(:,1),boundary_bottom(:,2),'-k');
    plot(boundary_top(:,1),boundary_top(:,2),'-k');
    title('Vector Plot of Velocity-Field')
    xlabel('x');
    ylabel('y');
    frame = getframe(gcf);
    im = frame2im(frame);
    imwrite(im,[params,'Vector_plot','.jpeg']);
    
    figure();
    u_avg = zeros(1,N_x);
    u_avg(1:N1_x) = -trapz(Y1(:,1),reshape(u(1:N1),size(X1))); % 1/h1_y missing, because it is simplified in massflux calc
    u_avg(N1_x+1:end) = -trapz(Y2(:,1),reshape(u(N1+N2_y+1:N),size(X2(:,2:end)))); %1/h2_y here
    massflux = density*u_avg;
    plot(massflux);
    axis([ 0 N1_x+N2_x 0 (max(massflux)*1.1)]);
    title('Massflux Plot')
    xlabel('x');
    ylabel('Massflux');
    frame = getframe(gcf);
    im = frame2im(frame);
    imwrite(im,[params,'Massflux','.jpeg']);
    
    figure();
    pcolor(X1,Y1,reshape(divergence(1:N1),size(X1)));
    hold on
    pcolor(X2,Y2,reshape(divergence(N1-N1_y+1:N),size(X2)));
    colorbar;
    plot(boundary_bottom(:,1),boundary_bottom(:,2),'-k');
    plot(boundary_top(:,1),boundary_top(:,2),'-k');
    axis([0 l_x 0 h_y]);
    title('Plot of Divergence Field')
    xlabel('x');
    ylabel('y');
    frame = getframe(gcf);
    im = frame2im(frame);
    imwrite(im,[params,'Divergence_Plot','.jpeg']);
    
    figure()
    semilogy(divergence_norm);
    title('Norm of Div(velocity) for each iteration');
    xlabel('iterations')
    ylabel('Div(velocity)');
    frame = getframe(gcf);
    im = frame2im(frame);
    imwrite(im,[params,'Divergence_Norm_Iterations','.jpeg']);
    
end