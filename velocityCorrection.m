% description: 2D Poisson equation solution using finite difference dscretization
% author: kfudali
% date: 14.03.2023

clear all
clc

%create domain mesh
Lx = 1;
Ly = 1;
nx = 50;
ny = 50;
n_elements = (nx-1)*(ny-1);

x = [0:(Lx/(nx-1)):Lx];
dx = Lx / nx;

y = [0:(Ly/(ny-1)):Lx];
dy = Ly / ny;

[x,y] = meshgrid(x,y);
points = [x(:), y(:)];

K_e = integrateLaplacian(dx,dy);
K_g = zeros(nx*ny);

%assemblacja
element_nodes_ids = zeros(n_elements, 4);
for i = 1:n_elements
    first_node = i + floor(i/(nx-1));
    element_nodes_ids(i,:) = [first_node, first_node + 1, first_node + nx + 1, first_node + nx];
end
element_nodes_ids((nx-1):(nx-1):end,:) = element_nodes_ids((nx-1):(nx-1):end,:) - 1;


omega = 1;
rho = dx^2;
x_vec = reshape(x,nx*ny,1);
y_vec = reshape(y,nx*ny,1);
W = zeros(nx*ny,1);
w = omega * ones(nx*ny,1);
c = ones(size(W));
for i = 1:n_elements
    K_g(element_nodes_ids(i,:),element_nodes_ids(i,:)) = K_g(element_nodes_ids(i,:),element_nodes_ids(i,:)) + K_e;
    w_local = integrateNeumannBC(x_vec(element_nodes_ids(i,:)), y_vec(element_nodes_ids(i,:)),w(element_nodes_ids(i,:)));
    W(element_nodes_ids(i,:)) = W(element_nodes_ids(i,:)) + w_local;
end
Ap = (K_g + rho/((W'*c)^2) * W*W');

K_g_case = [K_g, ones(nx*ny,1)];
K_g_case = [K_g_case; ones(nx*ny + 1,1)'];
K_g_case(end,end) = 0;

[grad_x, grad_y] = assembleGrad(nx,ny,dx,dy);
L = assembleLaplace(nx,ny,dx,dy);

L_rhs = assembleLaplace(nx+2,ny+2,dx,dy);

%Boundary conditions
top_ids = 1:nx;
bottom_ids = nx*(ny-1) + 1:nx*ny;
left_ids = 1:nx:nx*ny;
right_ids = nx:nx:nx*ny;

boundary_ids = unique([bottom_ids, top_ids, left_ids, right_ids]);
interior_ids = setdiff(1:nx*ny,boundary_ids);

ux = zeros(nx*ny,1);
uy = zeros(nx*ny,1);
p = zeros(nx*ny,1);
phi = zeros(nx*ny,1);

nu = 0.01;
t = 0;
t_step = 1;
dt = 0.01;
t_end = 3;
n_time_steps = t_end/dt;

ux_old = zeros(nx*ny,n_time_steps);
uy_old = zeros(nx*ny,n_time_steps);
p_old = zeros(nx*ny,n_time_steps);

x_vec = reshape(x,nx*ny,1);
y_vec = reshape(y,nx*ny,1);

integrated_rhs_x = zeros(size(ux));
integrated_rhs_y = zeros(size(uy));
 
ux(top_ids) = cos(x_vec(right_ids)*2 * pi) - 1;


x_vec = reshape(y,nx*ny,1);
y_vec = reshape(x,nx*ny,1);
xlabel('y');
ylabel('x');

% 
% [fx,fy,p_analytical] = getManufacturedSolution(x_vec,y_vec,1);

x_vec = reshape(x,nx*ny,1);
y_vec = reshape(y,nx*ny,1);

laplace_rhs_x = L_rhs * ux;
laplace_rhs_x = - laplace_rhs_x(interior_ids);
laplace_rhs_y = L_rhs * uy;
laplace_rhs_y = - laplace_rhs_y(interior_ids);

while t < t_end
    tic();

    %initialize
    if t_step < 4
        ux_old(:,t_step) = ux;
        uy_old(:,t_step) = uy;
        p_old(:,t_step) = p;    
        t_step = t_step + 1;
        t = t +dt;
    else
    
    %step 1
    rhs_x = (1/(2 * dt) * (7*ux_old(:,t_step-1) - 5*ux_old(:,t_step-2) + ux_old(:,t_step-3)));
    rhs_y = (1/(2 * dt) * (7*uy_old(:,t_step-1) - 5*uy_old(:,t_step-2) + uy_old(:,t_step-3)));
    integrated_rhs_x = zeros(size(ux));
    integrated_rhs_y = zeros(size(uy));

    for i = 1 :n_elements
        f_local = integrateUx(x_vec(element_nodes_ids(i,:)), y_vec(element_nodes_ids(i,:)),rhs_x(element_nodes_ids(i,:)));
        integrated_rhs_x(element_nodes_ids(i,:)) = integrated_rhs_x(element_nodes_ids(i,:)) + f_local;

        f_local = integrateUy(x_vec(element_nodes_ids(i,:)), y_vec(element_nodes_ids(i,:)),rhs_y(element_nodes_ids(i,:)));
        integrated_rhs_y(element_nodes_ids(i,:)) = integrated_rhs_y(element_nodes_ids(i,:)) + f_local;
    end

    F = (integrated_rhs_x + integrated_rhs_y);
    [phi, flag] = pcg(Ap, F, 1e-6, 1000, [], [],  phi);

    %step 2
    p = p_old(:, t_step - 1) + phi - nu * (grad_x * ux + grad_y * uy);

    %step 3
    lhs = 3/(2 * dt) * eye((nx-2)*(ny-2)) - nu * L;
    precond = ichol(sparse(lhs));

    Dxp = grad_x * p;
    Dyp = grad_y * p;
    
    rhs_x = -Dxp(interior_ids) + 1/(2 * dt) * (4 * ux_old(interior_ids, t_step - 1) - ux_old(interior_ids, t_step - 2)) + nu * laplace_rhs_x;
    rhs_y = -Dyp(interior_ids)  + 1/(2 * dt) * (4 * uy_old(interior_ids, t_step - 1) - uy_old(interior_ids, t_step - 2)) + nu * laplace_rhs_y;

   [ux(interior_ids) ,flag_x] = pcg(lhs, rhs_x, 1e-9, 200, precond, precond', ux_old(interior_ids, t_step - 1));
   [uy(interior_ids), flag_y] = pcg(lhs, rhs_y, 1e-9, 200, precond, precond', uy_old(interior_ids, t_step - 1));

   if flag_y ~= 0 || flag_x ~0
       disp(['pcg failed with flag = ', num2str(flag)]);
       break;
   end
%     ux_int = lhs\rhs_x;
%     ux_int= lhs\rhs_y;

    ux_old(:,t_step) = ux;
    uy_old(:,t_step) = uy;
    p_old(:,t_step) = p; 
    t_step = t_step + 1;
    t = t + dt;

    vel = (ux.^2 + uy.^2).^(0.5);
    contourf(x,y,reshape(p,nx,ny));

    title('vel')
    colorbar
    drawnow
%     iter_time = toc();
%     disp( ['Iteration time: ', num2str(iter_time)]);

    end
end
% 
% results.u = ux;
% results.v = uy;
% results.x = x;
% results.y = y;
% results.p = p;
% 
% save("results_stokes/results_nu_001.mat","results")