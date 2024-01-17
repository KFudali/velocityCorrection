% description: 2D Poisson equation solution using finite element dscretization
% author: kfudali
% date: 14.03.2023

clear all
clc

%create domain mesh
Lx = 1;
Ly = 1;
nx = 20;
ny = 20;
n_elements = (nx-1)*(ny-1);

x = [0:(Lx/(nx-1)):Lx];
dx = Lx / nx;
 
y = [0:(Ly/(ny-1)):Ly];
dy = Ly / ny;

[x,y] = meshgrid(x,y);
points = [x(:), y(:)];

K_e = integrateLaplacian(dx,dy);
K_g = zeros(nx*ny);

%assemblacja
element_nodes_ids = zeros(n_elements, 4);

for i = 1:n_elements
    first_node = i + floor(i/(nx-1));
    element_nodes_ids(i,:) = [first_node, first_node + 1,  first_node + nx + 1, first_node + nx];
end
element_nodes_ids((nx-1):(nx-1):end,:) = element_nodes_ids((nx-1):(nx-1):end,:) - 1;

for i = 1:n_elements
    K_g(element_nodes_ids(i,:),element_nodes_ids(i,:)) = K_g(element_nodes_ids(i,:),element_nodes_ids(i,:)) + K_e;
end


%Boundary conditions
bottom_ids = 1:nx;
top_ids = nx*(ny-1) + 1:nx*ny;
left_ids = 1:nx:nx*ny;
right_ids = nx:nx:nx*ny;

boundary_ids = unique([bottom_ids, top_ids, left_ids, right_ids]);
interior_ids = setdiff(1:nx*ny,boundary_ids);


%% Dirichlet
u = zeros(nx*ny,1);
u(top_ids) = 1;
rhs = -K_g * u;

K_g(boundary_ids,:) = [];
K_g(:,boundary_ids) = [];

rhs(boundary_ids) = [];
u_solved = K_g\rhs;
u(interior_ids) = u_solved;

surf(reshape(u,nx,ny))

%% Dirichlet BC + source term

f = - 2*pi^2 * sin(pi * x).* sin(pi * y);
f = reshape(f,nx*ny,1);

F = zeros(nx*ny,1);
x_vec = reshape(x,nx*ny,1);
y_vec = reshape(y,nx*ny,1);

for i = 1 :n_elements
    f_local = integrateSourceTerm(x_vec(element_nodes_ids(i,:)), y_vec(element_nodes_ids(i,:)),f(element_nodes_ids(i,:)));
    F(element_nodes_ids(i,:)) = F(element_nodes_ids(i,:)) + f_local;
end

K_g_case = K_g;

u = zeros(nx*ny,1);
rhs =K_g_case* u - F;

K_g_case(boundary_ids,:) = [];
K_g_case(:,boundary_ids) = [];

rhs(boundary_ids) = [];
u_solved = K_g_case\rhs;
u(interior_ids) = u_solved;

analytic_solution = sin(pi * x).* sin(pi * y);
contourf(x,y,reshape(u,nx,ny) - analytic_solution)
contourf(x,y,reshape(u,nx,ny))
colorbar
max(u)
%% Mixed boundary condition + source term
x_mix = x/2;
y_mix = y/2;

f = - 2*pi^2 * sin(pi * x_mix).* sin(pi * y_mix);
f = reshape(f,nx*ny,1);


u = zeros(nx*ny,1);
g = zeros(nx*ny,1);

F = zeros(nx*ny,1);
x_vec = reshape(x_mix,nx*ny,1);
y_vec = reshape(y_mix,nx*ny,1);

for i = 1 :n_elements
    f_local = integrateSourceTerm(x_vec(element_nodes_ids(i,:)), y_vec(element_nodes_ids(i,:)),f(element_nodes_ids(i,:)));
    g_local = integrateNeumannBC(x_vec(element_nodes_ids(i,:)), y_vec(element_nodes_ids(i,:)),g(element_nodes_ids(i,:)));
    F(element_nodes_ids(i,:)) = F(element_nodes_ids(i,:)) + f_local + g_local;
end
K_g_case = K_g;

rhs = K_g_case* u - F;

dir_bc_ids = [bottom_ids, left_ids];

K_g_case(dir_bc_ids,:) = [];
K_g_case(:,dir_bc_ids) = [];

rhs(dir_bc_ids) = [];
u_solved = K_g_case\rhs;

mix_interior_ids = setdiff(1:nx*ny,dir_bc_ids);

u(mix_interior_ids) = u_solved;

analytic_solution = sin(pi * x_mix).* sin(pi * y_mix);
contourf(x,y,reshape(u,nx,ny) - analytic_solution)
% contourf(x_mix,y_mix,reshape(u,nx,ny))
colorbar
max(u)


%% Non-zero Neumann boundary condition + source term

f = 2 * ones(nx*ny,1);
% f = reshape(f,nx*ny,1)

u = zeros(nx*ny,1);
g = zeros(nx*ny,1);

g(top_ids) = 0;
g(right_ids) = 0;
g(left_ids) = 0;
g(bottom_ids) = 0;


F = zeros(nx*ny,1);
x_vec = reshape(x, nx*ny,1);
y_vec = reshape(y, nx*ny,1);

for i = 1 :n_elements
    f_local = integrateSourceTerm(x_vec(element_nodes_ids(i,:)), y_vec(element_nodes_ids(i,:)),f(element_nodes_ids(i,:)));
    g_local = integrateNeumannBC(x_vec(element_nodes_ids(i,:)), y_vec(element_nodes_ids(i,:)),g(element_nodes_ids(i,:)));
    F(element_nodes_ids(i,:)) = F(element_nodes_ids(i,:)) + f_local + g_local;
end
K_g_case = K_g;

rhs = K_g_case* u - F;

rhs = [rhs;0];
u = [u;0];
K_g_case = [K_g_case, ones(nx*ny,1)];
K_g_case = [K_g_case; ones(nx*ny + 1,1)'];
K_g_case(end,end) = 0;


u = K_g_case\rhs;
alfa = u(end);
u(end) = [];


% analytic_solution = sin(pi * x).* sin(pi * y);
% contourf(x,y,reshape(u,nx,ny) - analytic_solution)
contourf(x,y,reshape(u,nx,ny))
surf(x,y,reshape(u,nx,ny))
colorbar
max(u)


%% Non-zero Neumann boundary condition + source term

f = - 2*pi^2 * sin(pi * x).* sin(pi * y);
f = reshape(f,nx*ny,1);

w = 1 * ones(nx*ny,1);

u = zeros(nx*ny,1);

F = zeros(nx*ny,1);
W = zeros(nx*ny,1);
x_vec = reshape(x, nx*ny,1);
y_vec = reshape(y, nx*ny,1);

for i = 1 :n_elements
    f_local = integrateSourceTerm(x_vec(element_nodes_ids(i,:)), y_vec(element_nodes_ids(i,:)),f(element_nodes_ids(i,:)));
    w_local = integrateNeumannBC(x_vec(element_nodes_ids(i,:)), y_vec(element_nodes_ids(i,:)),w(element_nodes_ids(i,:)));
    F(element_nodes_ids(i,:)) = F(element_nodes_ids(i,:)) + f_local;
    W(element_nodes_ids(i,:)) = W(element_nodes_ids(i,:)) + w_local;
end

c =ones(size(W));

Ap = (K_g + dx^2/((W'*c)^2) * W*W');
rhs = - F;

u = pcg(Ap, rhs, 1e-6,200);

K_g_case = K_g;



rhs = [rhs;0];
u = [u;0];
K_g_case = [K_g_case, ones(nx*ny,1)];
K_g_case = [K_g_case; ones(nx*ny + 1,1)'];
K_g_case(end,end) = 0;


u = K_g_case\rhs;
alfa = u(end);
u(end) = [];


% analytic_solution = sin(pi * x).* sin(pi * y);
% contourf(x,y,reshape(u,nx,ny) - analytic_solution)
contourf(x,y,reshape(u,nx,ny))
surf(x,y,reshape(u,nx,ny))
colorbar
max(u)
