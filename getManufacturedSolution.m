function [fx_vec,fy_vec, p_vec] = getManufacturedSolution(x_vec,y_vec,t_vec)

    Pi = sym(pi, 'd');
    syms u(x,y,t);
    syms v(x,y,t);
    syms p(x,y,t);
    syms fx(x,y,t);
    syms fy(x,y,t);
    syms A(t);
    syms nu;
    
    A(t) = 1;
    nu = 0.1;
    
    u(x,y,t) = A(t) * (1 - cos(2 * Pi * x)) * y * (2 - 3 * y);
    v(x,y,t) = -A(t) * 2 * Pi * sin(2 * Pi * x) * y^2 * (1 - y);
    p(x,y,t) = nu*((x-0.5)^2 + (y-0.5)^2);
%     (u*diff(u,x) + v*diff(u,y))
%     (u*diff(v,x) + v*diff(v,y))
    fx(x,y,t) = diff(u,t) + diff(p,x) - nu * laplacian(u,[x,y]);
    fy(x,y,t) = diff(v,t) + diff(p,y) - nu * laplacian(v,[x,y]);
    
    fx_vec = double(fx(x_vec,y_vec,t_vec));
    fy_vec = double(fy(x_vec, y_vec, t_vec));
    p_vec = double(p(x_vec,y_vec,t_vec));

end