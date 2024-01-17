function [Ix, Iy] = getShapeFunctionsRHS(dx,dy)

    syms n1(x,y,a,b);
    syms n2(x,y,a,b);
    syms n3(x,y,a,b);
    syms n4(x,y,a,b);
    
    n1(x,y) = (1 - x/a) * (1 - y/b);
    n2(x,y) = x/a * (1- y/b);
    n3(x,y) = x*y / (a*b);
    n4(x,y) = (1 - x/a) * y/b;
    
    N = {n1;n2;n3;n4};
    
    Nx = diff(N,x);
    Ny = diff(N,y);
    
    Nxx = [Nx * N(1), Nx * N(2), Nx * N(3), Nx * N(4)];
    Nyy = [Ny * N(1), Ny * N(2), Ny * N(3), Ny * N(4)];

    Ix = int(int(Nxx,x,[0,a]),y,[0,b]);
    Iy = int(int(Nyy,x,[0,a]),y,[0,b]);

    Ix = subs(Ix,[a,b],[dx,dy]);
    Iy = subs(Iy,[a,b],[dx,dy]);

end