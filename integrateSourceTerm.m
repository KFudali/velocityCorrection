function [f_integrated] = integrateSourceTerm(x,y,f)
    
    q_weights = [1.000000 1.000000 1.000000 1.000000];
    q_points = [-0.5774   -0.5774
               -0.5774    0.5774
                0.5774   -0.5774
                0.5774    0.5774];

    % The domain vertices
    xi = q_points(:,1);
    eta = q_points(:,2);

    N1 = @(xi,eta)(1-xi).*(1-eta)/4;
    N2 = @(xi,eta)(1+xi).*(1-eta)/4;
    N3 = @(xi,eta)(1+xi).*(1+eta)/4;
    N4 = @(xi,eta)(1-xi).*(1+eta)/4;

    evalN = [N1(xi,eta), N2(xi,eta), N3(xi,eta), N4(xi,eta)];

    qPointsLocal = evalN * [x,y];

    J = @(xi,eta) [-(1-eta)/4, (1-eta)/4, (1+eta)/4, -(1+eta)/4; ...
                     -(1-xi)/4,-(1+xi)/4, (1+xi)/4,  (1-xi)/4];
    for i=1:size(xi,1)
        evalDetJacb(i) = abs(det(J(xi(i),eta(i))*[x,y]));
    end

    x0 = x(1);
    y0 = y(1);
    dx = abs(x(1) - x(4));
    dy = abs(y(1) - y(2));

    n1 = @(x,y)(1 - (x-x0)/dx).*(1 - (y-y0)/dy);
    n2 = @(x,y)(x-x0)/dx.*(1- (y-y0)/dy);
    n3 = @(x,y)(x-x0).*(y-y0) / (dx*dy);
    n4 = @(x,y)(1 - (x-x0)/dx).* (y-y0)/dy;

    eval_n = f.*[n1(qPointsLocal(:,1),qPointsLocal(:,2)), n2(qPointsLocal(:,1),qPointsLocal(:,2)), n3(qPointsLocal(:,1),qPointsLocal(:,2)), n4(qPointsLocal(:,1),qPointsLocal(:,2))];
    % Finally, apply Gauss formula
    f_integrated =  zeros(4,1);
    for i = 1:4
         f_integrated = f_integrated + q_weights(i) * eval_n(i,:)' * evalDetJacb(i);
    end
end