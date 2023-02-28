function poly = TriangularBasis(i, der_xi, der_eta, xi, eta)
% To get the basis function or the derivative over the linear parent triangle
% on Page 121. Nodes: p1(-1, -1)   p2(1, -1)   p3(0, 1)
% Input:
%   i: the number of the basis function, 1 ~ 3.
%   der: the derivative required.
%       0 --- the value of basis function.
%       1 --- the 1st derivative of the basis function.
%       Supported choices:
%       (der_xi, der_eta) = 
%       (0, 0)  (1, 0)  (0, 1)  (1, 1)
%       If one of the der > 1, 0 would be output since it is a linear element.
%   xi, eta: the point position in the parent triangle.

if mod(der_xi, 1) ~= 0 || mod(der_eta, 1) ~=0
    disp("TriangularBasis: Please input appropriate 'der' groups.")
    return;
end

if abs(xi) > 1 || abs(eta) > 1
    disp("TriangularBasis: Please input appropriate (xi, eta).")
    return;
end

if i == 1
    if der_xi == 0 && der_eta == 0
        poly = 0.25 * (1 - xi) * (1 - eta);
    elseif der_xi == 1 && der_eta == 0
        poly = -0.25 * (1 - eta);
    elseif der_xi == 0 && der_eta == 1
        poly = -0.25 * (1 - xi);
    elseif der_xi == 1 && der_eta == 1
        poly = 0.25;
    elseif (der_xi > 1 || der_eta > 1) && (der_xi >=0 && der_eta >= 0)
        poly = 0;
    else
        disp("TriangularBasis: Please input appropriate 'der' groups.")
    end

elseif i == 2
    if der_xi == 0 && der_eta == 0
        poly = 0.25 * (1 + xi) * (1 - eta);
    elseif der_xi == 1 && der_eta == 0
        poly = 0.25 * (1 - eta);
    elseif der_xi == 0 && der_eta == 1
        poly = -0.25 * (1 + xi);
    elseif der_xi == 1 && der_eta == 1
        poly = -0.25;
    elseif (der_xi > 1 || der_eta > 1) && (der_xi >=0 && der_eta >= 0)
        poly = 0;
    else
        disp("TriangularBasis: Please input appropriate 'der' groups.")
    end
    
elseif i == 3
    if der_xi == 0 && der_eta == 0
        poly = 0.5 * (1 + eta);
    elseif der_xi == 0 && der_eta == 1
        poly = 0.5;
    elseif (der_xi > 0 || der_eta > 1) && (der_xi >=0 && der_eta >= 0)
        poly = 0;
    else
        disp("TriangularBasis: Please input appropriate 'der' groups.")
    end
else
    disp("TriangularBasis: Please input appropriate 'i'(the node number).")
end

end

