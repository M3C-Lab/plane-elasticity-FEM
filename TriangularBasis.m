function poly = TriangularBasis(i, der_x, der_y, x, y)
% To get the basis function or the derivative over a triangular element.
% The value of these basis functions are equivalent to the AREA COORDINATES.
% The nodes of default parent triangle: p1(0, 0)    p2(1, 0)    p3(0, 1).
%   You can choose another set of nodes, but this script can only be used
%   for LINEAR elements, and the nodes should be arranged COUNTERCLOCKWISE.
% Input:
%   i: the number of the basis function, 1 ~ 3.
%   der: the derivative required.
%       0 --- the value of basis function.
%       1 --- the 1st derivative of the basis function.
%       Supported choices:
%       (der_x, der_y) = 
%       (0, 0)  (1, 0)  (0, 1)
%       If one of the der > 1, 0 would be output since it is a linear element.
%   x, y: The points position in parent triangle.
% Output:
%   poly: The value of basis function or derivative.
% 关于直角坐标与面积坐标的转换关系，可参看清华大学出版社《有限单元法》第106页及第57页。

% The nodes of the parent triangular element (adjustable).
p1 = [0, 0]';
p2 = [1, 0]';
p3 = [0, 1]';

% The formulation in Page 57,58 of 《有限单元法》.
a1 = p2(1) * p3(2) - p3(1) * p2(2);
b1 = p2(2) - p3(2);
c1 = -p2(1) + p3(1);
a2 = p3(1) * p1(2) - p1(1) * p3(2);
b2 = p3(2) - p1(2);
c2 = -p3(1) + p1(1);
a3 = p1(1) * p2(2) - p2(1) * p1(2);
b3 = p1(2) - p2(2);
c3 = -p1(1) + p2(1);
coeffient_matrix = [1, p1(1), p1(2); 1, p2(1), p2(2); 1, p3(1), p3(2)];
A = 0.5 * det(coeffient_matrix);

if i == 1
    if der_x == 0 && der_y == 0
        poly = 0.5 * (a1 + b1 * x + c1 * y) / A;
    elseif der_x == 1 && der_y == 0
        poly = 0.5 * b1 / A;
    elseif der_x == 0 && der_y == 1
        poly = 0.5 * c1 / A;
    else
        disp("TriangularBasis: Please input appropriate 'der' groups.")
    end
    
elseif i == 2
    if der_x == 0 && der_y == 0
        poly = 0.5 * (a2 + b2 * x + c2 * y) / A;
    elseif der_x == 1 && der_y == 0
        poly = 0.5 * b2 / A;
    elseif der_x == 0 && der_y == 1
        poly = 0.5 * c2 / A;
    else
        disp("TriangularBasis: Please input appropriate 'der' groups.")
    end
    
elseif i == 3
    if der_x == 0 && der_y == 0
        poly = 0.5 * (a3 + b3 * x + c3 * y) / A;
    elseif der_x == 1 && der_y == 0
        poly = 0.5 * b3 / A;
    elseif der_x == 0 && der_y == 1
        poly = 0.5 * c3 / A;
    else
        disp("TriangularBasis: Please input appropriate 'der' groups.")
    end
    
else
    disp("TriangularBasis: Please input appropriate 'i'(the node number).")
end

end

