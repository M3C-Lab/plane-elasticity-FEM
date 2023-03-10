function poly = TriangularBasis(i, der_x, der_y, x, y, phys2rst)
% To get the basis function or the derivative over a triangular element.
% The value of these basis functions are equivalent to the AREA COORDINATES.
% The nodes of default parent triangle.
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
%   x, y: The points position in the triangular element.
%   phys2rst: The mapping matrix given by Mapping_p2rst.m .
% Output:
%   poly: The value of basis function or derivative.
% 关于直角坐标与面积坐标的转换关系，可参看清华大学出版社《有限单元法》第57页及第106页。

if i == 1
    if der_x == 0 && der_y == 0
        poly = phys2rst(1, 1) * x + phys2rst(1, 2) * y + phys2rst(1, 3);
    elseif der_x == 1 && der_y == 0
        poly = phys2rst(1, 1);
    elseif der_x == 0 && der_y == 1
        poly = phys2rst(1, 2);
    else
        disp("TriangularBasis: Please input appropriate 'der' groups.")
    end
    
elseif i == 2
    if der_x == 0 && der_y == 0
        poly = phys2rst(2, 1) * x + phys2rst(2, 2) * y + phys2rst(2, 3);
    elseif der_x == 1 && der_y == 0
        poly = phys2rst(2, 1);
    elseif der_x == 0 && der_y == 1
        poly = phys2rst(2, 2);
    else
        disp("TriangularBasis: Please input appropriate 'der' groups.")
    end
    
elseif i == 3
    if der_x == 0 && der_y == 0
        poly = phys2rst(3, 1) * x + phys2rst(3, 2) * y + phys2rst(3, 3);
    elseif der_x == 1 && der_y == 0
        poly = phys2rst(3, 1);
    elseif der_x == 0 && der_y == 1
        poly = phys2rst(3, 2);
    else
        disp("TriangularBasis: Please input appropriate 'der' groups.")
    end
    
else
    disp("TriangularBasis: Please input appropriate 'i'(the node number).")
end

end

