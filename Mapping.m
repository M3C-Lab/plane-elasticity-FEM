function [para2phys,phys2para] = Mapping(n1, n2, n3, p1, p2, p3)
% To get mapping matrice between the physical space and parametric space.
% Input:
%   n1, n2, n3: The points in the physical space.
%   p1, p2, p3: The points in the parametric space.
%       Each of them should be a array: [x, y] or [xi, eta].
%       If three points are on a line (cannot create a triangle) in a space,
%       this formulation would fail (Matrix A or B is singular).
% Output:
%   para2phys: The mapping matrix from parametric space to physical space.
%       x = para2phys(1,1) * xi + para2phys(1,2) * eta + para2phys(1,3)
%       y = para2phys(2,1) * xi + para2phys(2,2) * eta + para2phys(2,3)
%   phys2para: The mapping matrix from physical space to parametric space.
%       xi = phys2para(1,1) * x + phys2para(1,2) * y + phys2para(1,3)
%       eta = phys2para(2,1) * x + phys2para(2,2) * y + phys2para(2,3)

% Hence, dx_dxi = para2phys(1,1), dx_deta = para2phys(1,2),
%        dy_dxi = para2phys(2,1), dy_deta = para2phys(2,2),
%       and dxi_dx, dxi_dy, deta_dx, deta_dy are what we need.

A = [n1(1), n2(1), n3(1);
     n1(2), n2(2), n3(2);
     1, 1, 1];

B = [p1(1), p2(1), p3(1);
     p1(2), p2(2), p3(2);
     1, 1, 1];
 
% A = para2phys * B
para2phys = A / B;

% B = phys2para * A;
phys2para = B / A;

end

