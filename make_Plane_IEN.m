function Plane_IEN = make_Plane_IEN(msh)
% Input:
%   msh: The imported msh info.
% Output:
%   plane_IEN: The IEN matrix of the triangular elements on the plane.
%       plane_IEN(*, e): The the 'e'th element.
%       plane_IEN(1, e), plane_IEN(2, e), plane_IEN(3, e): 
%           The nodes of the 'e'th element.

Plane_IEN = zeros(3, msh.nbTriangles);

% Construct IEN
temp = 1;
for jj = msh.nbLines + 1 : msh.nbElm
    Plane_IEN(1, temp) = msh.ELE_NODES(jj, 1);
    Plane_IEN(2, temp) = msh.ELE_NODES(jj, 2);
    Plane_IEN(3, temp) = msh.ELE_NODES(jj, 3);
    temp = temp + 1;
end

end