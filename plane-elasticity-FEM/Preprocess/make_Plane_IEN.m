function Plane_IEN = make_Plane_IEN(msh, etype)
% Input:
%   msh: The imported msh info.
%   etype: The type of elements.
%       3 --- Triangle  4 --- Quadrangle
% Output:
%   plane_IEN: The IEN matrix of the triangular elements on the plane.
%       plane_IEN(*, e): The the 'e'th element.
%       plane_IEN(1, e), plane_IEN(2, e), plane_IEN(3, e), plane_IEN(4, e): 
%           The nodes of the 'e'th element.

if etype == 3
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

end