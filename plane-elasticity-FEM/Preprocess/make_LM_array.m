function LM_array = make_LM_array(Plane_IEN, ID_array)
% To create a LM_array
% Input:
%   Plane_IEN: The IEN array of triangular elements.
%   ID_array: The ID array modified with Dirichlet boundary conditions.
% Output:
%   LM_array: The LM array.

nEle = size(Plane_IEN, 2);
LM_array = zeros(6, nEle);
for ii = 1 : nEle
    for jj = 1 : 3
        for kk = 1 : 2
            LM_array(2*(jj-1)+kk , ii) = ID_array(kk, Plane_IEN(jj, ii));
        end
    end
end

end

