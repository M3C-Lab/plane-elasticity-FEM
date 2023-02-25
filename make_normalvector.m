function normalvector = make_normalvector(Neum_IEN, POS)
% To get the external normal vectors of Neumann boundaries
% Input:
%   Neum_IEN: The IEN matrix of the line elements on Neumann boundaries
%               (as a cell).
%       Neum_IEN{i}: The IEN matrix of the 'i'th Neumann boundary
%       Neum_IEN{i}(*, e): The 'e'th element.
%       Neum_IEN{i}(1, e), Neum_IEN{i}(2, e): The nodes of the 'e'th element.
%       Neum_IEN{i}(3, e): Which triangular element the 'e'th line element 
%                          is located at.
%       Neum_IEN{i}(4, e): The other node of the triangular element.
%   POS: the coordinates [x,y,z] (positions) of the nodes
% Output:
%   normalvector: The external normal vectors with respect to Neum_IEN.
%               (as a cell).
%       normalvector{i}: The normal vectors of the 'i'th Neumann boundary
%       normalvector{i}(*, e): The 'e'th element.
%       normalvector{i}(1, e): The x-direction component of the vector
%       normalvector{i}(2, e): The y-direction component of the vector
%       normalvector{i}(3, e): Which triangular element the 'e'th line element 
%                          is located at, checked with Neum_IEN{i}(3, e)

normalvector = cell(1, length(Neum_IEN));
Rotation_matrix = [0, -1;
                   1,  0];

for ii = 1 : length(Neum_IEN)
    nLineEle = size(Neum_IEN{ii}, 2);
    nvs = zeros(3, nLineEle);
    for kk = 1 : nLineEle
        nvs(3, kk) = Neum_IEN{ii}(3, kk);
        node1 = Neum_IEN{ii}(1, kk);
        node2 = Neum_IEN{ii}(2, kk);
        x_comp = POS(node1, 1) - POS(node2, 1);
        y_comp = POS(node1, 2) - POS(node2, 2);
        arc = sqrt(x_comp ^2 + y_comp ^2);
        tv = [x_comp/arc, y_comp/arc]';
        nv = Rotation_matrix * tv;
        
        node3 = Neum_IEN{ii}(4, kk);
        test_v_x = POS(node1, 1) - POS(node3, 1);
        test_v_y = POS(node1, 2) - POS(node3, 2);
        if test_v_x * nv(1) + test_v_y * nv(2) < 0
            nv = -nv;
        end
        
        nvs(1, kk) = nv(1);
        nvs(2, kk) = nv(2);
    end
    
    normalvector{ii} = nvs;
    
end

end

