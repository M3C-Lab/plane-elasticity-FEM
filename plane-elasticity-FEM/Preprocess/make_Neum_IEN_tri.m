function Neum_IEN = make_Neum_IEN_tri(msh, num_Diri, num_Neum, Plane_IEN)
% To create IEN arrays of line elements on Neumann Boundaries
% when the plane elements are triangular.
% Input:
%   msh: The imported msh info.
%       The physical groups should be set with the following order:
%       1 ~ m --- the Dirichlet boundary different with each other
%       m+1 ~ m+n --- the Neumann boundary different with each other
%       m+n+1 --- the plane surface
%   num_Diri: The number of Dirichlet boundaries i.e the previous 'm'.
%   num_Neum: The number of Neumann boundaries i.e the previous 'n'.
%   Plane_IEN: The constructed IEN of triangular elements.
% Output:
%   Neum_IEN: The IEN matrix of the line elements on Neumann boundaries
%               (as a cell).
%       Neum_IEN{i}: The IEN matrix of the 'i'th Neumann boundary
%       Neum_IEN{i}(*, e): The 'e'th element.
%       Neum_IEN{i}(1, e), Neum_IEN{i}(2, e): The nodes of the 'e'th element.
%       Neum_IEN{i}(3, e): Which triangular element the 'e'th line element 
%                          is located at.
%       Neum_IEN{i}(4, e): The other node of the triangular element.

Neum_IEN = cell(1, num_Neum);
% Search for Neumann boundaries
for ii = num_Diri + 1 : num_Diri + num_Neum
    nb_lineEle = 0;
% Search for element number on a Dirichlet boundary
    for jj = 1 : msh.nbLines
        if msh.LINES(jj, 3) == ii
            nb_lineEle = nb_lineEle + 1;
        end      
    end
    
% Construct IEN
    N_IEN = zeros(4, nb_lineEle);
    temp = 1;
    for jj = 1 : msh.nbLines
        if msh.LINES(jj, 3) == ii
            % The nodes of the line element
            N_IEN(1, temp) = msh.LINES(jj, 1);
            N_IEN(2, temp) = msh.LINES(jj, 2);
            
            % Search the location of each line element
            node_sum = N_IEN(1, temp) + N_IEN(2, temp);
            node_diff = abs(N_IEN(1, temp) - N_IEN(2, temp));
            % Compare 'node1 + node2' and '|node1 - node2|' respectively
            for kk = 1 : msh.nbTriangles       
                for ll = 1 : 3
                    choices = [1, 2, 3];
                    choices(ll) = [ ];
                    sum_test = ...
                        Plane_IEN(choices(1), kk) + Plane_IEN(choices(2), kk);
                    if sum_test == node_sum
                        diff_test = ...
                            abs(Plane_IEN(choices(1), kk) - Plane_IEN(choices(2), kk));
                        if diff_test == node_diff
                            N_IEN(3, temp) = kk;
                            % The 'temp'th line element is the boundary of
                            % the 'kk'th triangular element
                            N_IEN(4, temp) = Plane_IEN(ll, kk);
                            % 'Plane_IEN(ll, kk)' is the other node
                            break;
                        end
                    end
                end
                
                if N_IEN(3, temp) ~= 0
                    break;
                end
                
            end
            temp = temp + 1;
        end       
    end
    
% Put it into the cell 
    Neum_IEN{ii - num_Diri} = N_IEN; 
end

end