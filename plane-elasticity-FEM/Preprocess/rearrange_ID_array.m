function New_ID_array = rearrange_ID_array(ID_array, Diri_Node_array, direction)
% To insert Dirichlet boundary conditions to the ID_array
% Input:
%   ID_array: The original ID_array.
%   Diri_Node_array: The nodes on the Dirichlet boundary as a vector. 
%   direction: It should be input as 'x' or 'y', but not both of them.
% Output:
%   New_ID_array: The rearranged ID_array, and the info of Dirichlet 
%       boundary nodes would be set to 0 according to the degree of freedom.

nbNod = size(ID_array, 2);
nb_Diri_Node = length(Diri_Node_array);

if direction == 'x'
    for ii = 1 : nb_Diri_Node
        if ID_array(1, Diri_Node_array(ii)) ~= 0
            ID_array(1, Diri_Node_array(ii)) = 0;
            if ID_array(2, Diri_Node_array(ii)) ~= 0
                ID_array(2, Diri_Node_array(ii)) = ...
                    ID_array(2, Diri_Node_array(ii)) - 1;
            end
            for jj = Diri_Node_array(ii) + 1 : nbNod
                for kk = 1 : 2
                    if ID_array(kk, jj) ~= 0
                        ID_array(kk, jj) = ID_array(kk, jj) - 1;
                    end
                end
            end
        end
    end
    
elseif direction == 'y'
    for ii = 1 : nb_Diri_Node
        if ID_array(2, Diri_Node_array(ii)) ~= 0
            ID_array(2, Diri_Node_array(ii)) = 0;
            for jj = Diri_Node_array(ii) + 1 : nbNod
                for kk = 1 : 2
                    if ID_array(kk, jj) ~= 0
                        ID_array(kk, jj) = ID_array(kk, jj) - 1;
                    end
                end
            end
        end
    end
    
else
    disp("Please input a direction ('x' or 'y') as the third parameter.")
    return;
end

New_ID_array = ID_array;
end
