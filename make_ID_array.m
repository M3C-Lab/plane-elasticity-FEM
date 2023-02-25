function Initial_IDarray = make_ID_array(nbNod)
% To create a unmodified ID_array with 2 degree of freedom
% Input:
%   nbNod: The number of nodes.
% Output:
%   Initial_IDarray: The unmodified ID_array
Initial_IDarray = zeros(2, nbNod);

for ii = 1 : nbNod
    Initial_IDarray(1, ii) = 2 * ii - 1;
    Initial_IDarray(2, ii) = 2 * ii;
end

end

