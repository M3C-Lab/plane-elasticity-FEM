function Diri_Nodes = make_Diri_Nodes(msh, num_Diri)
% Input:
%   msh: The imported msh info.
%       The physical groups should be set with the following order:
%       1 ~ m --- the Dirichlet boundary different with each other
%       m+1 ~ m+n --- the Neumann boundary different with each other
%       m+n+1 --- the plane surface
%   num_Diri: The number of Dirichlet boundaries i.e the previous 'm'.
% Output:
%   Diri_nodes: The nodes infomation of Dirichlet boundaries (as a cell).
%       Diri_nodes{i}(j): The 'j'the node of the 'i'th Dirichlet boundary.

Diri_Nodes = cell(1, num_Diri);
% Search for Dirichlet boundaries
for ii = 1 : num_Diri
    nb_lineEle = 0;
% Search for element number on a Dirichlet boundary
    for jj = 1 : msh.nbLines
        if msh.LINES(jj, 3) == ii
            nb_lineEle = nb_lineEle + 1;
        end      
    end
    
% Search for element info
    lineEle = zeros(1, nb_lineEle);
    temp = 1;
    for jj = 1 : msh.nbLines
        if msh.LINES(jj, 3) == ii
            lineEle(temp) = jj;
            temp = temp + 1;
        end
    end
% Search for node info
    temp = 1;
    lineNode = zeros(1, 2* nb_lineEle);
    for kk = 1 : nb_lineEle
        for jj = 1 : 2
            lineNode(temp) = msh.LINES(lineEle(kk), jj);
            temp = temp + 1;
        end
    end
    lineNode = unique(lineNode);
    
% Put it in to the cell
    Diri_Nodes{ii} = lineNode;
end

end