function [Diri_Nodes, Neum_IEN, plane_IEN] = Preprocess(file_name, num_Diri, num_Neum)
% The preprocessor for plane-strain problem
% Input:
%   file_name: The Gmsh file with '.msh' format (exactly msh2)in which
%               the physical groups should be set as the following order:
%               1 ~ m --- the Dirichlet boundary different with each other
%               m+1 ~ m+n --- the Neumann boundary different with each other
%               m+n+1 --- the plane surface
%   num_Diri: The number of Dirichlet boundaries i.e the previous 'm'.
%   num_Neum: The number of Neumann boundaries i.e the previous 'n'.
% Output:
%   Diri_nodes: The nodes infomation of Dirichlet boundaries (as a cell).
%       Diri_nodes{i}(j): the 'j'the node of the 'i'th Dirichlet boundary
%   Neum_IEN: The IEN matrix of the linear elements on Neumann boundaries
%               (as a cell).
%       Neum_IEN{i}: The IEN matrix of the 'i'th Neumann boundary
%       Neum_IEN{i}(*, e): the 'e'th element
%       Neum_IEN{i}(1, e), Neum_IEN{i}(2, e): the nodes of the 'e'th element
%       Neum_IEN{i}(3, e): the globel number of the 'e'th element
%   plane_IEN: The IEN matrix of the triangular elements on the plane.
%       plane_IEN(*, e): the the 'e'th element
%       plane_IEN(1, e), plane_IEN(2, e), plane_IEN(3, e): 
%           the nodes of the 'e'th element
%       plane_IEN(4, e): the globel number of the 'e'th element

msh = load_gmsh2(file_name);

% Part 1 ---------- To obtain Diri_Nodes ----------
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

% Part 2 ---------- To obtain Neum_IEN ----------
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
    N_IEN = zeros(3, nb_lineEle);
    temp = 1;
    for jj = 1 : msh.nbLines
        if msh.LINES(jj, 3) == ii
            N_IEN(1, temp) = msh.LINES(jj, 1);
            N_IEN(2, temp) = msh.LINES(jj, 2);
            N_IEN(3, temp) = jj;
            temp = temp + 1;
        end
    end
    
% Put it into the cell 
    Neum_IEN{ii - num_Diri} = N_IEN; 
end

% Part 3 ---------- To obtain plane_IEN ----------
plane_IEN = zeros(4, msh.nbTriangles);

% Construct IEN
temp = 1;
for kk = msh.nbLines + 1 : msh.nbElm
    plane_IEN(1, temp) = msh.ELE_NODES(kk, 1);
    plane_IEN(2, temp) = msh.ELE_NODES(kk, 2);
    plane_IEN(3, temp) = msh.ELE_NODES(kk, 3);
    plane_IEN(4, temp) = kk;
    temp = temp + 1;
end

end

