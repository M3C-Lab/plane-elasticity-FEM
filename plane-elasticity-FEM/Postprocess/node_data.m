function [X, Y, XX, YY, strain_node, stress_node] = ...
    node_data(msh, Plane_IEN, element_type, trial_solution)
% To compute strain and stress on each node of the mesh
% Input:
%   msh: the Gmsh info imported.
%   Plane_IEN: The IEN array of plane elment.
%   element_type: 3 --- triangle
%   trial_solution: The solution data.
% Output:
%   X: Orignal X-coordinate of each node.
%   Y: Orignal Y-coordinate of each node.
%   XX: X after the deformation.
%   YY: Y after the deformation.
%   strain_node: [汍_xx, 汍_yy, 污_xy] (averaged for elements)
%       strain_node(i, *): The 'i'th node.
%   stress_node: [考_xx, 考_yy, 而_xy] (averaged for elements)
%       stress_node(i, *): The 'i'th node.

global D;

X = msh.POS(:, 1);
Y = msh.POS(:, 2); 
XX = zeros(msh.nbNod, 1);
YY = zeros(msh.nbNod, 1);
for nn = 1 : msh.nbNod
    XX(nn) = X(nn) + trial_solution(2*nn - 1); 
    YY(nn) = Y(nn) + trial_solution(2*nn);
end

strain_node = zeros(3, msh.nbNod); % [汍_xx, 汍_yy, 污_xy]' of the node.
stress_node = zeros(3, msh.nbNod); % [考_xx, 考_yy, 而_xy]' of the node.

counter_node = zeros(1, msh.nbNod);
% To count how many elements contain this node.

for ee = 1 : msh.nbElm - msh.nbLines
    p1 = [msh.POS(Plane_IEN(1, ee), 1), msh.POS(Plane_IEN(1, ee), 2)]';
    p2 = [msh.POS(Plane_IEN(2, ee), 1), msh.POS(Plane_IEN(2, ee), 2)]';
    p3 = [msh.POS(Plane_IEN(3, ee), 1), msh.POS(Plane_IEN(3, ee), 2)]';
    
    [J, phys2rst] = Mapping_p2rst(p1, p2, p3);
    
    d_ele = [trial_solution(2*Plane_IEN(1, ee) - 1);
             trial_solution(2*Plane_IEN(1, ee));
             trial_solution(2*Plane_IEN(2, ee) - 1);
             trial_solution(2*Plane_IEN(2, ee));
             trial_solution(2*Plane_IEN(3, ee) - 1);
             trial_solution(2*Plane_IEN(3, ee))];
         
    if element_type == 3
        for nn = 1 : 3
            node_x = msh.POS(Plane_IEN(nn, ee), 1);
            node_y = msh.POS(Plane_IEN(nn, ee), 2);
        
            dN1_dx = TriangularBasis(1, 1, 0, node_x, node_y, phys2rst);
            dN2_dx = TriangularBasis(2, 1, 0, node_x, node_y, phys2rst);
            dN3_dx = TriangularBasis(3, 1, 0, node_x, node_y, phys2rst);
            dN1_dy = TriangularBasis(1, 0, 1, node_x, node_y, phys2rst);
            dN2_dy = TriangularBasis(2, 0, 1, node_x, node_y, phys2rst);
            dN3_dy = TriangularBasis(3, 0, 1, node_x, node_y, phys2rst); 
            B_node = [dN1_dx, 0, dN2_dx, 0, dN3_dx, 0;
                  0, dN1_dy, 0, dN2_dy, 0, dN2_dy;
                  dN1_dy, dN1_dx, dN2_dy, dN2_dx, dN3_dy, dN3_dx];
              
            epsilon_node = B_node * d_ele;
        
            strain_node(1, Plane_IEN(nn, ee)) = ...
                strain_node(1, Plane_IEN(nn, ee)) + epsilon_node(1);
            strain_node(2, Plane_IEN(nn, ee)) = ...
                strain_node(2, Plane_IEN(nn, ee)) + epsilon_node(2);
            strain_node(3, Plane_IEN(nn, ee)) = ...
                strain_node(3, Plane_IEN(nn, ee)) + epsilon_node(3);
        
            sigma_node = D * epsilon_node;
        
            stress_node(1, Plane_IEN(nn, ee)) = ...
                stress_node(1, Plane_IEN(nn, ee)) + sigma_node(1);
            stress_node(2, Plane_IEN(nn, ee)) = ...
                stress_node(2, Plane_IEN(nn, ee)) + sigma_node(2);
            stress_node(3, Plane_IEN(nn, ee)) = ...
                stress_node(3, Plane_IEN(nn, ee)) + sigma_node(3);
        
            counter_node(Plane_IEN(nn, ee)) = ...
                counter_node(Plane_IEN(nn, ee)) + 1;
        end        
    end   
end

for nn = 1 : msh.nbNod
    strain_node(:, nn) = strain_node(:, nn) ./ counter_node(nn);    
    stress_node(:, nn) = stress_node(:, nn) ./ counter_node(nn);
end

% For using patch function conveniently.
strain_node = strain_node';
stress_node = stress_node'; 

end

