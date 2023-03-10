function [X, Y, XX, YY, strain_sampling, stress_sampling] = ...
    sampling_data(nbSampling, msh, Plane_IEN, element_type, trial_solution)
% To compute strain and stress on some sampling points in each element.
% Input:
%   nbSampling: The number of sampling points in each element.
%   msh: the Gmsh info imported.
%   Plane_IEN: The IEN array of plane elment.
%   element_type: 3 --- triangle
%   trial_solution: The solution data.
% Output:
%   X: Orignal X-coordinate of each sampling point.
%   Y: Orignal Y-coordinate of each sampling point.
%   XX: X after the deformation.
%   YY: Y after the deformation.
%   strain_node: [¦Å_xx, ¦Å_yy, ¦Ã_xy]
%       strain_node(i, *): The 'i'th sampling point.
%   stress_node: [¦Ò_xx, ¦Ò_yy, ¦Ó_xy]
%       stress_node(i, *): The 'i'th sampling point.

global D;

rng('shuffle');
totalSampling = nbSampling * (msh.nbElm - msh.nbLines);

X = zeros(totalSampling, 1);
Y = zeros(totalSampling, 1); % [x, y]

d_x_sampling = zeros(totalSampling, 1);
d_y_sampling = zeros(totalSampling, 1);   % [d_x, d_y]

strain_sampling = zeros(1, totalSampling); % [¦Å_xx, ¦Å_yy, ¦Ã_xy]'

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
    
    temp = 1;
    if element_type == 3
        for ss = 1 : nbSampling
            sampling = Random_tricoor();
        
            samp_x = sampling(1) * p1(1) + sampling(2) * p2(1) + sampling(3) * p3(1);
            samp_y = sampling(1) * p1(2) + sampling(2) * p2(2) + sampling(3) * p3(2);
            X(nbSampling*(ee-1) + temp) = samp_x;
            Y(nbSampling*(ee-1) + temp) = samp_y;
        
            samp_d_x = sampling(1) * d_ele(1) + sampling(2) * d_ele(3) + ...
                    sampling(3) * d_ele(5);
            samp_d_y = sampling(1) * d_ele(2) + sampling(2) * d_ele(4) + ...
                    sampling(3) * d_ele(6);        
            d_x_sampling(nbSampling*(ee-1) + temp) = samp_d_x;
            d_y_sampling(nbSampling*(ee-1) + temp) = samp_d_y;
        
            dN1_dx = TriangularBasis(1, 1, 0, samp_x, samp_y, phys2rst);
            dN2_dx = TriangularBasis(2, 1, 0, samp_x, samp_y, phys2rst);
            dN3_dx = TriangularBasis(3, 1, 0, samp_x, samp_y, phys2rst);
            dN1_dy = TriangularBasis(1, 0, 1, samp_x, samp_y, phys2rst);
            dN2_dy = TriangularBasis(2, 0, 1, samp_x, samp_y, phys2rst);
            dN3_dy = TriangularBasis(3, 0, 1, samp_x, samp_y, phys2rst); 
            B_samp = [dN1_dx, 0, dN2_dx, 0, dN3_dx, 0;
                    0, dN1_dy, 0, dN2_dy, 0, dN2_dy;
                    dN1_dy, dN1_dx, dN2_dy, dN2_dx, dN3_dy, dN3_dx];
        
            epsilon_samp = B_samp * d_ele;
            strain_sampling(1, nbSampling*(ee-1) + temp) = epsilon_samp(1);
            strain_sampling(2, nbSampling*(ee-1) + temp) = epsilon_samp(2);
            strain_sampling(3, nbSampling*(ee-1) + temp) = epsilon_samp(3);
        
            temp = temp + 1;
        end
    end
end

stress_sampling = D * strain_sampling;

strain_sampling = strain_sampling';
stress_sampling = stress_sampling';

XX = X + d_x_sampling;
YY = Y + d_y_sampling;

end

