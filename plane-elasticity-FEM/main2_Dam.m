clc; clear all; close all;
addpath('Preprocess', 'Mapping & Quadrature', 'Postprocess')
global D;
% ---------- Data given ----------
tic;
disp('1) Data given');

msh = load_gmsh2('Dam2.msh');

E = 3.0e10;   % Young's modulus.
nu = 0.2;   % Poisson's ratio.

plane_strain_or_plane_stress = 0; 
% 0:plane strain problem / 1:plane stress problem

num_Diri = 1; % The number of Dirichlet boundaries.
Diri_BC = cell(1, 1);

Diri_BC{1} = @(x, y) 0;
Diri_BC{2} = @(x, y) 0;
% Dirichlet boundary conditions, givens as scalar field functions.

num_Neum = 3; % The number of Neumann boundaries.
Neum_BC = cell(4, 1);

Neum_BC{1} = @(x, y) [-0.5 * sqrt(2) * 9800 * (125 - y),... 
    -0.5 * sqrt(2) * 9800 * (125 - y)]';
Neum_BC{2} = @(x, y) [0, 0]';
Neum_BC{3} = @(x, y) [0, 0]';
Neum_BC{4} = @(x, y) [0, 0]';
% Neumann boundary conditions as vector field functions in physical space.

f_x = @(x, y) 0 * x * y; % Body force field of x-component.
f_y = @(x, y) - 2400 * 9.8; % Body force field of y-component.

element_type = 3; % Triangular element.
Quad_degree = 4; % The degree of precision of numerical quadrature.
toc;

% ---------- Preprocess ----------
tic;
disp('2) Preprocess');

Plane_IEN = make_Plane_IEN(msh, element_type);
Neum_IEN = make_Neum_IEN_tri(msh, num_Diri, num_Neum, Plane_IEN);
ID_array = make_ID_array(msh.nbNod);

% Input Dirichlet BC.
Diri_Nodes = make_Diri_Nodes(msh, num_Diri);
% One line for one Diri_BC {
ID_array = rearrange_ID_array(ID_array, Diri_Nodes{1}, 'x');
ID_array = rearrange_ID_array(ID_array, Diri_Nodes{1}, 'y');
% } One boundary can hold 2 Diri_BC by the 2 degree of freedom of each nodes.

% Set the value of trial solution on Dirichlet BC.
trial_solution = zeros(2 * msh.nbNod, 1);
% One loop for one Diri_BC {
for ii = 1 : length(Diri_Nodes{1})
    trial_solution(2 * Diri_Nodes{1}(ii) -1) = ...
        Diri_BC{1}(msh.POS(Diri_Nodes{1}(ii), 1), msh.POS(Diri_Nodes{1}(ii), 2));
    % trial_solution(2 * Diri_Nodes{1}(ii) - 1) means the x-component of the node.
end

for ii = 1 : length(Diri_Nodes{1})
    trial_solution(2 * Diri_Nodes{1}(ii)) = ...
        Diri_BC{2}(msh.POS(Diri_Nodes{1}(ii), 1), msh.POS(Diri_Nodes{1}(ii), 2));
    % trial_solution(2 * Diri_Nodes{1}(ii)) means the y-component of the node.
end 
% } The data structure here remained to be improved.

LM_array = make_LM_array(Plane_IEN, ID_array);

Neum_normalvector = make_normalvector_tri(Neum_IEN, msh.POS);
toc;

% ---------- Construct stiffness K and load F ---------- 
tic;
disp('3) Construct stiffness K and load F');
% The physical equations for plane strain/stress problem: σ = D * ε .
% σ = [σ_xx, σ_yy, τ_xy]'
% ε = [ε_xx, ε_yy, γ_xy]'
% Note: γ_xy = 2 * ε_xy
% ε_ij = 0.5 * (ui_xj + uj_xi)
% This formulation are available on Page 83.
lamda = nu * E / ((1 + nu) * (1 -2 * nu));
mu = 0.5 * E / (1 + nu);
if plane_strain_or_plane_stress == 1
    lamda = 2 * lamda * mu / (lamda + 2 * mu);
    disp('This is a plane stress problem.');
else
    disp('This is a plane strain problem.');
end
D = [lamda + 2 * mu, lamda, 0;
    lamda, lamda + 2 * mu, 0;
    0, 0, mu];

n_eq = ID_array(2, msh.nbNod);  % The number of equations.
K = sparse(n_eq, n_eq); % The global stiffness matrix K.
F = zeros(n_eq, 1); % The global load vector F.

    
% Quadrature rule for parent triangle.
[qp, wq, nqp] = TriangularQuad(Quad_degree);
    
% Quadrature rule for parent line element (on the Neumann boundaries).
[lqp, lwq] = Gauss(Quad_degree + 1, -1, 1);
    
dof_e = 2 * element_type;   % The degree of freedom of the nodes of an element.

% Monitors.
% % % k_cell = cell(msh.nbTriangles, 1);
% % % f_cell = cell(msh.nbTriangles, 1);
% % % Nm_cell = cell(msh.nbNod, 1);
% % % for ll = 1 : msh.nbNod
% % %     Nm_cell{ll} = zeros(2, 1);
% % % end

for ee = 1 : msh.nbTriangles
    k_ele = zeros(dof_e, dof_e);
    f_ele = zeros(dof_e, 1);
        
    % The position of the nodes of an element in the physical space.
    p1 = [msh.POS(Plane_IEN(1, ee), 1), msh.POS(Plane_IEN(1, ee), 2)]';
    p2 = [msh.POS(Plane_IEN(2, ee), 1), msh.POS(Plane_IEN(2, ee), 2)]';
    p3 = [msh.POS(Plane_IEN(3, ee), 1), msh.POS(Plane_IEN(3, ee), 2)]';
    % Note: In gmsh file the nodes 1,2,3 of an triangular element are
    % automatically marked counterclockwise. 
        
    % The det(Jacobian matrix) for quadrature,
    %  and the mapping matrix for basis functions.
    [J, phys2rst] = Mapping_p2rst(p1, p2, p3);      
        
    % Construct k and f.
    for qua = 1 : nqp
% The term in the weak form: ε(w)' * D * ε(u)
%        [d_dx, 0;
% ε(u) =  0, d_dy;      *   [u, v] = L * u
%        d_dy, d_dx]
% u = [N1, 0, N2, 0, N3, 0;     *   [u1, v1, u2, v2, u3, v3]' = N * d
%       0, N1, 0, N2, 0, N3]
% B = L * N
% ε(w)' * D * ε(u) = B' * D * B
% 见于《有限单元法》第59、60页。
% If we use mapping from physical space to parametric space, we should
%   consider the general differential relations,
%   and Jacobian J should be multiplied with the infinitesimal volume for
%   the Integral.

        qua_x = qp(1,qua) * p1(1) + qp(2,qua) * p2(1) + qp(3,qua) * p3(1);
        qua_y = qp(1,qua) * p1(2) + qp(2,qua) * p2(2) + qp(3,qua) * p3(2);
            
        dN1_dx = TriangularBasis(1, 1, 0, qua_x, qua_y, phys2rst);
        dN2_dx = TriangularBasis(2, 1, 0, qua_x, qua_y, phys2rst);
        dN3_dx = TriangularBasis(3, 1, 0, qua_x, qua_y, phys2rst);
        dN1_dy = TriangularBasis(1, 0, 1, qua_x, qua_y, phys2rst);
        dN2_dy = TriangularBasis(2, 0, 1, qua_x, qua_y, phys2rst);
        dN3_dy = TriangularBasis(3, 0, 1, qua_x, qua_y, phys2rst); 

        B_qua = [dN1_dx, 0, dN2_dx, 0, dN3_dx, 0;
                 0, dN1_dy, 0, dN2_dy, 0, dN3_dy;
                 dN1_dy, dN1_dx, dN2_dy, dN2_dx, dN3_dy, dN3_dx];
            
        k_ele = k_ele + J * wq(qua) * B_qua' * D * B_qua; 
            
        for aa = 1 : dof_e
            node = ceil(aa/2);
            if mod(aa, 2) ~= 0
                f_ele(aa) = f_ele(aa) + J * wq(qua) * ...
                TriangularBasis(node, 0, 0, qua_x, qua_y, phys2rst) * f_x(qua_x, qua_y);
                else
                f_ele(aa) = f_ele(aa) + J * wq(qua) * ...
                TriangularBasis(node, 0, 0, qua_x, qua_y, phys2rst) * f_y(qua_x, qua_y);
            end
        end
    end
        
% % %     k_cell{ee} = k_ele;
% % %     f_cell{ee} = f_ele;
        
    % f = f + Neumman BC        
    for ii = 1 : num_Neum
        for jj = 1 : size(Neum_IEN{ii}, 2)
            % Search for triangular elements adjoining  Neumann boundaries.
            if ee == Neum_IEN{ii}(3, jj)
                % Search for nodes of line elements on Neumann boundaries.
                supported_N = [1, 2, 3];
                for kk = 1 : 3
                    if Plane_IEN(kk, ee) == Neum_IEN{ii}(4, jj)
                        supported_N(kk) = [ ];
                        % 'kk'th node is the other node of the triangular
                        % element, which is not lying on the Neumann boundary.
                    end
                end
                % 'Supported' means the value of basis function is not
                %   zero on (Neumann_node1, Neumann_node2).
                    
                % Neumann_node1 & Neumann_node2.
                n1 = [msh.POS(Neum_IEN{ii}(1, jj), 1), msh.POS(Neum_IEN{ii}(1, jj), 2)]';
                n2 = [msh.POS(Neum_IEN{ii}(2, jj), 1), msh.POS(Neum_IEN{ii}(2, jj), 2)]';
                    
                % The Jacobian for line quadrature,
                %  and the quadrature points in the physical space.
                [j, plqp] = Mapping_lineqp(lqp, n1, n2);
                    
                for node = 1 : 2
                    Nah_node = [0, 0]';
                    for qua = 1 : size(lqp)
                        Na_qua = TriangularBasis(supported_N(node), 0, 0,...
                            plqp(1, qua), plqp(2, qua), phys2rst);
                            
                        % If vector BC functions are given. {
                        Neum_BC_qua = Neum_BC{ii}(plqp(1, qua), plqp(2, qua));  % }
                            
                        % If scalar BC functions are given. {
%                         normalv = [Neum_normalvector{ii}(1, jj);
%                                    Neum_normalvector{ii}(2, jj)];
%                         Neum_BC_qua = Neum_BC{ii}(lqp(qua)) * normalv; % }
                            
                        Nah_node  = Nah_node + j * lwq(qua) * Na_qua * Neum_BC_qua;
                    end
                     
                    f_ele(2 * supported_N(node) - 1) = ...
                        f_ele(2 * supported_N(node) - 1) + Nah_node(1);
                    f_ele(2 * supported_N(node)) = ...
                        f_ele(2 * supported_N(node)) + Nah_node(2);
                        
% % %                     Nm_cell{Plane_IEN(supported_N(node), ee)} = ...
% % %                         Nm_cell{Plane_IEN(supported_N(node), ee)} + Nah_node;
                end
% % %                     f_cell{ee} = f_ele;
            end
        end
    end 
        
    % Assembly K and F
    for aa = 1 : dof_e
        LM_a = LM_array(aa, ee);
        if LM_a > 0
            F(LM_a) = F(LM_a) + f_ele(aa);
            for bb = 1 : dof_e
                LM_b = LM_array(bb, ee);
                if LM_b > 0
                    K(LM_a, LM_b) = K(LM_a, LM_b)+ k_ele(aa, bb);
                else
                    % F = F + Dirichlet BC.
                    F(LM_a) = F(LM_a) - k_ele(aa, bb) * trial_solution(bb);
                end
            end
        end
    end 
        
end
toc;

% ---------- Solve the equations ----------
tic;
disp('4) Solve the equations');

d = K \ F;
for nn = 1 : msh.nbNod
    for mm = 1 : 2
        if ID_array(mm, nn) ~= 0
            trial_solution(2*(nn-1) + mm, 1) = d(ID_array(mm, nn), 1);
        end
    end
end
toc;

% ---------- Postprocess ----------
tic;
disp('5) Postprocess');

% Data of nodes
[X1, Y1, XX1, YY1, strain_1, stress_1] = ...
    node_data(msh, Plane_IEN, element_type, trial_solution);

% Data of Interior sampling.
% If sample: {
nbSampling = 10; % The number of sampling points in each element.
[X2, Y2, XX2, YY2, strain_2, stress_2] = ...
    sampling_data(nbSampling, msh, Plane_IEN, element_type, trial_solution);
% }

X = X1; % 1:node / 2:sampling
Y = Y1;

sigma_xx = stress_1(:, 1); % 1:node / 2:sampling
sigma_yy = stress_1(:, 2);
tao_xy   = stress_1(:, 3);

max_sigma_xx = 0;
max_sigma_yy = 0;
max_tao_xy = 0;

% Find the most dangerous points.
for ii = 1 : length(X)
    if abs(sigma_xx(ii)) > max_sigma_xx
        max_sigma_xx = sigma_xx(ii);
        dangerous_point_xx = [X(ii), Y(ii)];
    end
    
    if abs(sigma_yy(ii)) > max_sigma_yy
        max_sigma_yy = sigma_yy(ii);
        dangerous_point_yy = [X(ii), Y(ii)];
    end
    
    if abs(tao_xy(ii)) > max_tao_xy
        max_tao_xy = tao_xy(ii);
        dangerous_point_xy = [X(ii), Y(ii)];
    end
end
disp('The Maximum stress:');
fprintf('max_σxx = %f    at [%f, %f]\n', ...
    max_sigma_xx, dangerous_point_xx(1), dangerous_point_xx(2));
fprintf('max_σyy = %f    at [%f, %f]\n', ...
    max_sigma_yy, dangerous_point_yy(1), dangerous_point_yy(2));
fprintf('max_τxy = %f    at [%f, %f]\n', ...
    max_tao_xy, dangerous_point_xy(1), dangerous_point_xy(2));
toc;

% ---------- Plot ----------
tic;
disp('6) Plot');

SHP = alphaShape(X, Y, 20, 'HoleThreshold', 0.000001);
% The third parameter has to be adjust by your element size!!!!!!!!!!!!!
TRI = alphaTriangulation(SHP);

% Deformed = alphaShape(XX1, YY1, 0.8, 'HoleThreshold', 0.000001);
% TRI = alphaTriangulation(Deformed);

MESH = alphaShape(msh.POS(:, 1), msh.POS(:, 2), 20, 'HoleThreshold', 0.000001);

figure(1)
patch('Faces', TRI, 'Vertices', [X, Y], 'facevertexCdata', sigma_xx, ...
    'edgecolor', 'none', 'facecolor', 'interp'); % 'flat'/'interp'
title('σ_x_x', 'fontsize', 16)
xlabel('X - axis (m)', 'fontsize', 13);
ylabel('Y - axis (m)', 'fontsize', 13);
colormap('jet');
col1 = colorbar;
set(get(col1, 'title'), 'string', '(pa)', 'fontsize', 12);
set(col1, 'Fontsize', 12);
set(gcf, 'unit', 'centimeters', 'position', [3 20 20 17.5]);

figure(2)
patch('Faces', TRI, 'Vertices', [X, Y], 'facevertexCdata', sigma_yy, ...
    'edgecolor', 'none', 'facecolor', 'interp');
title('σ_y_y', 'fontsize', 16);
xlabel('X - axis (m)', 'fontsize', 13);
ylabel('Y - axis (m)', 'fontsize', 13);
colormap('jet');
col2 = colorbar;
set(get(col2, 'title'), 'string', '(pa)', 'fontsize', 12);
set(col2, 'Fontsize', 12);
set(gcf, 'unit', 'centimeters', 'position', [10 20 20 17.5]);

figure(3)
patch('Faces', TRI, 'Vertices', [X, Y], 'facevertexCdata', tao_xy, ...
    'edgecolor', 'none', 'facecolor', 'interp');
title('τ_x_y', 'fontsize', 16);
xlabel('X - axis (m)', 'fontsize', 13);
ylabel('Y - axis (m)', 'fontsize', 13);
colormap('jet');
col3 = colorbar;
set(get(col3, 'title'), 'string', '(pa)', 'fontsize', 12);
set(col3, 'Fontsize', 12);
set(gcf, 'unit', 'centimeters', 'position', [17 20 20 17.5]);

figure(4)
plot(MESH);
title('MESH', 'fontsize', 16)
xlabel('X - axis (m)', 'fontsize', 13);
ylabel('Y - axis (m)', 'fontsize', 13);
axis([-35 205 -5 155]); % Limit the plot domain. You should adjust it!!!!!!!!!!
% axis([mix_x max_x min_y max_y])
set(gcf, 'unit', 'centimeters', 'position', [24 20 20 17.5]);

toc;
disp('Done!');





