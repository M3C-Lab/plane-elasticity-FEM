clc; clear all;
% ---------- Data given ----------
tic;
disp('1) Data given');

msh = load_gmsh2('example4-2.msh');

num_Diri = 2; % The number of Dirichlet boundaries.

Diri_BC = cell(2, 1);
Diri_BC{1} = @(x, y) 0;
Diri_BC{2} = @(x, y) 0;
% Dirichlet boundary conditions, givens as scalar field functions.

num_Neum = 3; % The number of Neumann boundaries.

Neum_BC = cell(3, 1);

Neum_BC{1} = @(x, y) [-10, 0]';
Neum_BC{2} = @(x, y) [0, 0]';
Neum_BC{3} = @(x, y) [0, 0]';
% Neumann boundary conditions as vector field functions in physical space.

% Neum_BC{1} = @(x, y) 10;
% Neum_BC{2} = @(x, y) 0;
% Neum_BC{3} = @(x, y) 0;
% Neumann boundary conditions as scalar field functions on [-1,1].

E = 1.0e5;   % Young's modulus.
nu = 0.3;   % Poisson's ratio.

f_x = @(x, y) 0 * x * y; % Body force field of x-component.
f_y = @(x, y) 0 * x * y; % Body force field of y-component.

element_type = 3; %The type of elements is triangle.
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
ID_array = rearrange_ID_array(ID_array, Diri_Nodes{1}, 'x');
ID_array = rearrange_ID_array(ID_array, Diri_Nodes{2}, 'y');
% One boundary can hold 2 Diri_BC by the 2 degree of freedom of each nodes.

% Set the value of trial solution on Dirichlet BC.
trial_solution = zeros(2 * msh.nbNod, 1);
% One loop for one Diri_BC
for ii = 1 : length(Diri_Nodes{1})
    trial_solution(2 * Diri_Nodes{1}(ii) -1) = ...
        Diri_BC{1}(msh.POS(Diri_Nodes{1}(ii), 1), msh.POS(Diri_Nodes{1}(ii), 2));
    % trial_solution(2 * Diri_Nodes{1}(ii) - 1) means the x-component of the node.
end

for ii = 1 : length(Diri_Nodes{2})
    trial_solution(2 * Diri_Nodes{2}(ii)) = ...
        Diri_BC{2}(msh.POS(Diri_Nodes{2}(ii), 1), msh.POS(Diri_Nodes{2}(ii), 2));
    % trial_solution(2 * Diri_Nodes{2}(ii)) means the y-component of the node.
end 
% The data structure here remained to be improved.

LM_array = make_LM_array(Plane_IEN, ID_array);

Neum_normalvector = make_normalvector_tri(Neum_IEN, msh.POS);
toc;

% ---------- Construct stiffness K and load F ---------- 
tic;
disp('3) Construct stiffness K and load F');
% The physical equations for plane-strain problem: σ = D * ε .
% σ = [σ_xx, σ_yy, τ_xy]'
% ε = [ε_xx, ε_yy, γ_xy]'
% Note: γ_xy = 2 * ε_xy
% ε_ij = 0.5 * (ui_xj + uj_xi)
% This formulation are available on Page 83.
D = [nu*E/((1+nu)*(1-2*nu))+E/(1+nu), nu*E/(1+nu)/(1-2*nu), 0;
    nu*E/(1+nu)/(1-2*nu), nu*E/((1+nu)*(1-2*nu))+E/(1+nu), 0;
    0, 0, 0.5*E/(1+nu)];

n_eq = ID_array(2, msh.nbNod);  % The number of equations.
K = zeros(n_eq, n_eq); % The global stiffness matrix K.
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
                 0, dN1_dy, 0, dN2_dy, 0, dN2_dy;
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
    for ii = 1 : 1
        for jj = 1 : size(Neum_IEN{ii}, 2)
            % Search for triangular elements adjoining  Neumann boundaries.
            if ee == Neum_IEN{ii}(3, jj)
                % Search for nodes of line elements on Neumann boundaries.
                supported_N = [1, 2, 3];
                for kk = 1 : 3
                    if Plane_IEN(kk, ee) == Neum_IEN{ii}(4, jj)
                        supported_N(kk) = [ ];
                        % 'kk'th node is the other node of the riangular
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

nbSampling = 1;
totalSampling = nbSampling * (msh.nbElm - msh.nbLines);    
% The number of sampling points in each element.

position_Sampling = zeros(2, totalSampling);   % [x, y]'
displacement_Sampling = zeros(2, totalSampling);   % [d_x, d_y]'

strain_Sampling = zeros(3, totalSampling); % [ε_xx, ε_yy, γ_xy]'
stress_Sampling = zeros(3, totalSampling); % [σ_xx, σ_yy, τ_xy]'
max_sigma_xx = 0;
dangerous_point = [0,0];

rng('shuffle');
for ee = 1 : msh.nbElm - msh.nbLines
    n1 = Plane_IEN(1, ee);
    p1 = [msh.POS(n1, 1), msh.POS(n1, 2)]';
    n2 = Plane_IEN(2, ee);
    p2 = [msh.POS(n2, 1), msh.POS(n2, 2)]';
    n3 = Plane_IEN(3, ee);
    p3 = [msh.POS(n3, 1), msh.POS(n3, 2)]';
    
    [J, phys2rst] = Mapping_p2rst(p1, p2, p3);
    
    temp = 1;
    for ss = 1 : nbSampling
        sampling = Random_tricoor();
        
        samp_x = sampling(1) * p1(1) + sampling(2) * p2(1) + sampling(3) * p3(1);
        samp_y = sampling(1) * p1(2) + sampling(2) * p2(2) + sampling(3) * p3(2);
        position_Sampling(1, nbSampling*(ee-1) + temp) = samp_x;
        position_Sampling(2, nbSampling*(ee-1) + temp) = samp_y;
        
        node_d = [trial_solution(2*n1-1);
                  trial_solution(2*n1);
                  trial_solution(2*n2-1);
                  trial_solution(2*n2);
                  trial_solution(2*n3-1);
                  trial_solution(2*n3)];
        samp_d_x = sampling(1) * node_d(1) + sampling(2) * node_d(3) + ...
                   sampling(3) * node_d(5);
        samp_d_y = sampling(1) * node_d(2) + sampling(2) * node_d(4) + ...
                   sampling(3) * node_d(6);        
        displacement_Sampling(1, nbSampling*(ee-1) + temp) = samp_d_x;
        displacement_Sampling(2, nbSampling*(ee-1) + temp) = samp_d_y;
        
        dN1_dx = TriangularBasis(1, 1, 0, samp_x, samp_y, phys2rst);
        dN2_dx = TriangularBasis(2, 1, 0, samp_x, samp_y, phys2rst);
        dN3_dx = TriangularBasis(3, 1, 0, samp_x, samp_y, phys2rst);
        dN1_dy = TriangularBasis(1, 0, 1, samp_x, samp_y, phys2rst);
        dN2_dy = TriangularBasis(2, 0, 1, samp_x, samp_y, phys2rst);
        dN3_dy = TriangularBasis(3, 0, 1, samp_x, samp_y, phys2rst); 
        B_samp = [dN1_dx, 0, dN2_dx, 0, dN3_dx, 0;
                 0, dN1_dy, 0, dN2_dy, 0, dN2_dy;
                 dN1_dy, dN1_dx, dN2_dy, dN2_dx, dN3_dy, dN3_dx];
        
        epsilon_samp = B_samp * node_d;
        strain_Sampling(1, nbSampling*(ee-1) + temp) = epsilon_samp(1);
        strain_Sampling(2, nbSampling*(ee-1) + temp) = epsilon_samp(2);
        strain_Sampling(3, nbSampling*(ee-1) + temp) = epsilon_samp(3);
        
        sigma_samp = D * epsilon_samp;
        stress_Sampling(1, nbSampling*(ee-1) + temp) = sigma_samp(1);
        stress_Sampling(2, nbSampling*(ee-1) + temp) = sigma_samp(2);
        stress_Sampling(3, nbSampling*(ee-1) + temp) = sigma_samp(3);
        
        if max_sigma_xx < sigma_samp(1)
            max_sigma_xx = sigma_samp(1);
            dangerous_point = [samp_x, samp_y]';
        end
        
        temp = temp + 1;
    end
end

% New [x, y]' after the deformation.
newposition_sampling = position_Sampling + displacement_Sampling;

toc;





