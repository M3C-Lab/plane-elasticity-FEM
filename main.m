clc; clear all;

% ---------- Data given ----------
msh = load_gmsh2('example4-2.msh');

num_Diri = 2; % The number of Dirichlet boundaries.

Diri_BC = cell(2, 1);
Diri_BC{1} = @(x, y) 0 * x * y;
Diri_BC{2} = @(x, y) 0 * x * y;
% Dirichlet boundary conditions, givens as field functions.

num_Neum = 3; % The number of Neumann boundaries.

Neum_BC = cell(3, 1);
Neum_BC{1} = @(x, y) [-10, 0];
Neum_BC{2} = @(x, y) [0, 0];
Neum_BC{3} = @(x, y) [0, 0];
% Neumann boundary conditions, givens as field functions.

E = 1.0e5;   % Young's modulus.
nu = 0.3;   % Poisson's ratio.

f_x = @(x, y) 0 * x * y; % Body force field of x-component.
f_y = @(x, y) 0 * x * y; % Body force field of y-component.

element_type = 3; %The type of elements is triangle.
Quad_degree = 4; % The degree of precision of numerical quadrature.

% ---------- Preprocess ----------

Plane_IEN = make_Plane_IEN(msh, element_type);

Neum_IEN = make_Neum_IEN_tri(msh, num_Diri, num_Neum, Plane_IEN);

ID_array = make_ID_array(msh.nbNod);
% Input Dirichlet BC.
Diri_Nodes = make_Diri_Nodes(msh, num_Diri);
ID_array = rearrange_ID_array(ID_array, Diri_Nodes{1}, 'x');
ID_array = rearrange_ID_array(ID_array, Diri_Nodes{2}, 'y');
% Set the value of trial solution on Dirichlet BC.
trial_solution = zeros(2 * msh.nbNod, 1);
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

% ---------- Construct stiffness K and load F ---------- 

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

if element_type == 3
    
    % Quadrature rule for parent triangle.
    [qp, wq, nqp] = TriangularQuad(Quad_degree);
    
    % Quadrature rule for parent line element (on the Neumann boundaries).
    [qpl, wql] = Gauss(Quad_degree, -1, 1);
    
    dof_e = 2 * element_type;   % The degree of freedom of the nodes of an element.
    k_cell = cell(msh.nbTriangles, 1);
    f_cell = cell(msh.nbTriangles, 1);
    
    for ee = 1 : msh.nbTriangles
        k_ele = zeros(dof_e, dof_e);
        f_ele = zeros(dof_e, 1);
        
        % The position of the nodes of an element in the physical space.
        p1 = [msh.POS(Plane_IEN(1, ee), 1), msh.POS(Plane_IEN(1, ee), 2)]';
        p2 = [msh.POS(Plane_IEN(2, ee), 1), msh.POS(Plane_IEN(2, ee), 2)]';
        p3 = [msh.POS(Plane_IEN(3, ee), 1), msh.POS(Plane_IEN(3, ee), 2)]';
        % Note: In gmsh file the nodes 1,2,3 of an triangular element are
        % automatically marked counterclockwise, which is suitable for our
        % parent triangle.
        
        [J, phys2rst] = Mapping_p2rst(p1, p2, p3);      
        % The det(Jacobian matrix) for quadrature,
        % and the mapping matrix for basis functions.
        
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
        
        k_cell{ee} = k_ele;
        f_cell{ee} = f_ele;
        
        
        for aa = 1 : dof_e
            LM_a = LM_array(aa, ee);
            if LM_a > 0
                F(LM_a) = F(LM_a) + f_ele(aa);
                for bb = 1 : dof_e
                    LM_b = LM_array(bb, ee);
                    if LM_b > 0
                        K(LM_a, LM_b) = K(LM_a, LM_b)+ k_ele(aa, bb);
                    else
                        % For Dirichlet BC.
                        F(LM_a) = F(LM_a) - k_ele(aa, bb) * trial_solution(bb);
                    end
                end
            end
        end
        
%         for ii = num_Neum
%             for jj = 1 : size(Neum_IEN{ii}, 1)
%                 % Search for triangular elements adjoining  Neumann boundaries.
%                 if ee == Neum_IEN{ii}(3, jj)
%                     % Search for nodes of line elements on Neumann boundaries.
%                     
%                 end
%             end
%         end
        
        
      
    end
    
end




