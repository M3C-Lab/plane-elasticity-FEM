clc; clear all;

% ---------- Data given ----------
msh = load_gmsh2('example4-2.msh');

num_Diri = 2;
num_Neum = 3;

E = 1.0e5;   % Young's modulus.
nu = 0.3;   % Poisson's ratio.

etype = 3; %The type of elements is triangle.
Quad_degree = 4; % The degree of precision of numerical quadrature.

% ---------- Preprocess ----------
Diri_Nodes = make_Diri_Nodes(msh, num_Diri);

Plane_IEN = make_Plane_IEN(msh, etype);

Neum_IEN = make_Neum_IEN_tri(msh, num_Diri, num_Neum, Plane_IEN);

ID_array = make_ID_array(msh.nbNod);
% Input Dirichlet BC.
ID_array = rearrange_ID_array(ID_array, Diri_Nodes{1}, 'x');
ID_array = rearrange_ID_array(ID_array, Diri_Nodes{2}, 'y');

LM_array = make_LM_array(Plane_IEN, ID_array);

Neum_normalvector = make_normalvector_tri(Neum_IEN, msh.POS);

% ---------- Construct stiffness K ---------- 

% The physical equations for plane-strain problem: sigma = D * epsilon .
% sigma = [sigma_xx, sigma_yy, tau_xy]'
% epsilon = [epsilon_xx, epsilon_yy, gama_xy]'
% Note: gama_xy = 2 * epsilon_xy
% epsilon_ij = 0.5 * (ui_xj + uj_xi)
inv_D = [(1-nu^2)/E,  -nu/(1-nu),   0;
         -nu/(1-nu),  (1-nu^2)/E,   0;
         0,            0,    2*(1+nu)/E];
D = inv(inv_D);

n_eq = ID_array(2, msh.nbNod);  % The number of equations.
K = sparse(n_eq, n_eq); % The global stiffness matrix K.
F = zeros(n_eq, 1); % The global load vector F.

if etype == 3
    % The nodes of the parent triangular element.
    p1 = [-1,-1]';
    p2 = [1, -1]';
    p3 = [0, 1]';
    % Quadrature rule for parent triangle.
    [qp, wq, nqp] = TriangularQuad(Quad_degree, p1, p2, p3);
    
    dof_e = 6;   % The degree of freedom of the nodes of an element.

    for ee = 1 : msh.nbTriangles
        k_ele = zeros(dof_e, dof_e);
        f_ele = zeros(dof_e, 1);
        
        % The position of the nodes of an element in the physical space.
        x_ele = zeros(2, etype);
        for aa = 1 : etype
            x_ele(1, aa) = msh.POS(Plane_IEN(aa, ee), 1);
            x_ele(2, aa) = msh.POS(Plane_IEN(aa, ee), 2);
        end
        % Note: In gmsh file the nodes 1,2,3 of an triangular element are
        % automatically marked counterclockwise, which is suitable for our
        % parent triangle.
        
        
    end
    
end




