clc; clear all;

% ---------- Data given ----------
msh = load_gmsh2('example4-2.msh');

num_Diri = 2;
num_Neum = 3;

E = 1.0e5;   % Young's modulus.
nu = 0.3;   % Poisson's ratio.

% ---------- Preprocess ----------
Diri_Nodes = make_Diri_Nodes(msh, num_Diri);

Plane_IEN = make_Plane_IEN(msh);

Neum_IEN = make_Neum_IEN(msh, num_Diri, num_Neum, Plane_IEN);

ID_array = make_ID_array(msh.nbNod);
% Input Dirichlet BC
ID_array = rearrange_ID_array(ID_array, Diri_Nodes{1}, 'x');
ID_array = rearrange_ID_array(ID_array, Diri_Nodes{2}, 'y');

LM_array = make_LM_array(Plane_IEN, ID_array);

Neum_normalvector = make_normalvector(Neum_IEN, msh.POS);

% ---------- Construct stiffness K ---------- 
n_eq = ID_array(2, msh.nbNod);
K = zeros(n_eq, n_eq);

% The physical equations for plane-strain problem: sigma = D * epsilon .
% sigma = [sigma_xx, sigma_yy, tau_xy]'
% epsilon = [epsilon_xx, epsilon_yy, gama_xy]'
% Note: gama_xy = 2 * epsilon_xy
% epsilon_ij = 0.5 * (ui_xj + uj_xi)
inv_D = [(1-nu^2)/E,  -nu/(1-nu),   0;
         -nu/(1-nu),  (1-nu^2)/E,   0;
         0,            0,    2*(1+nu)/E];
D = inv(inv_D);

[qp, wq, nqp] = TriangularQuad(7);
