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

% The physical equations for plane-strain problem: �� = D * �� .
% �� = [��_xx, ��_yy, ��_xy]'
% �� = [��_xx, ��_yy, ��_xy]'
% Note: ��_xy = 2 * ��_xy
% ��_ij = 0.5 * (ui_xj + uj_xi)
% This formulation are available on Page 83.
D = [nu*E/((1+nu)*(1-2*nu))+E/(1+nu), nu*E/(1+nu)/(1-2*nu), 0;
    nu*E/(1+nu)/(1-2*nu), nu*E/((1+nu)*(1-2*nu))+E/(1+nu), 0;
    0, 0, 0.5*E/(1+nu)];

n_eq = ID_array(2, msh.nbNod);  % The number of equations.
K = sparse(n_eq, n_eq); % The global stiffness matrix K.
F = zeros(n_eq, 1); % The global load vector F.

if etype == 3
%     % The nodes of the parent triangular element in the (xi, eta) space.
%     p1 = [-1, 0]';
%     p2 = [1, 0]';
%     p3 = [0, 1]';
%     % Note: The parametric mapping is given here:
%     % x = (x2 - x1) * xi + (x3 - x1) * eta + x1
%     % y = (y2 - y1) * xi + (y3 - y1) * eta + y1
%     % It makes p1 --> n1 , p2 --> n2  ,  p3 --> n3
%     % n1 = (x1, y1) , n2 = (x2, y2) , n3 = (x3, y3) in the physical space.
%     % It makes point p in triangle-p1p2p3 --> n in triangle-n1n2n3.
%     % This formulation is also applied for Mapping.m
%     
%     % Quadrature rule for parent triangle.
%     [qp, wq, nqp] = TriangularQuad(Quad_degree, p1, p2, p3);
    
    dof_e = 2 * etype;   % The degree of freedom of the nodes of an element.

    k_cell = cell(msh.nbTriangles, 1);
    for ee = 1 : msh.nbTriangles
        k_ele = zeros(dof_e, dof_e);
        f_ele = zeros(dof_e, 1);
        
        % The position of the nodes of an element in the physical space.
        n1 = [msh.POS(Plane_IEN(1, ee), 1), msh.POS(Plane_IEN(1, ee), 2)]';
        n2 = [msh.POS(Plane_IEN(2, ee), 1), msh.POS(Plane_IEN(2, ee), 2)]';
        n3 = [msh.POS(Plane_IEN(3, ee), 1), msh.POS(Plane_IEN(3, ee), 2)]';
        x_ele = [n1, n2, n3];
        % Note: In gmsh file the nodes 1,2,3 of an triangular element are
        % automatically marked counterclockwise, which is suitable for our
        % parent triangle.
        [qp, wq, nqp] = TriangularQuad(Quad_degree, n1, n2, n3);
        
%         [para2phys,phys2para] = Mapping(n1, n2, n3, p1, p2, p3);
        
        % Jacobian matrix for mapping.
%         Jacobian_matrix = [para2phys(1, 1), para2phys(1, 2);
%                            para2phys(2, 1), para2phys(2, 2)];
%         J = det(Jacobian_matrix);
        
        B = zeros(3, dof_e);
        for qua = 1 : nqp
% The term in the weak form: ��(w)' * D * ��(u)
%        [d_dx, 0;
% ��(u) =  0, d_dy;      *   [u, v] = L * u
%        d_dy, d_dx]
% u = [N1, 0, N2, 0, N3, 0;     *   [u1, v1, u2, v2, u3, v3]' = N * d
%       0, N1, 0, N2, 0, N3]
% B = L * N
% ��(w)' * D * ��(u) = B' * D * B
% ���ڡ����޵�Ԫ������59��60ҳ��
% If we use mapping from physical space to parametric space, we should
%   consider the general differential relations,
%   and Jacobian J should be multiplied with the infinitesimal volume for
%   the Integral.
%             dN1_dxi = TriangularBasis(1,1,0,qp(1,nqp),qp(2,nqp), p1, p2, p3);
%             dN2_dxi = TriangularBasis(2,1,0,qp(1,nqp),qp(2,nqp), p1, p2, p3);
%             dN3_dxi = TriangularBasis(3,1,0,qp(1,nqp),qp(2,nqp), p1, p2, p3);
%             dN1_deta = TriangularBasis(1,0,1,qp(1,nqp),qp(2,nqp), p1, p2, p3);
%             dN2_deta = TriangularBasis(2,0,1,qp(1,nqp),qp(2,nqp), p1, p2, p3);
%             dN3_deta = TriangularBasis(3,0,1,qp(1,nqp),qp(2,nqp), p1, p2, p3); 
%             dN1_dx = dN1_dxi * phys2para(1, 1) + dN1_deta * phys2para(2, 1);
%             dN2_dx = dN2_dxi * phys2para(1, 1) + dN2_deta * phys2para(2, 1);
%             dN3_dx = dN3_dxi * phys2para(1, 1) + dN3_deta * phys2para(2, 1);
%             dN1_dy = dN1_dxi * phys2para(1, 2) + dN1_deta * phys2para(2, 2);
%             dN2_dy = dN2_dxi * phys2para(1, 2) + dN2_deta * phys2para(2, 2);
%             dN3_dy = dN3_dxi * phys2para(1, 2) + dN3_deta * phys2para(2, 2);
            dN1_dx = TriangularBasis(1,1,0,qp(1,nqp),qp(2,nqp), n1, n2, n3);
            dN2_dx = TriangularBasis(2,1,0,qp(1,nqp),qp(2,nqp), n1, n2, n3);
            dN3_dx = TriangularBasis(3,1,0,qp(1,nqp),qp(2,nqp), n1, n2, n3);
            dN1_dy = TriangularBasis(1,0,1,qp(1,nqp),qp(2,nqp), n1, n2, n3);
            dN2_dy = TriangularBasis(2,0,1,qp(1,nqp),qp(2,nqp), n1, n2, n3);
            dN3_dy = TriangularBasis(3,0,1,qp(1,nqp),qp(2,nqp), n1, n2, n3);
            
            B_qua = [dN1_dx, 0, dN2_dx, 0, dN3_dx, 0;
                    0, dN1_dy, 0, dN2_dy, 0, dN2_dy;
                    dN1_dy, dN1_dx, dN2_dy, dN2_dx, dN3_dy, dN3_dx];
            %B = B + wq(qua) * B_qua;
            k_ele = k_ele + wq(qua) * B_qua' * D * B_qua; 
        end
        %k_ele = B' * D * B;
        k_cell{ee, 1} = k_ele;
        
    end
end