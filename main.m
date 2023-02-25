clc; clear all;

msh = load_gmsh2('example4-2.msh');

num_Diri = 2;
num_Neum = 3;

Diri_Nodes = make_Diri_Nodes(msh, num_Diri);

Plane_IEN = make_Plane_IEN(msh);

Neum_IEN = make_Neum_IEN(msh, num_Diri, num_Neum, Plane_IEN);

ID_array = make_ID_array(msh.nbNod);
ID_array = rearrange_ID_array(ID_array, Diri_Nodes{1}, 'x');
ID_array = rearrange_ID_array(ID_array, Diri_Nodes{2}, 'y');

LM_array = make_LM_array(Plane_IEN, ID_array);

Neum_normalvector = make_normalvector(Neum_IEN, msh.POS);
