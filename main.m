clc; clear all;

msh = load_gmsh2('example4-2.msh');

num_Diri = 2;
num_Neum = 3;

Diri_Nodes = make_Diri_Nodes(msh, num_Diri);

Plane_IEN = make_Plane_IEN(msh);

Neum_IEN = make_Neum_IEN(msh, num_Diri, num_Neum, Plane_IEN);
