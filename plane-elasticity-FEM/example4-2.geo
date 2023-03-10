// Gmsh project created on Thu Feb 16 14:44:35 2023
SetFactory("OpenCASCADE");
//+
Mesh.MshFileVersion = 2;
//+
Point(1) = {0, 0, 0, 0.1};
//+
Point(2) = {-1, 0, 0, 0.1};
//+
Point(3) = {0, 1, 0, 0.1};
//+
Point(4) = {0, 4, 0, 0.1};
//+
Point(5) = {-4, 4, 0, 0.1};
//+
Point(6) = {-4, 0, 0, 0.1};
//+
Line(1) = {3, 4};
//+
Line(2) = {4, 5};
//+
Line(3) = {5, 6};
//+
Line(4) = {6, 2};
//+
Circle(5) = {2, 1, 3};
//+
Curve Loop(1) = {3, 4, 5, 1, 2};
//+
Plane Surface(1) = {1};
//+
Physical Curve("Diri-1") = {1};
//+
Physical Curve("Diri-2") = {4};
//+
Physical Curve("Neum-1") = {3};
//+
Physical Curve("Neum-2") = {2};
//+
Physical Curve("Neum-3") = {5};
//+
Physical Surface("plane-1") = {1};