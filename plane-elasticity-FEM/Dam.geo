// Gmsh project created on Fri Mar 10 22:37:09 2023
SetFactory("OpenCASCADE");
//+
Mesh.MshFileVersion = 2;
//+
Point(1) = {-30, 0, 0, 1.0};
//+
Point(2) = {200, 0, 0, 1.0};
//+
Point(3) = {50, 150, 0, 1.0};
//+
Point(4) = {0, 150, 0, 1.0};
//+//+
Point(5) = {75, 125, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 5};
//+
Line(3) = {5, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {4, 1};
//+
//+
Curve Loop(1) = {2, 3, 4, 5, 1};
//+
Curve Loop(2) = {5, 1, 2, 3, 4};
//+
Surface(1) = {2};
//+
Physical Curve("Diri-1") = {1};
//+
Physical Curve("Neum-1") = {2};
//+
Physical Curve("Neum-2") = {3};
//+
Physical Curve("Neum-3") = {4};
//+
Physical Curve("Neum-4") = {5};
//+
Physical Surface("Surf-1") = {1};
