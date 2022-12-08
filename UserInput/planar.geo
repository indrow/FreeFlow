// Gmsh project created on Fri Oct  2 06:22:57 2020
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {2, 0.2, 0, 1.0};
//+
Point(3) = {2, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Curve("Inlet") = {4};
//+
Physical Curve("Outlet") = {2};
//+
Physical Curve("Wall") = {1, 3};
//+
Transfinite Curve {4, 2} = 21 Using Progression 1;
//+
Transfinite Curve {3, 1} = 41 Using Progression 1;
//+
Transfinite Surface {1};
//+
Physical Surface("Fluid") = {1};
