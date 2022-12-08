// Gmsh project created on Tue Oct 13 09:04:38 2020
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0.3, 0, 0, 1.0};
//+
Point(3) = {1.0, 0.061, 0, 1.0};
//+
Point(4) = {1.0, 1, 0, 1.0};
//+
Point(5) = {0, 1, 0, 1.0};
//+
Point(6) = {0.3, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 6};
//+
Line(5) = {6, 5};
//+
Line(6) = {5, 1};
//+
Line(7) = {2, 6};
//+
Curve Loop(1) = {1, 7, 5, 6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, 3, 4, -7};
//+
Plane Surface(2) = {2};
//+
Physical Curve("Inlet") = {6};
//+
Physical Curve("Wall") = {1, 2, 4, 5};
//+
Physical Curve("Outlet") = {3};
//+
Physical Surface("Fluid") = {1, 2};
//+
Transfinite Curve {1, 5} = 30 Using Progression 1;
//+
Transfinite Curve {2, 4} = 70 Using Progression 1;
//+
Transfinite Curve {6, 7, 3} = 100 Using Progression 1;
//+
Recombine Surface {1, 2};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
