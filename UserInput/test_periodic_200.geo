// Gmsh project created on Fri Aug 07 17:07:15 2020
SetFactory("OpenCASCADE");
//+
Point(1) = {-1, -1, 0, 1.0};
//+
Point(2) = {1, -1, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {-1, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Periodic Curve {-1} = {3};
//+
Periodic Curve {-2} = {4};
//+
Plane Surface(1) = {1};
//+
Physical Curve("Periodic") = {1,2,3,4};
//+
Physical Surface("Fluid") = {1};
//+
Transfinite Curve {1} = 201 Using Progression 1;
Transfinite Curve {3} = 201 Using Progression 1;
Transfinite Curve {2} = 201 Using Progression 1;
Transfinite Curve {4} = 201 Using Progression 1;
//+
Transfinite Surface {1};
