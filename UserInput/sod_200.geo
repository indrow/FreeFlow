// Gmsh project created on Tue Oct 13 23:09:51 2020
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 0.3, 0, 1.0};
//+
Point(4) = {0, 0.3, 0, 1.0};
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
Plane Surface(1) = {1};
//+
Physical Curve("ZeroGrad") = {4, 2};
//+
Physical Curve("Wall") = {3, 1};
//+
Physical Surface("Fluid") = {1};
//+
Transfinite Curve {1, 3} = 201 Using Progression 1;
//+
Transfinite Curve {4, 2} = 61 Using Progression 1;
//+
Transfinite Surface {1};
