//+
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1.5, -4.0, 0, 1.0};
Point(4) = {2, 0, 0, 1.0};
Point(5) = {3, 0, 0, 1.0};
Point(6) = {3, 1, 0, 1.0};
Point(7) = {2, 1, 0, 1.0};
Point(8) = {1, 1, 0, 1.0};
Point(9) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Circle(2) = {2, 3, 4};
//+
Line(3) = {4, 5};
//+
Line(4) = {5, 6};
//+
Line(5) = {6, 7};
//+
Line(6) = {7, 8};
//+
Line(7) = {8, 9};
//+
Line(8) = {9, 1};
//+
Line(9) = {2, 8};
//+
Line(10) = {4, 7};
//+
Curve Loop(1) = {1, 9, 7, 8};
//+
Surface(1) = {1};
//+
Curve Loop(2) = {2, 10, 6, -9};
//+
Surface(2) = {2};
//+
Curve Loop(3) = {3, 4, 5, -10};
//+
Surface(3) = {3};
//+
Transfinite Curve {8, 1, 9, 7, 6, 2, 10, 3, 4, 5} = 101 Using Progression 1;
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Recombine Surface {1, 2, 3};
//+
Transfinite Surface {3};
//+
Recombine Surface {1, 2, 3};
//+
Physical Curve("Inlet") = {8};
//+
Physical Curve("Outlet") = {4};
//+
Physical Curve("Wall") = {1, 2, 3, 7, 6, 5};
//+
Physical Surface("Fluid") = {1, 2, 3};
