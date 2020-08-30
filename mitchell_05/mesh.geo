cl__1 = 1;
Point(1) = {0, 0, 0, 1};
Point(2) = {8.4, 0, 0, 1};
Point(3) = {8.4, 24, 0, 1};
Point(4) = {0, 24, 0, 1};
Point(5) = {0, 23.2, 0, 1};
Point(6) = {8, 23.2, 0, 1};
Point(7) = {8, 0.8, 0, 1};
Point(8) = {0, 0.8, 0, 1};
Point(9) = {0, 21.2, 0, 1};
Point(10) = {6.1, 21.2, 0, 1};
Point(11) = {6.1, 18.8, 0, 1};
Point(12) = {0, 18.8, 0, 1};
Point(13) = {0, 3.6, 0, 1};
Point(14) = {6.1, 3.6, 0, 1};
Point(15) = {6.1, 1.6, 0, 1};
Point(16) = {0, 1.6, 0, 1};
Point(17) = {6.1, 0.8, 0, 1};
Point(18) = {6.5, 0.8, 0, 1};
Point(19) = {6.5, 21.2, 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};
Line(9) = {5, 9};
Line(10) = {9, 10};
Line(11) = {10, 11};
Line(12) = {11, 12};
Line(13) = {12, 9};
Line(14) = {14, 11};
Line(15) = {13, 14};
Line(16) = {13, 12};
Line(17) = {15, 14};
Line(18) = {15, 16};
Line(19) = {16, 13};
Line(20) = {15, 17};
Line(21) = {17, 8};
Line(22) = {8, 16};
Line(23) = {10, 19};
Line(24) = {17, 18};
Line(25) = {7, 18};
Line(26) = {18, 19};

//Region 1
Line Loop(1) = {1, 2, 3, 4, 5, 6, 25, -24, 21, 8};
Plane Surface(1) = {1};
Recombine Surface {1};

//Region 2 (top)
Line Loop(2) = {-12, -11, -10, -13};
Plane Surface(2) = {2};
Recombine Surface {2};

//Region 2 (bottom)
Line Loop(3) = {-18, 17, -15, -19};
Plane Surface(3) = {3};

//Region 3
Line Loop(4) = {15, -16, 12, 14};
Plane Surface(4) = {4};

//Region 4
Line Loop(5) = {11, -14, -17, 20, 24, 26, -23};
Plane Surface(5) = {5};

//Region 5 (bottom left)
Line Loop(6) = {-21, -22, 18, -20};
Plane Surface(6) = {6};

//Region 5 (top)
Line Loop(7) = {9, -5, -6, -25, -26, 23, 10};
Plane Surface(7) = {7};

//Boundary surfaces
Physical Line(1) = {1};//bottom
Physical Line(2) = {2};//left
Physical Line(3) = {3};//top
Physical Line(4) = {4, 8, 9, 13, 16, 19, 22};


cl__1 = 1;
Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Surface(6) = {6};
Physical Surface(7) = {7};
