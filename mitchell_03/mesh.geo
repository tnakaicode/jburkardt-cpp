cl__1 = 1;

Point(1) = {1, 0, 0, 1};
Point(2) = {1, 1, 0, 1};
Point(3) = {-1, 1, 0, 1};
Point(4) = {-1, -1, 0, 1};
Point(5) = {1, -1, 0, 1};
Point(6) = {1, -1e-8, 0, 1};
Point(7) = {0, 0, 0, 1};

Line(1) = {7,1};
Line(2) = {1,2};
Line(3) = {2,3};
Line(4) = {3,4};
Line(5) = {4,5};
Line(6) = {5,6};
Line(7) = {6,7};

Line Loop(1) = {1, 2, 3, 4, 5, 6, 7};
Plane Surface(1) = {1};
Physical Surface(1) = {1};

Physical Line(0) = {1, 2, 3, 4, 5, 6, 7};
