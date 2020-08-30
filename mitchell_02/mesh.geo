cl__1 = 1;

DefineConstant [omega = {5*Pi/4}];

/******************************************************************************
*******************************************************************************/
If (omega >= 7*Pi/4)

If (omega == 7*Pi/4)
	omega = 1.000001*omega;
EndIf

If (omega == 2*Pi)
	omega = 0.999999*omega;
EndIf

Point(1) = {1, 0, 0, 1};
Point(2) = {1, 1, 0, 1};
Point(3) = {-1, 1, 0, 1};
Point(4) = {-1, -1, 0, 1};
Point(5) = {1, -1, 0, 1};
Point(6) = {1, Tan(omega), 0, 1};
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

EndIf

/******************************************************************************
*******************************************************************************/
If (omega >= 5*Pi/4 && omega < 7*Pi/4)

If (omega == 5*Pi/4)
	omega = 1.000001*omega;
EndIf

Point(1) = {1, 0, 0, 1};
Point(2) = {1, 1, 0, 1};
Point(3) = {-1, 1, 0, 1};
Point(4) = {-1, -1, 0, 1};
Point(5) = {-1/Tan(omega), -1, 0, 1};
Point(6) = {0, 0, 0, 1};

Line(1) = {6,1};
Line(2) = {1,2};
Line(3) = {2,3};
Line(4) = {3,4};
Line(5) = {4,5};
Line(6) = {5,6};

Line Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};
Physical Surface(1) = {1};

Physical Line(0) = {1, 2, 3, 4, 5, 6};

EndIf

/******************************************************************************
*******************************************************************************/
If (omega >= 3*Pi/4 && omega < 5*Pi/4)

If (omega == 3*Pi/4)
	omega = 1.00001*omega;
EndIf

Point(1) = {1, 0, 0, 1};
Point(2) = {1, 1, 0, 1};
Point(3) = {-1, 1, 0, 1};
Point(4) = {-1, -Tan(omega), 0, 1};
Point(5) = {0, 0, 0, 1};

Line(1) = {5,1};
Line(2) = {1,2};
Line(3) = {2,3};
Line(4) = {3,4};
Line(5) = {4,5};

Line Loop(1) = {1, 2, 3, 4, 5};
Plane Surface(1) = {1};
Physical Surface(1) = {1};

Physical Line(0) = {1, 2, 3, 4, 5};

EndIf

/******************************************************************************
*******************************************************************************/
If (omega >= Pi/4 && omega < 3*Pi/4)

If (omega == Pi/4)
	omega = 1.00001*omega;
EndIf

Point(1) = {1, 0, 0, 1};
Point(2) = {1, 1, 0, 1};
Point(3) = {1/Tan(omega), 1, 0, 1};
Point(4) = {0, 0, 0, 1};

Line(1) = {4,1};
Line(2) = {1,2};
Line(3) = {2,3};
Line(4) = {3,4};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Physical Surface(1) = {1};

Physical Line(0) = {1, 2, 3, 4};

EndIf

/******************************************************************************
*******************************************************************************/
If (omega < Pi/4)

Point(1) = {1, 0, 0, 1};
Point(2) = {1, Tan(omega), 0, 1};
Point(3) = {0, 0, 0, 1};

Line(1) = {3,1};
Line(2) = {1,2};
Line(3) = {2,3};

Line Loop(1) = {1, 2, 3};
Plane Surface(1) = {1};
Physical Surface(1) = {1};

Physical Line(0) = {1, 2, 3};

EndIf


