D = 0.1;
size = 1e-0;

Xi = 0.9; // um
L = 0.2; // um
 
x0 = Xi + L/2.0;
R = -0.001;   // um
f0 = -10.0;   // 0--1 

Point (1) = {10*D, 0.05*D + D/2, 0, size/25};
Point (2) = {10*D + D/2, 0.05*D + D/2, 0, size/25};
Point (3) = {10*D, 0.05*D + D, 0, size/25};
Point (4) = {10*D - D/2, 0.05*D + D/2, 0, size/25};
Point (5) = {10*D, 0.05*D, 0, size/25};
Point (6) = {0, 0, 0, size/25};
Point (7) = {10*D - D, 0, 0, size/25};
Point (8) = {10*D + D, 0, 0, size/25};
Point (9) = {30*D, 0, 0, size/25};
Point (10) = {30*D, 4*D, 0, size/25};
Point (11) = {0, 4*D, 0, size/25};

Circle (1) = {2, 1, 3};
Circle (2) = {3, 1, 4};
Circle (3) = {4, 1, 5};
Circle (4) = {5, 1, 2};
Line (5) = {6, 11};
Line (6) = {11, 10};
Line (7) = {6, 7};

pList[0] = 7; // First point label
nPoints = 10; // Number of discretization points
For i In {1 : nPoints}
  //x = (10*D - D) + L*i/(nPoints + 1);
  x = Xi + L*i/(nPoints + 1);
  pList[i] = newp;
  //Point(pList[i]) = {x, 0.25*Cos(Pi*x), 0, size/50};
  Point(pList[i]) = {x, ( R * (1 - f0/2 *(1 + Cos(2.0*Pi * (x-x0)/L) ) )), 0, size/25};
EndFor
pList[nPoints+1] = 8; // Last point label
 
Spline(8) = pList[];

//Line (8) = {7,8};
Line (9) = {8, 9};
Line (10) = {9, 10};

Physical Line (1) = {5}; //left
Physical Line (2) = {6}; //top
//Physical Line (3) = {7}; //bottom_left
Physical Line (3) = {7, 8, 9}; //bottom
//Physical Line (5) = {9}; //bottom_right
Physical Line (4) = {10}; //right
Physical Line (5) = {1, 2, 3, 4}; //circle

Line Loop (11) = {1, 2, 3, 4};
Line Loop (12) = {-5, -6, 7, 8, 9, 10};
Plane Surface (13) = {12, 11};
Physical Surface (14) = {13};

Field[1]=Attractor;
Field[1].EdgesList={3,4,8};

Field[2]=MathEval;
Field[2].F="0.005+0.5*F1";
Background Field=2;

View "comments" {T2(10, 15, 0){ StrCat("File created on ", Today) };};
