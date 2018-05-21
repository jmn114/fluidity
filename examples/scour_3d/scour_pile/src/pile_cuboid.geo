lc = 1e-1;
D = 0.5;

xSize = 5.66+1.5; 
ySize = 3.60; 
zSize = -0.4;

Point(1) = {0,0,0};
Point(2) = {0,ySize,0};
Point(3) = {xSize,ySize,0};
Point(4) = {xSize,0,0};

Point(5) = {3.85,1.8,0};
Point(6) = {3.85-D/2,1.8,0};
Point(7) = {3.85,1.8+D/2,0};
Point(8) = {3.85+D/2,1.8,0};
Point(9) = {3.85,1.8-D/2,0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Circle(5) = {6,5,7};
Circle(6) = {7,5,8};
Circle(7) = {8,5,9};
Circle(8) = {9,5,6};

Line Loop(1) = {5,6,7,8}; // interior loop
Line Loop(2) = {1,2,3,4}; // exterior loop
//Plane Surface(1) = {1}; // interior surface
Plane Surface(1) = {2,1}; // exterior surface (with a whole)

Physical Surface(1) = {21}; // inlet
Physical Surface(2) = {29}; // oulet
Physical Surface(3) = {33,25}; // walls
Physical Surface(4) = {1}; // top
Physical Surface(5) = {50}; // bottom
Physical Surface(6) = {37,41,45,49}; //pile

Extrude {0,0,zSize} {
  Surface{1};
}

// Define Volume
Physical Volume(1) = {1};  

// refine mesh locally using attractor field 
Field[1] = Attractor; 
//Field[1].NodesList = {6,7,8,9,26,38,33,28}; 
Field[1].FacesList = {37,41,45,49}; 

Field[2] = MathEval;
Field[2].F = "0.1";

//Field[2] = Threshold; 
//Field[2].IField = 1; 
//Field[2].LcMin = lc / 5; 
//Field[2].LcMin = lc; 
//Field[2].LcMax = 10*lc; 
//Field[2].DistMin = lc / 2;
//Field[2].DistMin = 0.1;  
//Field[2].DistMax = 2 * lc;
//Field[2].DistMax = 0.5;  

// Use minimum of all the fields as the background field 
Field[3] = Max; 
Field[3].FieldsList = {2,1}; 
//Background Field = 3; 

// Don't extend the elements sizes from the boundary inside the domain 

Mesh.CharacteristicLengthExtendFromBoundary = 0;

Field[4] = MathEvalAniso;
Field[4].m12 = "0.0";
Field[4].m13 = "0.0";
Field[4].m23 = "0.0";
//Field[4].m11 = "1/(0.01+0.1*(F3-0.1))^2";
//Field[4].m22 = "1/(0.01+0.1*(F3-0.1))^2";
Field[4].m11 = "1/(0.08+0.1*(F3-0.1))^2";
Field[4].m22 = "1/(0.08+0.1*(F3-0.1))^2";
Field[4].m33 = "1/(0.1+0.2*z)^2";

Background Field = 4;

Mesh.Algorithm = 7;
Mesh.Algorithm3D = 7;

