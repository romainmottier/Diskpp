squareSize = 4600.;
gridsize = squareSize/200;
air_thick = 600.;
solid_thick = squareSize-air_thick;

Point(1) = {-squareSize/2.,-solid_thick,0.,gridsize};
Point(2) = {squareSize/2.,-solid_thick,0.,gridsize};
Point(3) = {squareSize/2.,0.,0.,gridsize};
Point(4) = {-squareSize/2.,0.,0.,gridsize};
Point(5) = {-squareSize/2.,air_thick,0.,gridsize};
Point(6) = {squareSize/2.,air_thick,0.,gridsize};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {4,5};
Line(6) = {5,6};
Line(7) = {6,3};

Line Loop(1) = {1,2,3,4};
Line Loop(2) = {5,6,7,3};

Plane Surface(1)={1};
Plane Surface(2)={2};

//Transfinite Surface {0,1};
Transfinite Surface {1,2};
Recombine Surface "*";

Physical Surface(1) = {1};
Physical Surface(2) = {2};
