nk__lc = 0.005;

// Define the points
Point(1) = {-1.0,  0,   0, nk__lc};
Point(2) = {   0,  0,   0, nk__lc};
Point(3) = {   0,  1.0, 0, nk__lc};
Point(4) = {-1.0,  1.0, 0, nk__lc};
Point(5) = { 1.0,  0,   0, nk__lc};
Point(6) = { 1.0,  1.0, 0, nk__lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 2};
Line(6) = {6, 5};
Line(7) = {3, 6};

Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {2, 7, 6, 5};

// Define the rectangles
Plane Surface(1) = {1};
Plane Surface(2) = {2};

// Define the physical surfaces (optional)
Physical Surface("Rectangle1") = {1};
Physical Surface("Rectangle2") = {2};

