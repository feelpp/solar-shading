//+ Building example for shading masks

SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 2, 2, 2};
//+
Box(2) = {5, 0, 0, 3, 2, 5};
//+
Physical Surface("Bat1Sud", 25) = {3};
//+
Physical Surface("Bat2Sud", 26) = {9};
//+
Physical Surface("Bat1Ouest", 27) = {1};
//+
Physical Surface("Bat2Ouest", 28) = {7};
//+
Physical Surface("Bat2Nord", 29) = {10};
//+
Physical Surface("Bat1Nord", 30) = {4};
//+
Physical Surface("Bat1Est", 31) = {2};
//+
Physical Surface("Bat2Est", 32) = {8};
//+
Physical Surface("Bat2Top", 33) = {12};
//+
Physical Surface("Bat1Top", 34) = {6};
//+
// Physical Surface("Bat1Bottom", 35) = {5};
//+
// Physical Surface("Bat2Bottom", 36) = {11};
//+
Physical Volume("Bat1") = {1};
//+
Physical Volume("Bat2") = {2};
