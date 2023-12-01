//+ Building example for shading masks

SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 2, 2, 2};
//+
Box(2) = {5, 0, 0, 3, 2, 5};
//+
Physical Surface("Bat1_face_Sud", 25) = {3};
//+
Physical Surface("Bat2_face_Sud", 26) = {9};
//+
Physical Surface("Bat1_face_Ouest", 27) = {1};
//+
Physical Surface("Bat2_face_Ouest", 28) = {7};
//+
Physical Surface("Bat2_face_Nord", 29) = {10};
//+
Physical Surface("Bat1_face_Nord", 30) = {4};
//+
Physical Surface("Bat1_face_Est", 31) = {2};
//+
Physical Surface("Bat2_face_Est", 32) = {8};
//+
Physical Surface("Bat2_face_Top", 33) = {12};
//+
Physical Surface("Bat1_face_Top", 34) = {6};
//+
// Physical Surface("Bat1_face_Bottom", 35) = {5};
//+
// Physical Surface("Bat2_face_Bottom", 36) = {11};
//+
Physical Surface("Bat1") = {3,1,4,2,6,5};
//+
Physical Surface("Bat2") = {7,8,9,10,12,11};
Delete {
  Volume{1}; 
}
//+
Delete {
  Volume{2}; 
}
