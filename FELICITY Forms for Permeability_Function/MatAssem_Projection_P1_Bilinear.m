function MATS = MatAssem_Projection_P1_Bilinear()
Omega = Domain('tetrahedron');
Ball = Domain('tetrahedron') < Omega;
Vector_P1 = Element(Ball,lagrange_deg1_dim3,3); 
vv = Test(Vector_P1);
uu = Trial(Vector_P1);
Mass_Matrix = Bilinear(Vector_P1,Vector_P1);
Mass_Matrix = Mass_Matrix + Integral(Ball,vv.val' * uu.val);
Quadrature_Order = 3;
G1 = GeoElement(Ball);
MATS = Matrices(Quadrature_Order,G1);
MATS = MATS.Append_Matrix(Mass_Matrix);
end