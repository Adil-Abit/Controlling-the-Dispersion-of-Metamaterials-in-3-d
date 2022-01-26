function MATS = MatAssem_Projection_B_For_RHS_Blinear()
Omega = Domain('tetrahedron');
Ball = Domain('tetrahedron') < Omega;
Vector_P1 = Element(Ball,lagrange_deg1_dim3,3);
RT_Space = Element(Ball,raviart_thomas_deg0_dim3);
vv = Test(Vector_P1);
uu = Trial(RT_Space);
Matrix = Bilinear(Vector_P1,RT_Space);
Matrix =Matrix + Integral(Ball,vv.val' * uu.val);
Quadrature_Order = 3;
G1 = GeoElement(Ball);
MATS = Matrices(Quadrature_Order,G1);
MATS = MATS.Append_Matrix(Matrix);
end