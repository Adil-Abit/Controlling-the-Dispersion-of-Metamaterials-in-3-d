function MATS = MatAssem_Poisson_3D_Vector()
Omega = Domain('tetrahedron');
Vector_P1 = Element(Omega,lagrange_deg1_dim3,3); 
vv = Test(Vector_P1);
uu = Trial(Vector_P1);
Mass_Matrix = Bilinear(Vector_P1,Vector_P1);
Mass_Matrix = Mass_Matrix + Integral(Omega,vv.val' * uu.val);
Stiff_Matrix = Bilinear(Vector_P1,Vector_P1);
Stiff_Matrix = Stiff_Matrix + Integral(Omega,sum(sum(vv.grad .* uu.grad)));
Quadrature_Order = 3;
G1 = GeoElement(Omega);
MATS = Matrices(Quadrature_Order,G1);
MATS = MATS.Append_Matrix(Stiff_Matrix);
MATS = MATS.Append_Matrix(Mass_Matrix);
end