function MATS = MatAssem_Laplace_On_R()
Omega = Domain('tetrahedron');
P1_Space = Element(Omega,lagrange_deg1_dim3,1);
v = Test(P1_Space);
u = Trial(P1_Space);
Stiff_Matrix = Bilinear(P1_Space,P1_Space);
Stiff_Matrix = Stiff_Matrix + Integral(Omega, v.grad' * u.grad );
Quadrature_Order = 1;
G1 = GeoElement(Omega);
MATS = Matrices(Quadrature_Order,G1);
MATS = MATS.Append_Matrix(Stiff_Matrix);
end