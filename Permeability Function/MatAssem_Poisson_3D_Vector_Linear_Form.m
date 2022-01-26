function LINEAR = MatAssem_Poisson_3D_Vector_Linear_Form()
Omega = Domain('tetrahedron');
Ball = Domain('tetrahedron') < Omega;
Vector_P1 = Element(Omega,lagrange_deg1_dim3,3); 
vv = Test(Vector_P1);
Vec_P1_Ball_Space = Element(Ball, lagrange_deg1_dim3,3);
Lv = Coef(Vec_P1_Ball_Space);
RHS = Linear(Vector_P1);
RHS = RHS + Integral(Ball,sum(Lv.val'* vv.val));
Quadrature_Order = 3;
G1 = GeoElement(Omega);
LINEAR = Matrices(Quadrature_Order,G1);
LINEAR = LINEAR.Append_Matrix(RHS);
end