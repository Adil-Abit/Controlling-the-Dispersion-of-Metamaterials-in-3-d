function Right_Bilinear = MatAssem_Right_Bilinear()
Omega = Domain('tetrahedron');
Ball = Domain('tetrahedron') < Omega;
Vec_P1_Ball_Space = Element(Ball, lagrange_deg1_dim3,3);
P1_Ball_Space_test = Test(Vec_P1_Ball_Space);
P1_Ball_Space_trial = Trial(Vec_P1_Ball_Space);
Right_bi = Bilinear(Vec_P1_Ball_Space,Vec_P1_Ball_Space);
Right_bi = Right_bi + Integral(Ball,dot(P1_Ball_Space_test.val, P1_Ball_Space_trial.val));
Quadrature_Order = 3;
G1 = GeoElement(Omega);
Right_Bilinear = Matrices(Quadrature_Order,G1);
Right_Bilinear = Right_Bilinear.Append_Matrix(Right_bi);
end