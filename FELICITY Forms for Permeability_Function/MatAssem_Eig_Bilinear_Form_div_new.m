function Right_Matrix = MatAssem_Eig_Bilinear_Form_div_new()
Omega = Domain('tetrahedron');
Ball = Domain('tetrahedron') < Omega;
RT_Space = Element(Ball,raviart_thomas_deg0_dim3);
RT_test = Test(RT_Space);
RT_trial = Trial(RT_Space);
Right = Bilinear(RT_Space,RT_Space);
Right = Right + Integral(Ball,(RT_test.div)*( RT_trial.div));
Quadrature_Order = 3;
G1 = GeoElement(Ball);
Right_Matrix = Matrices(Quadrature_Order,G1);
Right_Matrix = Right_Matrix.Append_Matrix(Right);
end