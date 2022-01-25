function Q_PSI_Matrix_P = MatAssem_Grad_q_Grad_psi_P()
Omega = Domain('tetrahedron');
P_Region = Domain('tetrahedron') < Omega;
P1_Space = Element(Omega,lagrange_deg1_dim3,1);
v = Test(P1_Space);
u = Trial(P1_Space);
q_psi_Matrix_P = Bilinear(P1_Space,P1_Space);
q_psi_Matrix_P = q_psi_Matrix_P + Integral(P_Region, v.grad' * u.grad );
Quadrature_Order = 1;
G1 = GeoElement(Omega);
MATS = Matrices(Quadrature_Order,G1);
Q_PSI_Matrix_P = MATS.Append_Matrix(q_psi_Matrix_P);
end