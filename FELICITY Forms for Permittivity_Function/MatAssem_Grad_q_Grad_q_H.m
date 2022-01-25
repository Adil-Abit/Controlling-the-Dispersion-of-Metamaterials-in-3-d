function QQ_Matrix_H = MatAssem_Grad_q_Grad_q_H()
Omega = Domain('tetrahedron');
H_Region = Domain('tetrahedron') < Omega;
P1_Space = Element(Omega,lagrange_deg1_dim3,1);
v = Test(P1_Space);
u = Trial(P1_Space);
q_q_Matrix_H = Bilinear(P1_Space,P1_Space);
q_q_Matrix_H = q_q_Matrix_H + Integral(H_Region, v.grad' * u.grad );
Quadrature_Order = 1;
G1 = GeoElement(Omega);
MATS = Matrices(Quadrature_Order,G1);
QQ_Matrix_H = MATS.Append_Matrix(q_q_Matrix_H);
end