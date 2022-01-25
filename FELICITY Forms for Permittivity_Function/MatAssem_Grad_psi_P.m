function PSI_Matrix_P = MatAssem_Grad_psi_P()
Omega = Domain('tetrahedron');
P_Region = Domain('tetrahedron') < Omega;
P1_Space = Element(Omega,lagrange_deg1_dim3,1);
u = Test(P1_Space);
first_component = Linear(P1_Space);
second_component = Linear(P1_Space);
third_component = Linear(P1_Space);
Grad_U = u.grad;

first_component = first_component + Integral(P_Region,Grad_U(1));
second_component = second_component + Integral(P_Region,Grad_U(2));
third_component =third_component + Integral(P_Region,Grad_U(3));

Quadrature_Order = 3;
G1 = GeoElement(Omega);
PSI_Matrix_P = Matrices(Quadrature_Order,G1);
PSI_Matrix_P = PSI_Matrix_P.Append_Matrix(first_component);
PSI_Matrix_P = PSI_Matrix_P.Append_Matrix(second_component);
PSI_Matrix_P = PSI_Matrix_P.Append_Matrix(third_component);
end