function Left_P_Matrix = MatAssem_NEW_Left_P()
Omega = Domain('tetrahedron');
P_Region = Domain('tetrahedron')< Omega;
Lagrange_P1 = Element(Omega,lagrange_deg1_dim3,1); 
vv = Test(Lagrange_P1);
uu = Trial(Lagrange_P1);
Left_P = Bilinear(Lagrange_P1,Lagrange_P1);
Left_P =Left_P + Integral(P_Region,sum(vv.grad .* uu.grad));
Quadrature_Order = 3;
G1 = GeoElement(Omega);
Left_P_Matrix = Matrices(Quadrature_Order,G1);
Left_P_Matrix =Left_P_Matrix.Append_Matrix(Left_P);
end