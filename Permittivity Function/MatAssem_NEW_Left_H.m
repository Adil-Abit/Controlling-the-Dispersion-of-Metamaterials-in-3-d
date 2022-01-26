function Left_H_Matrix = MatAssem_NEW_Left_H()
Omega = Domain('tetrahedron');
H_Region = Domain('tetrahedron')< Omega;
Lagrange_P1 = Element(Omega,lagrange_deg1_dim3,1); 
vv = Test(Lagrange_P1);
uu = Trial(Lagrange_P1);
Left_H = Bilinear(Lagrange_P1,Lagrange_P1);
Left_H =Left_H + Integral(H_Region,sum(vv.grad .* uu.grad));
Quadrature_Order = 3;
G1 = GeoElement(Omega);
Left_H_Matrix = Matrices(Quadrature_Order,G1);
Left_H_Matrix =Left_H_Matrix.Append_Matrix(Left_H);
end