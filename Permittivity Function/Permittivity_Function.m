function  E_effective = Permittivity_Function()
format long
RR =0.3; % The radius of plasmonic spherical inclusion; Needs to be modified for a different geometry.
rr =0.12; % The radius of high dielectric spherical inclusion; Needs to be modified for a different geometry.
Unit_Vectors = [1,0,0;0,1,0;0,0,1];
MY_MESH = load_gmsh2('con0.3(0.12)'); %For a different geometry, choose a different geo file here!
Vertex = MY_MESH.POS;
All_Tets = MY_MESH.TETS;
Elements_Without_R = All_Tets(1:16384,1:4);% This line(and the next 4 lines) need to be manually modified for  different geometry.
P_Elements= zeros(9216,4);
P_Elements(:) = All_Tets(7169:16384,1:4);
H_Elements= zeros(7168,4);
H_Elements(:) = All_Tets(1:7168,1:4);
Omega = Domain('tetrahedron');
Mesh_Without_R = MeshTetrahedron(Elements_Without_R, Vertex, 'Omega');

Mesh_Without_R = Mesh_Without_R.Append_Subdomain('3D','P_Region',P_Elements);
P_Mesh = Mesh_Without_R.Output_Subdomain_Mesh('P_Region');
Mesh_Without_R = Mesh_Without_R.Append_Subdomain('3D','H_Region',H_Elements);
H_Mesh = Mesh_Without_R.Output_Subdomain_Mesh('H_Region');
Mesh_Without_R = Mesh_Without_R.Remove_Unused_Vertices;
Mesh_Without_R = Mesh_Without_R.Reorder();
integration_domains = {'Omega';'P_Region';'H_Region'};
Subdomain_Embed1 = Mesh_Without_R.Generate_Subdomain_Embedding_Data(integration_domains);

V_ind = transpose(1:1:Mesh_Without_R.Num_Vtx); 
Vtx = Vertex;
P1 = uint32(Mesh_Without_R.ConnectivityList);
Element = lagrange_deg1_dim3();
REF = ReferenceFiniteElement(Element);
Lag_FES = FiniteElementSpace('LAG',REF,Mesh_Without_R,[]);
Main_Dir8 = '/Users/Almire/Documents/MATLAB';
Create_DoF_Allocator(Element,'Lag_DoF',Main_Dir8);
Lag_Space_DoFmap = Lag_DoF(uint32(Mesh_Without_R.ConnectivityList));
Lag_FES = Lag_FES.Set_DoFmap(Mesh_Without_R,Lag_Space_DoFmap);
Lag_FES_DoFs = Lag_FES.Get_DoFs;

Bdy_Faces = Mesh_Without_R.freeBoundary;

m=length(Bdy_Faces(:,1));
Center_Coordinates_Index_1 = zeros(m,1);
Center_Coordinates_Index_2 = zeros(m,1);
Center_Coordinates_Index_3 = zeros(m,1);
Center_Coordinates_Index_1(:)= Bdy_Faces(:,1);
Center_Coordinates_Index_2(:)= Bdy_Faces(:,2);
Center_Coordinates_Index_3(:)= Bdy_Faces(:,3);
Center_Coordinates_1 = Mesh_Without_R.Points(Center_Coordinates_Index_1,:);
Center_Coordinates_2 = Mesh_Without_R.Points(Center_Coordinates_Index_2,:);
Center_Coordinates_3 = Mesh_Without_R.Points(Center_Coordinates_Index_3,:);

Center_Coordinates = zeros(m,3);
Distance = zeros(m,1);
for n=1:m
    Center_Coordinates(n,1) = ( Center_Coordinates_1(n,1)+ Center_Coordinates_2(n,1)+ Center_Coordinates_3(n,1)   )/3;
    Center_Coordinates(n,2) = ( Center_Coordinates_1(n,2)+ Center_Coordinates_2(n,2)+ Center_Coordinates_3(n,2)   )/3;
    Center_Coordinates(n,3) = ( Center_Coordinates_1(n,3)+ Center_Coordinates_2(n,3)+ Center_Coordinates_3(n,3)   )/3;
    Distance(n) = ((0.5-Center_Coordinates(n,1))^2 + (0.5-Center_Coordinates(n,2))^2 +(0.5-Center_Coordinates(n,3))^2)^0.5;
end
Bdy_Logic  = ( Distance(:) < RR);
R_Bdy_Faces = Bdy_Faces(Bdy_Logic,:);
 
Mesh_Without_R = Mesh_Without_R.Append_Subdomain('2D','R_Bdy',R_Bdy_Faces);
R_Bdy_P1_DoFs = Lag_FES.Get_DoFs_On_Subdomain(Mesh_Without_R,'R_Bdy');
R_Bdy_XC = Mesh_Without_R.Points( R_Bdy_P1_DoFs,:);
Convert_Form_Definition_to_MEX(@MatAssem_Laplace_On_R,{},'mex_Laplace_On_R_assemble');
FEM = mex_Laplace_On_R_assemble([], Mesh_Without_R.Points,uint32(Mesh_Without_R.ConnectivityList),[],[],Lag_Space_DoFmap);
Soln = zeros(Mesh_Without_R.Num_Vtx,1);
Soln(R_Bdy_P1_DoFs(:)) = R_Bdy_XC(:,1);
A   = FEM(1).MAT;
RHS = 0*Soln;
RHS = RHS - A * Soln;
All_Vtx_Indices = (1:1:Mesh_Without_R.Num_Vtx)';
FreeNodes = setdiff(All_Vtx_Indices,R_Bdy_P1_DoFs);
Soln(FreeNodes,1) = A(FreeNodes,FreeNodes) \ RHS(FreeNodes,1); 

Convert_Form_Definition_to_MEX(@MatAssem_NEW_Left_H,{},'mex_NEW_Left_H_assemble');
Left_H_Matrix = mex_NEW_Left_H_assemble([],Mesh_Without_R.Points,P1,[],Subdomain_Embed1,Lag_Space_DoFmap);
H_Matrix = Left_H_Matrix.MAT;
Periodic_H_Matrix = Periodic_Matrix_Producer(V_ind,Vtx,H_Matrix);
Convert_Form_Definition_to_MEX(@MatAssem_NEW_Left_P,{},'mex_NEW_Left_P_assemble');
Left_P_Matrix = mex_NEW_Left_P_assemble([],Mesh_Without_R.Points,P1,[],Subdomain_Embed1,Lag_Space_DoFmap);
P_Matrix = Left_P_Matrix.MAT;
Periodic_P_Matrix = Periodic_Matrix_Producer(V_ind,Vtx,P_Matrix);
Convert_Form_Definition_to_MEX(@MatAssem_NEW_Right,{},'mex_NEW_Right_assemble');
Right_Matrix = mex_NEW_Right_assemble([],Mesh_Without_R.Points,P1,[],[],Lag_Space_DoFmap);
Right = Right_Matrix.MAT;
Left = 0.5*(H_Matrix-P_Matrix);
Periodic_Left = 0.5*(Periodic_H_Matrix - Periodic_P_Matrix);
[Periodic_Right,R_keep,left_right,bottom_top,front_back] = Periodic_Matrix_Producer(V_ind,Vtx,Right);
B = Periodic_Right;
DIM = size(B,1);
A = Periodic_Left;
[Q,Lambda_B] =eigs(B,DIM);
Tilde_A = Q'*A*Q;
Hat_A = Tilde_A(1:DIM-1,1:DIM-1);
Hat_Lambda_B = Lambda_B(1:DIM-1,1:DIM-1);
Final_Matrix = inv(Hat_Lambda_B)*Hat_A;
[Eigenvectors, Eigenvalues,left_eigenvectors] = eig(Final_Matrix);
Full_Rank_Eigenvecs = zeros(DIM,DIM);
Full_Rank_Eigenvecs(1:DIM-1,1:DIM-1)=Eigenvectors;
sup = ones(1,DIM-1);
Full_Rank_Eigenvecs(DIM,1:DIM-1)=sup';
last_column = Q'*ones(DIM,1);
Full_Rank_Eigenvecs(:,DIM)=last_column;
Eigenvecs_Returned_To_X = Q*Full_Rank_Eigenvecs;
Original_Eigenvecs = Eigenvecs_Returned_To_X;
full_eigenvalues = zeros(DIM,DIM);
full_eigenvalues(1:DIM-1,1:DIM-1) = Eigenvalues; 
full_eigenvalues(DIM,1:DIM) = 0;
full_eigenvalues(1:DIM,DIM) =0;

Eigenvecs = zeros(length(V_ind),DIM);  
Eigenvecs(R_keep,:) = Original_Eigenvecs;
Eigenvecs(front_back(:,2),:) = Eigenvecs(front_back(:,1),:);
Eigenvecs(bottom_top(:,2),:) = Eigenvecs(bottom_top(:,1),:);  
Eigenvecs(left_right(:,2),:) = Eigenvecs(left_right(:,1),:);
[Eigenvalues_diag,ind] = sort(real(diag(full_eigenvalues)));
Eigenvalues_sorted = full_eigenvalues(ind,ind);
Eigenvecs_sorted = Eigenvecs(:,ind); 

Convert_Form_Definition_to_MEX(@MatAssem_Grad_q_Grad_q_H,{},'mex_MatAssem_Grad_q_Grad_q_H_assemble');
Convert_Form_Definition_to_MEX(@MatAssem_Grad_q_Grad_q_P,{},'mex_MatAssem_Grad_q_Grad_q_P_assemble');
QQ_Matrix_H = mex_MatAssem_Grad_q_Grad_q_H_assemble([], Mesh_Without_R.Points,uint32(Mesh_Without_R.ConnectivityList),[],Subdomain_Embed1,Lag_Space_DoFmap);
QQ_Matrix_P = mex_MatAssem_Grad_q_Grad_q_P_assemble([], Mesh_Without_R.Points,uint32(Mesh_Without_R.ConnectivityList),[],Subdomain_Embed1,Lag_Space_DoFmap);
QQ_Integral_Matrix_H = QQ_Matrix_H(1).MAT; 
QQ_Integral_Matrix_P = QQ_Matrix_P(1).MAT;

Convert_Form_Definition_to_MEX(@MatAssem_Grad_q_Grad_psi_H,{},'mex_MatAssem_Grad_q_Grad_psi_H_assemble');
Convert_Form_Definition_to_MEX(@MatAssem_Grad_q_Grad_psi_P,{},'mex_MatAssem_Grad_q_Grad_psi_P_assemble');
Q_PSI_Matrix_H = mex_MatAssem_Grad_q_Grad_psi_H_assemble([], Mesh_Without_R.Points,uint32(Mesh_Without_R.ConnectivityList),[],Subdomain_Embed1,Lag_Space_DoFmap);
Q_PSI_Matrix_P = mex_MatAssem_Grad_q_Grad_psi_P_assemble([], Mesh_Without_R.Points,uint32(Mesh_Without_R.ConnectivityList),[],Subdomain_Embed1,Lag_Space_DoFmap);
Q_PSI_Integral_Matrix_H = Q_PSI_Matrix_H(1).MAT; 
Q_PSI_Integral_Matrix_P = Q_PSI_Matrix_P(1).MAT;

Convert_Form_Definition_to_MEX(@MatAssem_Grad_psi_H,{},'mex_MatAssem_Grad_psi_H_assemble');
Convert_Form_Definition_to_MEX(@MatAssem_Grad_psi_P,{},'mex_MatAssem_Grad_psi_P_assemble');
PSI_Matrix_H = mex_MatAssem_Grad_psi_H_assemble([], Mesh_Without_R.Points,uint32(Mesh_Without_R.ConnectivityList),[],Subdomain_Embed1,Lag_Space_DoFmap);
PSI_Matrix_P = mex_MatAssem_Grad_psi_P_assemble([], Mesh_Without_R.Points,uint32(Mesh_Without_R.ConnectivityList),[],Subdomain_Embed1,Lag_Space_DoFmap);
PSI_Integral_Matrix_H_1 = PSI_Matrix_H(1).MAT;
PSI_Integral_Matrix_H_2 = PSI_Matrix_H(2).MAT;
PSI_Integral_Matrix_H_3 = PSI_Matrix_H(3).MAT;
PSI_Integral_Matrix_P_1 = PSI_Matrix_P(1).MAT;
PSI_Integral_Matrix_P_2 = PSI_Matrix_P(2).MAT;
PSI_Integral_Matrix_P_3 = PSI_Matrix_P(3).MAT;


Eig_Vec_Num = 400; %% This (as well as the value of "Starting_Index" below) needs to be modified(simultaneously) if used a different geometry!
Starting_Index = 1700;  
Eigenvals = diag(Eigenvalues_sorted);
Lambda = Eigenvals(Starting_Index:Starting_Index+Eig_Vec_Num);
mu = 0.5 - Lambda;
Theta_P = (4/3)*pi*(RR)^3-(4/3)*pi*(rr)^3;
Theta_H = 1-(4/3)*pi*(RR)^3;
Theta_Y_Minus_R = 1- (4/3)*pi*(rr)^3;
ratio = 0.00001:0.000005:1;
step_num = 0.00001:0.000005:1;
E_effective = zeros(length(step_num), 1);
QQ_Integral_H = Soln'* QQ_Integral_Matrix_H *Soln;
QQ_Integral_P = Soln'* QQ_Integral_Matrix_P *Soln;
A_n_Numerator_First_Integral = zeros(length(Lambda),3);
A_n_Numerator_First_Integral(:,1) = (PSI_Integral_Matrix_H_1'* Eigenvecs_sorted(:,Starting_Index:Starting_Index+Eig_Vec_Num))';
A_n_Numerator_First_Integral(:,2) = (PSI_Integral_Matrix_H_2'* Eigenvecs_sorted(:,Starting_Index:Starting_Index+Eig_Vec_Num))';
A_n_Numerator_First_Integral(:,3) = (PSI_Integral_Matrix_H_3'* Eigenvecs_sorted(:,Starting_Index:Starting_Index+Eig_Vec_Num))';

A_n_Numerator_Second_Integral = zeros(length(Lambda),3);
A_n_Numerator_Second_Integral(:,1) = (PSI_Integral_Matrix_P_1'* Eigenvecs_sorted(:,Starting_Index:Starting_Index+Eig_Vec_Num))';
A_n_Numerator_Second_Integral(:,2) = (PSI_Integral_Matrix_P_2'* Eigenvecs_sorted(:,Starting_Index:Starting_Index+Eig_Vec_Num))';
A_n_Numerator_Second_Integral(:,3) = (PSI_Integral_Matrix_P_3'* Eigenvecs_sorted(:,Starting_Index:Starting_Index+Eig_Vec_Num))';

B_n_Numerator_First_Integral = zeros(length(Lambda),1);
B_n_Numerator_First_Integral(:) = Soln'* Q_PSI_Integral_Matrix_H *Eigenvecs_sorted(:,Starting_Index:Starting_Index+Eig_Vec_Num);
B_n_Numerator_Second_Integral = zeros(length(Lambda),1);
B_n_Numerator_Second_Integral(:) = Soln'* Q_PSI_Integral_Matrix_P *Eigenvecs_sorted(:,Starting_Index:Starting_Index+Eig_Vec_Num);

for  Count =1:length(step_num)
     E_p = 1-1/(ratio(Count));
     First_Term  = E_p*Theta_P+Theta_H;
     Sum = 0;
     A_n = zeros(length(Lambda),1);
     B_n = zeros(length(Lambda),1);
     for sum_index = 1:Eig_Vec_Num+1
         A_n_Denominator = 1- mu(sum_index)+ E_p*mu(sum_index);
         A_n(sum_index) =(-dot(Unit_Vectors(1,:),( A_n_Numerator_First_Integral(sum_index,:) +  E_p*A_n_Numerator_Second_Integral(sum_index,:))))/A_n_Denominator;
         B_n_Denominator =  A_n_Denominator;
         B_n(sum_index) =-( B_n_Numerator_First_Integral(sum_index) + E_p*B_n_Numerator_Second_Integral(sum_index))/ B_n_Denominator;
         Sum = Sum + (1+mu(sum_index)*(E_p-1))*( A_n(sum_index)+B_n(sum_index))^2;
     end

     E_effective(Count) =  First_Term + Theta_Y_Minus_R*Sum + QQ_Integral_H + E_p* QQ_Integral_P;

end



end
