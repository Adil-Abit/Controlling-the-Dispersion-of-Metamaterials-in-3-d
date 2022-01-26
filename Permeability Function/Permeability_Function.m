function  MU_Effective = Permeability_Function()
format long
[Mesh,All_Tets, Vertex] = All_Mesh4();
P1_Space_DoFmap = uint32(Mesh.ConnectivityList); 
Ball_Elements= zeros(2688,4); % This line(and the next line) need to be manually modified for  different geometry.
Ball_Elements(:) = All_Tets(16385:19072,1:4);
Mesh = Mesh.Append_Subdomain('3D','Ball',Ball_Elements);
Ball_Mesh = Mesh.Output_Subdomain_Mesh('Ball');
integration_domains = {'Omega'; 'Ball'}; 
Subdomain_Embed = Mesh.Generate_Subdomain_Embedding_Data(integration_domains); 

Element = lagrange_deg1_dim3();
REF = ReferenceFiniteElement(Element);
Lag_FES = FiniteElementSpace('LAG',REF,Mesh,[]);
Main_Dir11 = '/Users/Almire/Documents/MATLAB'; % This needs to be  the local MATLAB folder  of your computer.
Create_DoF_Allocator(Element,'Lag_DoF',Main_Dir11);
Lag_Space_DoFmap = Lag_DoF(uint32(Mesh.ConnectivityList));
Lag_FES = Lag_FES.Set_DoFmap(Mesh,Lag_Space_DoFmap);
Lag_FES_DoFs = Lag_FES.Get_DoFs;
Ball_P1_DoFs = Lag_FES.Get_DoFs_On_Subdomain(Mesh,'Ball');
Ball_P1_Space_DoFmap = uint32(Ball_Mesh.ConnectivityList);

Elem    =raviart_thomas_deg0_dim3();
RE = ReferenceFiniteElement(Elem);
Ball_RT_FES = FiniteElementSpace('RT',RE,Ball_Mesh,[]);
Main_Dir00 = '/Users/Almire/Desktop';
Create_DoF_Allocator(Elem,'RT_DoF',Main_Dir00);
Ball_RT_Space_DoFmap = RT_DoF(uint32(Ball_Mesh.ConnectivityList));
Ball_RT_FES = Ball_RT_FES.Set_DoFmap(Ball_Mesh,Ball_RT_Space_DoFmap);
Ball_RT_DoFs = Ball_RT_FES.Get_DoFs;
Ball_Bdy_Faces = Ball_Mesh.freeBoundary; 
Ball_Mesh = Ball_Mesh.Append_Subdomain('2D','Ball_Bdy',Ball_Bdy_Faces);
Ball_RT_Bdy_DoFs = Ball_RT_FES.Get_DoFs_On_Subdomain(Ball_Mesh,'Ball_Bdy'); 
Ball_Internal_DoFs = setdiff(Ball_RT_DoFs,Ball_RT_Bdy_DoFs);
[~, Ball_Orient] = Mesh.Get_Facet_Info(Ball_Bdy_Faces); 
[~, Ball_FE_Orient] = Ball_Mesh.Get_Facet_Info(Ball_Bdy_Faces);

Convert_Form_Definition_to_MEX(@MatAssem_Eig_Bilinear_Form_div_new,{},'mex_Eig_Bilinear_Form_div_new_assemble');
Bilinear = mex_Eig_Bilinear_Form_div_new_assemble([],Ball_Mesh.Points,Ball_P1_Space_DoFmap,Ball_FE_Orient,[],Ball_RT_Space_DoFmap);
LL1 = Bilinear.MAT;
LL = LL1 (Ball_Internal_DoFs,Ball_Internal_DoFs);
[V,D] =eigs(LL,length(Ball_Internal_DoFs)); 
leng =0; 
for zer = 1:length(Ball_Internal_DoFs)
    if D(zer,zer)<10^(-7)
        leng =leng+1;
    end
end
RT_basis_vectors = V(:,length(Ball_Internal_DoFs)-leng+1:end);
N  = length(Ball_RT_DoFs);
W = zeros(N, leng);
W(Ball_Internal_DoFs,:) = RT_basis_vectors;

Convert_Form_Definition_to_MEX(@MatAssem_Projection_P1_Bilinear,{},'mex_Projection_P1_Bilinear_assemble');
Projection_P1_Bilinear = mex_Projection_P1_Bilinear_assemble([],Ball_Mesh.Points,Ball_P1_Space_DoFmap,[],[],Ball_P1_Space_DoFmap);
P1_Psi_Psi = Projection_P1_Bilinear(1).MAT;

Convert_Form_Definition_to_MEX(@MatAssem_Projection_B_For_RHS_Blinear,{},'mex_Projection_B_For_RHS_Blinear_assemble');
Projection_P1_RT_Bilinear = mex_Projection_B_For_RHS_Blinear_assemble([],Ball_Mesh.Points,Ball_P1_Space_DoFmap,Ball_FE_Orient,[],Ball_RT_Space_DoFmap,Ball_P1_Space_DoFmap);
P1_RT_For_B = Projection_P1_RT_Bilinear(1).MAT;

T = P1_RT_For_B*W;
UU =zeros(length(P1_Psi_Psi(1,:)), length(T(1,:)));
for j= 1:length(T(1,:))
UU(:,j) = P1_Psi_Psi \ T(:,j);
end
[lin_independent_cols,idx]=licols(UU,1e-10);

Convert_Form_Definition_to_MEX(@MatAssem_Poisson_3D_Vector,{},'mex_Poisson_3D_Vector_assemble');
FEM = mex_Poisson_3D_Vector_assemble([],Mesh.Points,P1_Space_DoFmap,[],[],P1_Space_DoFmap);
V_ind = transpose(1:1:Mesh.Num_Vtx); 
Vtx = Vertex;
A = FEM(2).MAT;
M = FEM(1).MAT;

Convert_Form_Definition_to_MEX(@MatAssem_Poisson_3D_Vector_Linear_Form,{},'mex_Poisson_3D_Vector_Linear_Form_assemble');
for j =1:length(lin_independent_cols(1,:))
    tt = length(lin_independent_cols(:,1))/3;
    Lv = zeros(tt,3);

    for r=1:tt
        Lv(r,1) = lin_independent_cols(r,j);
        Lv(r,2) = lin_independent_cols(tt+r,j);
        Lv(r,3) = lin_independent_cols(2*tt+r,j);
    end
    LINEAR = mex_Poisson_3D_Vector_Linear_Form_assemble([],Mesh.Points,P1_Space_DoFmap,[],Subdomain_Embed,Ball_P1_Space_DoFmap,P1_Space_DoFmap,Lv);
    R = LINEAR(1).MAT;
    [U, Ux, Uy, Uz] = Periodic_Solve(V_ind,Vtx,A,R,M);
    S(j).U = U;
end
eigen_prob_first_term = zeros(length(lin_independent_cols(1,:)),length(lin_independent_cols(1,:)));
for n=1:length(lin_independent_cols(1,:))
    for t = 1:length(lin_independent_cols(1,:))
    eigen_prob_first_term(n,t) = S(n).U'*A*S(t).U;
    end
end
psi = zeros(length(A(1,:)),length(lin_independent_cols(1,:)));
for j=1:length(lin_independent_cols(1,:))
    psi_carrier = S(j).U;
    psi(1:length(A(1,:)), j) = psi_carrier(1:length(A(1,:)),1);
end

Convert_Form_Definition_to_MEX(@MatAssem_Eig_Left_Cross_Product_Linear_Form_With_P1Basis,{},'mex_Eig_Left_Cross_Product_Linear_Form_With_P1Basis_assemble');
x_cross_f_Linear = mex_Eig_Left_Cross_Product_Linear_Form_With_P1Basis_assemble([],Mesh.Points,P1_Space_DoFmap,[],Subdomain_Embed,Ball_P1_Space_DoFmap);
Chi_First = x_cross_f_Linear(1).MAT;
Chi_Second = x_cross_f_Linear(2).MAT;
Chi_Third = x_cross_f_Linear(3).MAT;
Second_Term_System_Matrix = Chi_First*Chi_First' + Chi_Second*Chi_Second' + Chi_Third*Chi_Third';
eigen_prob_second_term = lin_independent_cols'*Second_Term_System_Matrix*lin_independent_cols;
eig_prob_left_matrix_non_symmetric = eigen_prob_first_term + 0.25* eigen_prob_second_term;

Convert_Form_Definition_to_MEX(@MatAssem_Right_Bilinear,{},'mex_Right_Bilinear_assemble');
Right_Bilinear = mex_Right_Bilinear_assemble([],Mesh.Points,P1_Space_DoFmap,[],Subdomain_Embed,Ball_P1_Space_DoFmap);
original_right_matrix =Right_Bilinear.MAT;
Right_Matrix = lin_independent_cols'*original_right_matrix*lin_independent_cols;
eig_prob_right_matrix_non_symmetric = Right_Matrix;
eig_prob_left_matrix = 0.5*(eig_prob_left_matrix_non_symmetric'+eig_prob_left_matrix_non_symmetric);
eig_prob_right_matrix = 0.5*(eig_prob_right_matrix_non_symmetric'+eig_prob_right_matrix_non_symmetric);

Convert_Form_Definition_to_MEX(@MatAssem_Mu_Effective_Linear_Form_With_P1Basis,{},'mex_Mu_Effective_Linear_Form_With_P1Basis_assemble');
Mu_Matrices = mex_Mu_Effective_Linear_Form_With_P1Basis_assemble([],Mesh.Points,P1_Space_DoFmap,[],Subdomain_Embed,Ball_P1_Space_DoFmap);
First_Component = Mu_Matrices(1).MAT;
Second_Component = Mu_Matrices(2).MAT;
Third_Component = Mu_Matrices(3).MAT;
[Eigenvectors,Eigenvalues] =eig(eig_prob_left_matrix, eig_prob_right_matrix);
Eigenvectors = real(Eigenvectors);
Eigenvalues = real(Eigenvalues);

num_eig_vec = 3; % This (as well as the value of "s" below) needs to be modified(simultaneously) if used a different geometry!
Mu_Factor = zeros(num_eig_vec,1);
para = length(Eigenvectors(:,1));
s=1906;
Mu_effective_eigenvectors = zeros(length(lin_independent_cols(:,1)),num_eig_vec);
for filling_column_index = 1:num_eig_vec
    for filling_row_index = 1:para
        Mu_effective_eigenvectors(filling_row_index,filling_column_index) = Eigenvectors(filling_row_index,s+filling_column_index);
    end
end
for num_vector = 1:num_eig_vec
G = zeros(3,1);
G(1) = dot(First_Component, Mu_effective_eigenvectors(:,num_vector));
G(2) = dot(Second_Component, Mu_effective_eigenvectors(:,num_vector));
G(3) = dot(Third_Component, Mu_effective_eigenvectors(:,num_vector));
Mu_Factor(num_vector,1) = sqrt(G(1)^2+G(2)^2+G(3)^2);
end

First_Hundred_Eigenvalues = zeros(num_eig_vec,1);
for nn = 1:num_eig_vec
    First_Hundred_Eigenvalues(nn) = Eigenvalues(s+nn,s+nn);
end
step_num = 0.00001:0.000005:1; 
MU_Effective = zeros(length(step_num), 1);
for  h =1:length(step_num)
     coeff = 0;
     Sum_Term = zeros(num_eig_vec,1);
     Sum_Term_All = 0;
     for mm = 1:num_eig_vec
     coeff = (step_num(h))/(1/(3e15)-First_Hundred_Eigenvalues(mm)*step_num(h));
     Sum_Term(mm) =  coeff*Mu_Factor(mm);
     Sum_Term_All = Sum_Term_All + Sum_Term(mm);
     end
     MU = (0.25*Sum_Term_All)/3 + 1;
     MU_Effective(h) = MU;
end









end