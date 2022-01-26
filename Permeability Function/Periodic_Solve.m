function  [U, left_right,bottom_top,front_back] = Periodic_Solve(V_ind,Vtx,A,R,M)

[left_right,bottom_top,front_back] = new_aligner(Vtx, V_ind );
A_p = A;
m = length(V_ind);
twom = 2*m;

A_p(left_right(:,1),:) = A_p(left_right(:,1),:)+ A_p(left_right(:,2),:);
A_p(left_right(:,1)+m,m:twom) = A_p(left_right(:,1)+m,m:twom)+ A_p(left_right(:,2)+m,m:twom);
A_p(left_right(:,1)+twom,twom:end) = A_p(left_right(:,1)+twom,twom:end)+ A_p(left_right(:,2)+twom,twom:end);


A_p(:,left_right(:,1)) = A_p(:,left_right(:,1))+ A_p(:,left_right(:,2));
A_p(m:twom,left_right(:,1)+m) = A_p(m:twom,left_right(:,1)+m)+ A_p(m:twom,left_right(:,2)+m);
A_p(twom:end,left_right(:,1)+twom) = A_p(twom:end,left_right(:,1)+twom)+ A_p(twom:end,left_right(:,2)+twom);


A_p(bottom_top(:,1),:) = A_p(bottom_top(:,1),:)+ A_p(bottom_top(:,2),:);
A_p(bottom_top(:,1)+m,m:twom) = A_p(bottom_top(:,1)+m,m:twom)+ A_p(bottom_top(:,2)+m,m:twom);
A_p(bottom_top(:,1)+twom,twom:end) = A_p(bottom_top(:,1)+twom,twom:end)+ A_p(bottom_top(:,2)+twom,twom:end);


A_p(:,bottom_top(:,1)) = A_p(:,bottom_top(:,1))+ A_p(:,bottom_top(:,2));
A_p(m:twom,bottom_top(:,1)+m) = A_p(m:twom,bottom_top(:,1)+m)+ A_p(m:twom,bottom_top(:,2)+m);
A_p(twom:end,bottom_top(:,1)+twom) = A_p(twom:end,bottom_top(:,1)+twom)+ A_p(twom:end,bottom_top(:,2)+twom);


A_p(front_back(:,1),:) = A_p(front_back(:,1),:)+ A_p(front_back(:,2),:);
A_p(front_back(:,1)+m,m:twom) = A_p(front_back(:,1)+m,m:twom)+ A_p(front_back(:,2)+m,m:twom);
A_p(front_back(:,1)+twom,twom:end) = A_p(front_back(:,1)+twom,twom:end)+ A_p(front_back(:,2)+twom,twom:end);


A_p(:,front_back(:,1)) = A_p(:,front_back(:,1))+ A_p(:,front_back(:,2));
A_p(m:twom,front_back(:,1)+m) = A_p(m:twom,front_back(:,1)+m)+ A_p(m:twom,front_back(:,2)+m);
A_p(twom:end,front_back(:,1)+twom) = A_p(twom:end,front_back(:,1)+twom)+ A_p(twom:end,front_back(:,2)+twom);


RHS = R;
RHS(left_right(:,1)) = RHS(left_right(:,1)) + RHS(left_right(:,2));
RHS(left_right(:,1)+m) = RHS(left_right(:,1)+m) + RHS(left_right(:,2)+m);
RHS(left_right(:,1)+twom) = RHS(left_right(:,1)+twom) + RHS(left_right(:,2)+twom);

RHS(bottom_top(:,1)) = RHS(bottom_top(:,1)) + RHS(bottom_top(:,2));
RHS(bottom_top(:,1)+m) = RHS(bottom_top(:,1)+m) + RHS(bottom_top(:,2)+m);
RHS(bottom_top(:,1)+twom) = RHS(bottom_top(:,1)+twom) + RHS(bottom_top(:,2)+twom);

RHS(front_back(:,1)) = RHS(front_back(:,1)) + RHS(front_back(:,2));
RHS(front_back(:,1)+m) = RHS(front_back(:,1)+m) + RHS(front_back(:,2)+m);
RHS(front_back(:,1)+twom) = RHS(front_back(:,1)+twom) + RHS(front_back(:,2)+twom);



R_keep_one = setdiff(V_ind, union(left_right(:,2),bottom_top(:,2)));
R_keep = setdiff( R_keep_one,front_back(:,2) );
Free_nodes = setdiff(R_keep, R_keep(end));



x_block  = A_p(1:m,1:m);
y_block = A_p(m+1:twom,m+1:twom);
z_block  = A_p(twom+1:end,twom+1:end);

A_x = x_block(Free_nodes, Free_nodes);
A_y = y_block(Free_nodes, Free_nodes);
A_z = z_block(Free_nodes, Free_nodes);


RHS_x_block = RHS(1:m,1);
RHS_y_block = RHS(m+1:twom,1);
RHS_z_block = RHS(twom+1:end,1);


RHSx = RHS_x_block(Free_nodes);
RHSy = RHS_y_block(Free_nodes);
RHSz = RHS_z_block(Free_nodes);


Solnx = zeros(m,1);
Solny = zeros(m,1);
Solnz = zeros(m,1);

disp(' ');
TOL = 1e-8;
disp(['Solve linear system with AGMG with TOL = ', num2str(TOL,'%1.2G'), ':']);
Soln_AGMGx = Solnx;
tic
Soln_AGMGx(Free_nodes,1) = agmg(A_x,RHSx,[],TOL);
toc
Soln_AGMGy = Solny;
tic
Soln_AGMGy(Free_nodes,1) = agmg(A_y,RHSy,[],TOL);
toc
Soln_AGMGz = Solnz;
tic
Soln_AGMGz(Free_nodes,1) = agmg(A_z,RHSz,[],TOL);
toc
Solnx = Soln_AGMGx;
Solny = Soln_AGMGy;
Solnz = Soln_AGMGz;


Solnx(front_back(:,2)) = Solnx(front_back(:,1));
Solnx(bottom_top(:,2)) = Solnx(bottom_top(:,1));
Solnx(left_right(:,2)) = Solnx(left_right(:,1));



Solny(front_back(:,2)) = Solny(front_back(:,1));
Solny(bottom_top(:,2)) = Solny(bottom_top(:,1));
Solny(left_right(:,2)) = Solny(left_right(:,1));

Solnz(front_back(:,2)) = Solnz(front_back(:,1));
Solnz(bottom_top(:,2)) = Solnz(bottom_top(:,1));
Solnz(left_right(:,2)) = Solnz(left_right(:,1));


Soln_ghost_x = zeros(3*m,1);
Soln_ghost_x(V_ind) = Solnx(V_ind);
Ones = ones(1,3*m);
Ux = Solnx - Ones*M*Soln_ghost_x;
 

Soln_ghost_y = zeros(3*m,1);
Soln_ghost_y(V_ind) = Solny(V_ind);
Ones = ones(1,3*m);
Uy = Solny - Ones*M*Soln_ghost_y;
 

Soln_ghost_z = zeros(3*m,1);
Soln_ghost_z(V_ind) = Solnz(V_ind);
Ones = ones(1,3*m);
Uz = Solnz - Ones*M*Soln_ghost_z;
 
U = [Ux; Uy; Uz];
end

