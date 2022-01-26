function  [Periodic_Matrix,R_keep,left_right,bottom_top,front_back, LR,BT,FB] = Periodic_Matrix_Producer(V_ind,Vtx,testleft)

[left_right,bottom_top,front_back] = new_aligner(Vtx, V_ind );

Left_p = testleft;

Left_p(left_right(:,1),:) = Left_p(left_right(:,1),:)+ Left_p(left_right(:,2),:);
Left_p(:,left_right(:,1)) = Left_p(:,left_right(:,1))+ Left_p(:,left_right(:,2));

Left_p(bottom_top(:,1),:) = Left_p(bottom_top(:,1),:)+ Left_p(bottom_top(:,2),:);
Left_p(:,bottom_top(:,1)) = Left_p(:,bottom_top(:,1))+ Left_p(:,bottom_top(:,2));
 
Left_p(front_back(:,1),:) = Left_p(front_back(:,1),:)+ Left_p(front_back(:,2),:);
Left_p(:,front_back(:,1)) = Left_p(:,front_back(:,1))+ Left_p(:,front_back(:,2));


R_keep_one = setdiff(V_ind, union(left_right(:,2),bottom_top(:,2)));
 R_keep = setdiff( R_keep_one,front_back(:,2) );

Periodic_Matrix = Left_p(R_keep,R_keep);
end

