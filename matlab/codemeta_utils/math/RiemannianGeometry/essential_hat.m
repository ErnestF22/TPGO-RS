function v=essential_hat(Q,vVec)
v=zeros(6,3,size(vVec,2));
v(1:3,:,:)=multiprod(Q(1:3,:,:),rot_hat(eye(3),vVec(1:3,:)));
v(4:6,:,:)=multiprod(Q(4:6,:,:),rot_hat(eye(3),vVec(4:6,:)));
