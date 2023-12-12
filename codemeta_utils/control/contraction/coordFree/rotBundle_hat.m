function X=rotBundle_hat(R,x)
X=zeros(6,3);
X(1:3,:)=rot_hat(R,x(1:3));
X(4:6,:)=rot_hat(R,x(4:6));
