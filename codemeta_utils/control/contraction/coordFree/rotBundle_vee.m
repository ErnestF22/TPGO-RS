function x=rotBundle_vee(R,X)
x=zeros(6,1);
x(1:3)=rot_vee(R,X(1:3,:));
x(4:6)=rot_vee(R,X(4:6,:));
