function x=rot2euler(R)
x(1,:)=atan(R(3,2,:)./R(3,3,:));
x(2,:)=-asin(R(3,1,:));
x(3,:)=atan(R(2,1,:)./R(1,1,:));
