function vHat=rot3r3_hat(G,v)
N=size(v,2);
vHat=zeros(size(G));
vHat(1:3,1:3,:)=rot_hat(G2R(G),v(1:3,:));
vHat(1:3,4,:)=v(4:6,:);
