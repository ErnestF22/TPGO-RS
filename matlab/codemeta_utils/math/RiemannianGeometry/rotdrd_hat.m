function vHat=rotdrd_hat(G,v)
d=size(G,2)-1;
dR=rot_dim(eye(d));
vHat=zeros(size(G));
vHat(1:d,1:d,:)=rot_hat(G2R(G),v(1:dR,:));
vHat(1:d,d+1,:)=v(dR+1:end,:);
