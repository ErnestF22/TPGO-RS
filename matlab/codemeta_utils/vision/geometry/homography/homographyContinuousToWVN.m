function [w,v,nVec]=homographyContinuousToWVN(H)
[U,Lambda]=eig(H+H');
l=flipud(diag(Lambda));
U=fliplr(U);

vn1=(sqrt(2*l(1))*U(:,1)+sqrt(-2*l(3))*U(:,3))/2;
vn2=(sqrt(2*l(1))*U(:,1)-sqrt(-2*l(3))*U(:,3))/2;

v=zeros(3,2);
nVec=zeros(3,2);

if vn2(3)<0
    v(:,1)=vn1;
    nVec(:,1)=vn2;
else
    v(:,1)=-vn1;
    nVec(:,1)=-vn2;
end    

if vn1(3)<0
    v(:,2)=vn2;
    nVec(:,2)=vn1;
else
    v(:,2)=-vn2;
    nVec(:,2)=-vn1;
end

w=zeros(3,2);
for iw=1:2
    w(:,iw)=-vee3(H-v(:,iw)*nVec(:,iw)');
end
