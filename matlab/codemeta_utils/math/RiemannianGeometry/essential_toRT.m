function [R,T]=essential_toRT(Q)
R1=essential_getR1(Q);
R2=essential_getR2(Q);
N=size(Q,3);
R=zeros(3,3,N);
T=zeros(3,N);
for iR=1:N
    R(:,:,iR)=R1(:,:,iR)'*R2(:,:,iR);
    T(:,iR)=R1(:,:,iR)'*[0;0;1];
end
