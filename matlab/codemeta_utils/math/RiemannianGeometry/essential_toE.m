function E=essential_toE(Q)
e3hat=[0 -1 0; 1 0 0; 0 0 0];

R1=essential_getR1(Q);
R2=essential_getR2(Q);

NQ=size(Q,3);
E=zeros(3,3,NQ);
for iQ=1:NQ
    E(:,:,iQ)=R1(:,:,iQ)'*e3hat*R2(:,:,iQ);
end
