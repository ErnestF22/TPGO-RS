%this function finds K1 for distance measurments:
% y = L-x and u = Ky+k_added
function K = optFirstorderU_newCost(Cb,Cl,y,M,xc1,xc2)
%Wb:steering factor for cbf
%Wl:steering factor for clf
%Cb:cbf coefficoient
%Cl:clf coefficoient
%z:exit direction
%Ah:barrier matrix
%bh:barrier bias
%Au:velocity matrix
%bu:velocity bias
%Ax:convex cell matrix
%bx:convex cell bias
%Y: vectorized of landmarks
%d:dimestion

numCells = 2;
L = y;
I =kron(ones(size(L,2),1),eye(2));
L=reshape(L,[],1);
nLandmark = length(L);

Ax1 = M(1).Ax;
bx1 = M(1).bx;
Ah1 = M(1).Ah;
bh1 = M(1).bh;
z1 = M(1).z;
xe1 = M(1).xe;
Ax2 = M(2).Ax;
bx2 = M(2).bx;
Ah2 = M(2).Ah;
bh2 = M(2).bh;
z2 = M(2).z;
xe2 = M(2).xe;

nBarriers = size(Ah1,1);
nConvex = size(Ax1,1);


cvx_begin
variables  S1_l sb1 K1(2,nLandmark) k1_added(2,1)... 
           S2_l sb2 K2(2,nLandmark) k2_added(2,1)...
           lB1(nConvex,nBarriers) lL1(nConvex,1) ...
           lB2(nConvex,nBarriers) lL2(nConvex,1)
       
Sb1 = kron(ones(nBarriers,1),sb1);
Sb2 = kron(ones(nBarriers,1),sb2);
c1 = K1*(reshape((y-xc1*ones(1,size(nLandmark,2))),[],1))+k1_added;
c2 = K2*(reshape((y-xc1*ones(1,size(nLandmark,2))),[],1))+k2_added;
c3 = K1*(reshape((y-xc2*ones(1,size(nLandmark,2))),[],1))+k1_added;
c4 = K2*(reshape((y-xc2*ones(1,size(nLandmark,2))),[],1))+k2_added;

minimize (norm((c1-c2),2)+norm((c3-c4),2))

subject to
%constraints for safety:
bx1'*lB1 <= (Sb1+Ah1*K1*L+Cb*bh1)'+(Ah1*k1_added)'
Ax1'*lB1 == (Ah1*K1*I-Ah1*Cb)'

%constraints for stability:
bx1'*lL1<= (S1_l-z1'*K1*L+Cl*z1'*xe1)'-z1'*k1_added
Ax1'*lL1 == (-z1'*K1*I+z1'*Cl)'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%constraints for safety:
bx2'*lB2 <= (Sb2+Ah2*K2*L+Cb*bh2)'+(Ah2*k2_added)'
Ax2'*lB2 == (Ah2*K2*I-Ah2*Cb)'

%constraints for stability:
bx2'*lL2<= (S2_l-z2'*K2*L+Cl*z2'*xe2)'-z2'*k2_added
Ax2'*lL2 == (-z2'*K2*I+z2'*Cl)'

lB1>=0
lL1>=0
sb1<=0
S1_l<=0
sb1>=-0.1
S1_l>=-0.1

lB2>=0
lL2>=0
sb2<=0
S2_l<=0
sb2>=-0.1
S2_l>=-0.1

k1_added<=2
k2_added<=2

cvx_end

K(1).K = K1;
K(2).K = K2;
K(1).added = k1_added;
K(2).added = k2_added;

end



