%this function finds K1
function [K] = optFirstorderWithS(Wb,Wl,Cb,Cl,z,Ah,Ax,bx,bh,L,xe,sMin,sMax)
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
%d:dimestion 0.068* 0.1*

nLandmark = length(L);
I=kron(ones(size(L,2),1),eye(2));
L1 = blkdiag(L(:,1),L(:,2),L(:,3),L(:,4));
L2 = reshape(L,[],1);
nBarriers = size(Ah,1);
nConvex = size(Ax,1);

cvx_begin

variables  S_l K(2,nLandmark*2) sb ...
           lB(nConvex,nBarriers) lL(nConvex,1)...
           P1(4,3) P2(4,3) P3(4,3) P4(4,3) P5(4,3) P6(4,3)...
           P7(4,1) P8(4,1) P9(4,1) P10(4,1) P11(4,1) P12(4,1)
        

Sb = [sb;sb;sb];
minimize(Wb'*Sb+Wl*S_l)

subject to
% constraints for safety:
[-sMin' sMax']*[P1;P2] <= (Sb+Cb*bh)'-bx'*lB
[-eye(nLandmark) eye(nLandmark)]*[P1;P2]>=(-Ah*K*L1)'

for i=1:nBarriers
    ah = Ah(i,:)';
    [-sMin' sMax']*[P3(:,i);P4(:,i)]==[1 0]*(-ah*Cb-Ax'*lB(:,i))
    [-eye(nLandmark) eye(nLandmark)]*[P3(:,i);P4(:,i)]>=-([1 0]*K*kron(eye(nLandmark),ah))'
end

for i=1:nBarriers
    ah = Ah(i,:)';
    [-sMin' sMax']*[P5(:,i);P6(:,i)]==[0 1]*(-ah*Cb-Ax'*lB(:,i))
    [-eye(nLandmark) eye(nLandmark)]*[P5(:,i);P6(:,i)]>=-([0 1]*K*kron(eye(nLandmark),ah))'
end

%constraints for stability:
[-sMin' sMax']*[P7;P8] <= (S_l+Cl*z'*xe)-bx'*lL
[-eye(nLandmark) eye(nLandmark)]*[P7;P8]>= (z'*K*L1)'

[-sMin' sMax']*[P9;P10]==[1 0]*(z*Cl-Ax'*lL)
[-eye(nLandmark) eye(nLandmark)]*[P9;P10]>=([1 0]*K*kron(eye(nLandmark),z))'

[-sMin' sMax']*[P11;P12]==[0 1]*(z*Cl-Ax'*lL)
[-eye(nLandmark) eye(nLandmark)]*[P11;P12]>=([0 1]*K*kron(eye(nLandmark),z))'
% %

lB>=0
lL>=0
sb<=0
S_l<=0
sb>=-0.1
S_l>=-0.1

P1>=0
P2>=0
P3>=0
P4>=0
P5>=0
P6>=0
P7>=0
P8>=0
P9>=0
P10>=0
P11>=0
P12>=0


cvx_end
end

%%
% bx'*lB <= (Sb+Ah*KS*L+Cb*bh)'
% Ax'*lB == (Ah*KS*I-Ah*Cb)'
% 
% %constraints for stability:
% bx'*lL<= (S_l-z'*KS*L+Cl*z'*xe)'
% Ax'*lL == (-z'*KS*I+z'*Cl)'
% 
% (-ah*K*L1)'*s<=(sb+Cb*bh(i))'-bx'*lB(:,i)
% -([1 0]*K*kron(eye(nLandmark),ah))'*s==[1 0]*(-ah*Cb-Ax'*lB(:,i))
% -([0 1]*K*kron(eye(nLandmark),ah))'s==[0 1]*(-ah*Cb-Ax'*lB(:,i))
% 
% (z'*K*L1)'*s <= (S_l+Cl*z'*xe)-bx'*lL
% ([1 0]*K*kron(eye(nLandmark),z))'*s==[1 0]*(z*Cl-Ax'*lL)
% ([0 1]*K*kron(eye(nLandmark),z))'*s==[0 1]*(z*Cl-Ax'*lL)