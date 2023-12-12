%this function finds K1
function K = optNewForm_minmax(Wb,Wl,Cb,Cl,z,Ah,Ax,bx,bh,L,xe,S,sMin,sMax)
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

I=kron(ones(size(L,2),1),eye(2));
nBarriers = size(Ah,1);
nConvex = size(Ax,1);
L1 = blkdiag(L(:,1),L(:,2),L(:,3),L(:,4));
L2 = reshape(L,[],1);
nLandmark = length(L);

cvx_begin

variables  deltaL K(2,nLandmark*2) deltab lambdaB(nConvex,nBarriers) ...
           lambdaL(nConvex,1) ...
           P1(4,3) P2(4,3)

deltaB = [deltab;deltab;deltab];
minimize(Wb'*deltaB+Wl*deltaL)

subject to
%constraints for safety:
% -(Ah*K*S*L2)' <= (deltaB+Cb*bh)'-bx'*lambdaB 
% % 
[-sMin' sMax']*[P1;P2] <= (deltaB+Cb*bh)'-bx'*lambdaB
[-eye(nLandmark) eye(nLandmark)]*[P1;P2]==(-Ah*K*L1)'
P1>=0
P2>=0

Ax'*lambdaB == (Ah*K*S*I-Ah*Cb)'

%constraints for stability:
bx'*lambdaL<= (deltaL-z'*K*S*L2+Cl*z'*xe)'
Ax'*lambdaL == (-z'*K*S*I+z'*Cl)'


lambdaB<=0
lambdaL<=0
deltab<=0
deltaL<=0
deltab>=-0.1
deltaL>=-0.1 

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