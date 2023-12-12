%this function finds K1 for distance measurments:
% y = L-x and u = Ky+k_added
% this code is for the cell which we want to stop in the center of the cell
% all sides of the cell are barriers and xe is at the center 
function [K,k_added] = optFirstorderWithU_center(Wb,Wl,Cb,Ah,Ax,bx,bh,L,xe)
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
L=reshape(L,[],1);

nBarriers = size(Ah,1);
nConvex = size(Ax,1);
nLandmark = length(L);

cvx_begin

variables  S_l sb K(2,nLandmark) k_added(2,1)... 
           lB(nConvex,nBarriers) lL(nConvex,1)

Sb = kron(ones(nBarriers,1),sb);
minimize(Wb'*Sb+Wl*S_l)

subject to
%constraints for safety:
bx'*lB <= (Sb+Ah*K*L+Cb*bh)'+(Ah*k_added)'
Ax'*lB == (Ah*K*I-Ah*Cb)'

%constraints for stability:
% bx'*lL<= (S_l-z'*K*L+Cl*z'*xe)'-z'*k_added
% Ax'*lL == (-z'*K*I+z'*Cl)'

% we change the above CLF constrain:
abs(K*L - K*kron(ones(4,1),eye(2))*xe+k_added)<=10e-4

% % 
% K*kron(ones(4,1),eye(2)) == [Cb 0;0 Cb]

lB>=0
lL>=0
sb<=0
S_l<=0
sb>=-0.1
S_l>=-0.1

cvx_end
end



