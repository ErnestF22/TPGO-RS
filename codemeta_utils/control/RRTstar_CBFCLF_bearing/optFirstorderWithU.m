%this function finds K1
function [K1] = optFirstorderWithU(Wb,Wl,Cb,Cl,z,Ah,Au,Ax,bx,bh,bu,y,xe,b)
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
Y=zeros(2*size(y,2),1);
for i=1:size(y,2)
    Y(2*i-1)=y(1,i);
    Y(2*i)=y(2,i);
end
nBarriers = size(Ah,1);
nConvex = size(Ax,1);
nVelocity = size(Au,1); 
nLandmark=length(Y);
I=kron(ones(size(y,2),1),eye(2));
cvx_begin

    variables Sb(nBarriers,1) S_l Su k1 k2 k3 k4 k5 k6 ...%K1(2,nLandmark) ...
              lB(nConvex,nBarriers) lL(nConvex,1) lU(nConvex,nVelocity)
    K1 = [k1 0 k2 0 k3 0;0 k4 0 k5 0 k6];
%     Au = -K1*I;
%     bu = Bu - K1*Y;
    minimize(Wb'*Sb+Wl*S_l)

    subject to
    
    %constraints for safety:
    lB'*bx <= Sb+Cb*bh+Ah*K1*Y
    Ax'*lB == (Ah*K1*I-Ah*Cb)'
    
    %constraints for stability:
    lL'*bx<= S_l-z'*K1*Y+Cl*z'*xe+Cl*b
    Ax'*lL == (-z'*K1*I+z'*Cl)'
    
%     %constraints for actuator:
%     lU'*bx<= Su+bu-Au*K1*Y+5
%     Ax'*lU == (-Au*K1*I)'
    
    

        lB>=0
        lL>=0
%         lU>=0
        Sb<=0
        S_l<=0
%         Su<=0
        
cvx_end
end