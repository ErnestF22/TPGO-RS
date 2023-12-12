%this function finds K1
function [K] = optFirstorder(Wb,Wl,Cb,Cl,exitDir,Ah,Ax,bx,bh,xe,L)
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
%d:dimestion
%y:bering measurments  EI = diag(reshape(E,[],1));
%L:landmark matrix
nBarriers = size(Ah,1);
nConvex = size(Ax,1);
nLandmarks = size(L,2);

% we need to define A and b such that AX=b where X=[x;y;z]
vecL = reshape(L,[],1);
ay = eye(2*nLandmarks);
al = blkdiag([1;1],[1;1],[1;1],[1;1]);
az = -diag(vecL)*al;
ax = repmat(eye(2),4,1);
%ax ay aL;
%zeros(2*nLandmarks,1);
A = [Ax zeros(4,2*nLandmarks) zeros(4,nLandmarks);
    ax ay az;
    zeros(nLandmarks,2) zeros(nLandmarks,2*nLandmarks) eye(nLandmarks)
    ];

b = [bx;zeros(2*nLandmarks,1);ones(nLandmarks,1)];

A = A';
A1 = A(1:2,:);
A2 = A(3:10,:);
A3 = A(11:14,:);

CH1 = -Cb*Ah;
CH3 = zeros(nBarriers,nLandmarks);
CL1 = Cl*exitDir';
CL3 = zeros(1,nLandmarks);



cvx_begin

    variables Sb(nBarriers,1) S_l K(2,2*nLandmarks) ...
             Pb1(4,3) Pb2(2*nLandmarks,3) Pb3(nLandmarks,3)...
             Pl1(4,1) Pl2(2*nLandmarks,1) Pl3(nLandmarks,1)
CH2 = -Ah*K;
CL2 = exitDir'*K;


minimize(Wb'*Sb+Wl*S_l)

%constraints for safety:
[Pb1;Pb2;Pb3]'*b <= Sb+Cb*bh
A1*[Pb1;Pb2;Pb3]==CH1'
A2*[Pb1;Pb2;Pb3]==CH2'
A3*[Pb1;Pb2;Pb3]>=CH3'
Pb1>=0
Pb2>=0


%constraints for stability:
[Pl1;Pl2;Pl3]'*b <= S_l+Cl*exitDir'*xe
A1*[Pl1;Pl2;Pl3]==CL1'
A2*[Pl1;Pl2;Pl3]==CL2'
A3*[Pl1;Pl2;Pl3]>=CL3'
Pl1>=0
Pl2>=0




Sb<=0
S_l<=0

cvx_end
end