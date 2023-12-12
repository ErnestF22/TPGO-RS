%this function finds K1
function K = optFirstorderwithBearing_twoInnerOpt(Wb,Wl,Cb,Cl,exitDir,Ah,Ax,bx,bh,xe,L,Zmax,E)
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
%y:bering measurments  
%EI = diag(reshape(E,[],1));
%L:landmark matrix

nBarriers = size(Ah,1);
nConvex = size(Ax,1);
nLandmarks = size(L,2);

Ix = repmat(eye(2),nLandmarks,1);

ETy = blkdiag(E',E',E',E');
ETL = blkdiag(E',E',E',E')*reshape(L,[],1);

%A=[Ax<=b, eTy=1, zmaxP=L-x, P-y<=0];

A = [Ax zeros(nConvex,2*nLandmarks) zeros(nConvex,nLandmarks);
    Ix zeros(nLandmarks,2*nLandmarks) Zmax*eye(nLandmarks);
    zeros(nLandmarks,2)  -eye(2*nLandmarks) eye(nLandmarks)
    ];

b = [bx;ETL;zeros(nLandmarks,1)];

A = A';
A1 = A(1:2,:);
A2 = A(3:2+2*nLandmarks,:);
A3 = A(3+2*nLandmarks:end,:);

%
cvx_begin

    variables Sb(nBarriers,1) S_l K(2,2*nLandmarks) ...
             Pb1(nConvex,nBarriers) Pb2(nLandmarks,nBarriers)...
             Pb3(nLandmarks,nBarriers) Pb4(nLandmarks,nBarriers)  ...
             Pl1(nConvex,1) Pl2(nLandmarks,1) Pl3(nLandmarks,1)...
             Pl4(nLandmarks,1)
             
         
CH1 = -Cb*Ah;
CH2 = -Ah*K;
CH3 = zeros(nBarriers,nLandmarks);
CL1 = Cl*exitDir';
CL2 = exitDir'*K;
CL3 = zeros(1,nLandmarks);


minimize(1)%Wb'*Sb+Wl*S_l

%constraints for safety:
[Pb1;Pb2;Pb3;Pb4]'*b <= -1+Cb*bh
A1*[Pb1;Pb2;Pb3;Pb4]==CH1'
A2*[Pb1;Pb2;Pb3;Pb4]==CH2'
A3*[Pb1;Pb2;Pb3;Pb4]>=CH3'
Pb1>=0
Pb4>=0



%constraints for stability:
[Pl1;Pl2;Pl3;Pl4]'*b <= -1+Cl*exitDir'*xe
A1*[Pl1;Pl2;Pl3;Pl4]==CL1'
A2*[Pl1;Pl2;Pl3;Pl4]==CL2'
A3*[Pl1;Pl2;Pl3;Pl4]>=CL3'
Pl1>=0
Pl4>=0

K*ones(8,1)==[1;1]
K*reshape(L,[],1) == Cb*xe

Sb<0
S_l<0

cvx_end
end