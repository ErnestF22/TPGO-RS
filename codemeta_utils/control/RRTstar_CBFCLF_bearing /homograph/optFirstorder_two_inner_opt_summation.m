%this function finds K1
function K = optFirstorder_two_inner_opt_summation(Wb,Wl,Cb,Cl,exitDir,Ah,Ax,bx,bh,L,xe,Zmax,eCamera)
Zmax = Zmax'*blkdiag([1 1 0 0],[1 1 0 0],[1 1 0 0],[1 1 0 0]);
L=reshape(L,[],1);
Ie = blkdiag(eCamera',eCamera',eCamera',eCamera');
Ix = repmat(eye(2),4,1);
digagL = blkdiag([1 1],[1 1],[1 1],[1 1]);
Ap = [1 0 1 -1;0 1 -1 1];
Ap = blkdiag(Ap,Ap,Ap,Ap);
nBarriers = size(Ah,1);
nConvex = size(Ax,1);
nLandmarks = length(L)/2;

A1 = [Ax zeros(nConvex,nLandmarks*2) zeros(nConvex,nLandmarks*4)];
A2 = [zeros(nLandmarks,2) Ie zeros(nLandmarks,nLandmarks*4)];
A3 = [digagL*repmat(eye(2),nLandmarks,1) zeros(nLandmarks,nLandmarks*2) Zmax];
A4 = [zeros(nLandmarks*2,2) -eye(nLandmarks*2) Ap];

A = [A1;A2;A3;A4];
A = A';
A1 = A(1:2,:);
A2 = A(3:2+nLandmarks*2,:);
A3 = A(3+nLandmarks*2:end,:);

b1 = bx;
b2 = ones(nLandmarks,1);
b3 = digagL*L;
b4 = zeros(nLandmarks*2,1);

cvx_begin

variables sb S_L K(2,nLandmarks*2) ...
    Pb1(nConvex,nBarriers) Pb2(nLandmarks,nBarriers)...
    Pb3(nLandmarks,nBarriers) Pb4(nLandmarks*2,nBarriers)  ...
    Pl1(nConvex,1) Pl2(nLandmarks,1) Pl3(nLandmarks,1) Pl4(nLandmarks*2,1)

S_B = [sb;sb;sb];      

CH1 = -Cb*Ah;
CH2 = -Ah*K;
CH3 = zeros(nBarriers,nLandmarks*4);


CL1 = Cl*exitDir'; 
CL2 = exitDir'*K; 
CL3 = zeros(1,nLandmarks*4);

b = [b1;b2;b3;b4];

minimize(Wb'*S_B+Wl*S_L)

% constraints for safety:
[Pb1;Pb2;Pb3;Pb4]'*b <= S_B+Cb*bh
A1*[Pb1;Pb2;Pb3;Pb4]<=CH1'
A2*[Pb1;Pb2;Pb3;Pb4]<=CH2'
A3*[Pb1;Pb2;Pb3;Pb4]<=CH3'
Pb1>=0
Pb4>=0

% constraints for stability:
[Pl1;Pl2;Pl3;Pl4]'*b <= S_L+Cl*exitDir'*xe
A1*[Pl1;Pl2;Pl3;Pl4]<=CL1'
A2*[Pl1;Pl2;Pl3;Pl4]<=CL2'
A3*[Pl1;Pl2;Pl3;Pl4]<=CL3'
Pl1>=0
Pl4>=0


sb<=0
S_L<=0
sb>=-0.1
S_L>=-0.1
% 
% K*kron(ones(4,1),eye(2)) == [Cb 0;0 Cb]
% K*reshape(L,[],1) == K*kron(ones(4,1),eye(2))*xe

cvx_end
end


