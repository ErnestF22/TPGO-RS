%this function finds K1
% this function does not work as it multiples two variable which makes it
% nonlinear! 
function K = optFirstorder_two_inner_opt_withoutPb(Wb,Wl,Cb,Cl,exitDir,Ah,Ax,bx,bh,L,xe,Zmax,eCamera)
Zmax = Zmax'*blkdiag([1 1],[1 1],[1 1],[1 1]); %Zmax*Pb
L=reshape(L,[],1);
Ie = blkdiag(eCamera',eCamera',eCamera',eCamera');
Ix = repmat(eye(2),4,1);
digagL = blkdiag([1 1],[1 1],[1 1],[1 1]);
Ap = [1 0 1 -1;0 1 -1 1];
Ap = blkdiag(Ap,Ap,Ap,Ap); %AP*Pb
nBarriers = size(Ah,1);
nConvex = size(Ax,1);
nLandmarks = length(L)/2;

A1 = [Ax zeros(nConvex,nLandmarks*2)]; %AxX<=bx
A2 = [digagL*repmat(eye(2),nLandmarks,1) zeros(nLandmarks,nLandmarks*2)];%X=L-zP
A3 = [zeros(nLandmarks*2,2) -eye(nLandmarks*2)];%-yimg<=-P
A4 = [zeros(nLandmarks,2) Ie];%e*yimg=1

A = [A1;A2;A3;A4];
A = A';

cvx_begin

variables sb S_L K(2,nLandmarks*2) Pb(2*nLandmarks,1)...
    lambda_1(nConvex,nBarriers) lambda_2(nLandmarks,nBarriers)...
    lambda_3(nLandmarks,nBarriers) lambda_4(nLandmarks*2,nBarriers)  ...
    alpha_1(nConvex,1) alpha_2(nLandmarks,1) alpha_3(nLandmarks,1)...
    alpha_4(nLandmarks*2,1)
lambda = [lambda_1;lambda_2;lambda_3;lambda_4];
alpha = [alpha_1;alpha_2;alpha_3;alpha_4];
S_B = [sb;sb;sb];      
b1 = bx;
b2 = digagL*L-Zmax*Pb;
b3 = -Pb;
b4 = ones(nLandmarks,1);
b = [b1;b2;b3;b4];

CH1 = -Cb*Ah;
CH2 = -Ah*K;

CL1 = Cl*exitDir'; 
CL2 = exitDir'*K; 

minimize(Wb'*S_B+Wl*S_L)

% constraints for safety:
lambda'*b <= S_B+Cb*bh
A*lambda<=[CH1;CH2]'
Pb1>=0
Pb3>=0

% constraints for stability:
alpha'*b <= S_L+Cl*exitDir'*xe
A*alpha<=[CL1;CL2]'
Pl1>=0
Pl3>=0


sb<=0
S_L<=0
sb>=-0.1
S_L>=-0.1
Pb<=0

% K*kron(ones(4,1),eye(2)) == [Cb 0;0 Cb]
% K*reshape(L,[],1) == K*kron(ones(4,1),eye(2))*xe

cvx_end

 end



