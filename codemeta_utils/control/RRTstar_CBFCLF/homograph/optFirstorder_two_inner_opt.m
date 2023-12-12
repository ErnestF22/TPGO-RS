%this function finds K1
%his does not work as z1 and z2 are the same. 
function K = optFirstorder_two_inner_opt(Wb,Wl,Cb,Cl,exitDir,Ah,Ax,bx,bh,L,xe,Zmax,eCamera)
Zmax = diag(Zmax); %Zmax*Pb

L=reshape(L,[],1);
Ie = blkdiag(eCamera',eCamera',eCamera',eCamera');
diagE = diag([eCamera;eCamera;eCamera;eCamera]);
nBarriers = size(Ah,1);
nConvex = size(Ax,1);
nLandmarks = length(L)/2;
Ix = kron(ones(size(L,1)/2,1),eye(2));

A1 = [Ax zeros(nConvex,nLandmarks*2) zeros(nConvex,nLandmarks*2)]; %AxX<=bx
A2 = [diagE*Ix zeros(2*nLandmarks,nLandmarks*2) Zmax];%zP=e*(L-x)
A3 = [zeros(nLandmarks*2,2) -eye(nLandmarks*2) eye(nLandmarks*2)];%P-yimg<=0
A4 = [zeros(nLandmarks,2) Ie zeros(nLandmarks,nLandmarks*2)];%e*yimg=1

A = [A1;A2;A3;A4];
A = A';
A1 = A(1:2,:); %x
A2 = A(3:2+nLandmarks*2,:); %yimg
A3 = A(3+nLandmarks*2:end,:); %P

b1 = bx;
b2 = diagE*L;
b3 = zeros(nLandmarks*2,1);
b4 = ones(nLandmarks,1);
b = [b1;b2;b3;b4];

cvx_begin
variables sb S_L K(2,nLandmarks*2) Pb(2*nLandmarks,1)...
    lambda_1(nConvex,nBarriers) lambda_2(2*nLandmarks,nBarriers)...
    lambda_3(2*nLandmarks,nBarriers) lambda_4(nLandmarks,nBarriers)  ...
    alpha_1(nConvex,1) alpha_2(2*nLandmarks,1) alpha_3(2*nLandmarks,1)...
    alpha_4(nLandmarks,1)

lambda = [lambda_1;lambda_2;lambda_3;lambda_4];
alpha = [alpha_1;alpha_2;alpha_3;alpha_4];
S_B = [sb;sb;sb]; 

CH1 = -Cb*Ah;
CH2 = -Ah*K;
CH3 = zeros(nBarriers,2*nLandmarks);


CL1 = Cl*exitDir'; 
CL2 = exitDir'*K; 
CL3 = zeros(1,2*nLandmarks);

% 
minimize(Wb'*S_B+Wl*S_L)

%constraints for safety:
lambda'*b <= S_B+Cb*bh
A1*lambda==CH1'
A2*lambda==CH2'
A3*lambda<=CH3'
lambda_1>=0
lambda_4>=0

%constraints for stability:
alpha'*b <= S_L+Cl*exitDir'*xe
A1*alpha ==CL1'
A2*alpha ==CL2'
A3*alpha <=CL3'
alpha_1>=0
alpha_4>=0

% K>=-3
% K<=3
sb<=0
S_L<=0
sb>=-0.1
S_L>=-0.1

% sum(K(1,2:2:end))==0;
% sum(K(2,1:2:end))==0;
% sum(K(1,1:2:end))==20;
% sum(K(2,2:2:end))==20;


cvx_end

 end



