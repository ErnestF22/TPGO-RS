function K=optFirstorder_yimg(Wb,Wl,Cb,Cl,Ah,Ax,bx,bh,L,xe,eCamera,d_exit,zMax)
Cl = 0.01;
Cb = 100;
nLandmark = size(L,2);
I=kron(ones(size(L,2),1),eye(2));
Iz = kron(eye(nLandmark),[1 -1]);
nBarriers = size(Ah,1);
nConvex = size(Ax,1);
diagE = kron(eye(nLandmark),eCamera');
L=reshape(L,[],1);
cvx_begin

variables  delta_clf K(2,2*nLandmark) delta_b ... 
           P1(nConvex,nBarriers) P2(nLandmark,nBarriers)...
           P3(nLandmark*2,nBarriers)...
           Q1(nConvex,1) Q2(nLandmark,1) Q3(nLandmark*2,1) ...
           L_cbf(12,24) lambda_clf(12,8)
           

delta_cbf = [delta_b;delta_b;delta_b];
minimize(Wb'*delta_cbf+Wl*delta_clf)
subject to
%CBF:
% [-Ah*K -Cb*Ah]*[yimg;x]<= delta_cbf+Cb*bh
% P1: Ax*x<=bx
% P2: diagE'*yimg == ones(nLandmark,1)
% P3: diag(z)*yimg == L-1*x
%CLF:
% [d_exit'*K Cl*d_exit']*[yimg;x]<= delta_cbf+Cl*d_exit'*xe
% Q1: Ax*x<=bx
% Q2: diagE'*yimg == ones(nLandmark,1)
% Q3: diag(z)*yimg == L-1*x
%form the problem as :
% max:       C*[yimg;x]
%subject to: A*[yimg;x]<=b, x>=0
%%
C_cbf1 = (-Ah*K)';
C_cbf2 = -Cb*Ah';
C_clf1 = (d_exit'*K)';
C_clf2 = Cl*d_exit;
% A1 = [zeros(8,4) diagE' ones(8,8)];
A2 = [Ax' zeros(2,nLandmark) I'];
b = [bx;ones(nLandmark,1);L];

%% CBF
b'*[P1;P2;P3]<=delta_cbf'+Cb*bh'
A2*[P1;P2;P3] >= C_cbf2
% A1*[P1;P2;P3] >= C_cbf1, for all z:
Az = [eye(8);Iz];
bz = [zMax;zeros(4,1)];
RHS = diagE'*P2-C_cbf1;
for i=1:3
    bz'*L_cbf(:,(i-1)*8+1:i*8)<=RHS(:,i)'
    Az'*L_cbf(:,(i-1)*8+1:i*8)>=-diag(P3(:,i))
end
%% CLF
b'*[Q1;Q2;Q3]<=delta_clf+Cl*d_exit'*xe
A2*[Q1;Q2;Q3] >= C_clf2
RHS = diagE'*Q2-C_clf1;
bz'*lambda_clf<=RHS'
Az'*lambda_clf>=-diag(Q3)

%%
P1>=0
Q1>=0
L_cbf(1:8,:)>=0
lambda_clf(1:8,:)>=0
delta_b<=0
delta_b>=-0.1
delta_clf<=0
delta_clf>=-0.1
cvx_end
end