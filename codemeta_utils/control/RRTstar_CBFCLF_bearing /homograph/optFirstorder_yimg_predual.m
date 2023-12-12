function [K,k_added]=optFirstorder_yimg_predual(Wb,Wl,Cb,Cl,Ah,Ax,bx,bh,L,xe,eCamera,d_exit,zMax,zMin)
% Cl = 0.001;
% Cb = 500;
nLandmark = size(L,2);
I=kron(ones(size(L,2),1),eye(2));
nBarriers = size(Ah,1);
nConvex = size(Ax,1);
diagE = kron(eye(nLandmark),eCamera');
L=reshape(L,[],1);
cvx_begin

variables  delta_clf K(2,2*nLandmark) delta_b ...
    P1(nConvex,nBarriers) P2(nLandmark,nBarriers)...
    P3(nLandmark*2,nBarriers)...
    Q1(nConvex,1) Q2(nLandmark,1) Q3(nLandmark*2,1)...
    k_added(2,1) ...
    lambda_cbf(2*nLandmark*2*nLandmark,nBarriers) lambda_clf(2*nLandmark*2*nLandmark,1)


delta_cbf = [delta_b;delta_b;delta_b];
minimize(Wb'*delta_cbf+Wl*delta_clf)
subject to
%CBF:
% [-Ah*K -Cb*Ah]*[yimg;x]<= delta_cbf+Cb*bh
% P1: Ax*x<=bx
% P2: diagE'*yimg == ones(nLandmark,1)
% P3: diag(z)*yimg == L-1*x
%CLF:
% [d_exit'*K Cl*d_exit']*[yimg;x]<= delta_clf+Cl*d_exit'*xe
% Q1: Ax*x<=bx
% Q2: diagE'*yimg == ones(nLandmark,1)
% Q3: diag(z)*yimg == L-1*x
%form the problem as :
% max:       C*[yimg;x]
%subject to: A*[yimg;x]<=b, x>=0
%%
C_cbf_y = (-Ah*K)';
C_cbf_x = -Cb*Ah';
C_clf_y = (d_exit'*K)';
C_clf_x = Cl*d_exit;
% A1 = [zeros(nLandmark*2,nConvex) diagE' eye(nLandmark*2)]; %y_img free
A_x = [Ax' zeros(2,nLandmark) I'];%x>=0
b = [bx;ones(nLandmark,1);L];

%% CBF
b'*[P1;P2;P3]<=delta_cbf'+Cb*bh'+(Ah*k_added)'
% A1*[P1;P2;P3] == C_cbf1
A_x*[P1;P2;P3] >= C_cbf_x %x>=0 if free then ==

RHS = -diagE'*P2+C_cbf_y;
Az = [eye(nLandmark);-eye(nLandmark)];
bz = [zMax;-zMin];
Iz = kron(eye(nLandmark),[1 1]);
for i=1:nBarriers
    rhs = (Iz*diag(P3(:,i)))';
    for t=1:2*nLandmark
        bz'*lambda_cbf((t-1)*(2*nLandmark)+1:t*(2*nLandmark),i)==RHS(t,i)'
        Az'*lambda_cbf((t-1)*(2*nLandmark)+1:t*(2*nLandmark),i)>= rhs(t)
    end
end
%% CLF
b'*[Q1;Q2;Q3]<=delta_clf+Cl*d_exit'*xe-(d_exit'*k_added)'
% A1*[Q1;Q2;Q3] == C_clf1
A_x*[Q1;Q2;Q3] >= C_clf_x %x>=0 if free then ==

RHS = -diagE'*Q2+C_clf_y;
rhs = (Iz*diag(Q3))';
for t=1:2*nLandmark
    bz'*lambda_cbf((t-1)*(2*nLandmark)+1:t*(2*nLandmark),i)==RHS(t,1)'
    Az'*lambda_cbf((t-1)*(2*nLandmark)+1:t*(2*nLandmark),i)>= rhs(t) 
    %if ==: the PV is free, if >=: PV>=0, if <=: PV<=0
end

%%
P1>=0
Q1>=0
delta_b<=0
delta_b>=-0.1
delta_clf<=0
delta_clf>=-0.1
lambda_cbf>=0
lambda_clf>=0
cvx_end
end