%this function finds K1Max
function [K,k_added] = optFirstorderWithS_test_inequality_trace(Wb,Wl,Cb,Cl,z,Ah,Ax,bx,bh,L,xe,sMin,sMax)

ones_diag = kron(eye(size(L,2)),[1;1]');
I=kron(ones(size(L,2),1),eye(2));
L=reshape(L,[],1);
nBarriers = size(Ah,1);
nConvex = size(Ax,1);
nLandmark = length(L);


cvx_begin

variables  delta_clf K(2,nLandmark) k_added(2,1) delta_b ...
    lambda_cbf(nConvex,nBarriers) lambda_clf(nConvex,1)...
    P1(nLandmark,nBarriers) P2(nLandmark,nBarriers) P3(nLandmark,nBarriers)...
    Q1(nLandmark,1) Q2(nLandmark,1) Q3(nLandmark,1)


dual variables y1 y2 y3 y4 y5 y6

delta_cbf = [delta_b;delta_b;delta_b];
minimize(Wb'*delta_cbf+Wl*delta_clf)

subject to
%%
%constraints for safety:
% bx'*lambda_cbf <= (delta_cbf+Ah*K*L+Cb*bh)'
% Ax'*lambda_cbf >= (Ah*K*I-Ah*Cb)'
%%
% for i=1:nBarriers
%     -S'*((diag(L*Ah(i,:)*K))'*ones_diag')'<=-bx'*lambda_cbf(:,i)+delta_cbf(i)+Cb*bh(i)'
%     [eye(nLandmark/2);-eye(nLandmark/2)]*S<=[sMax;-sMin]
%     RHS = Ax'*lambda_cbf(:,i)+Cb*Ah(i,:)';
%     S'*((diag(I(:,1)*Ah(i,:)*K))'*ones_diag')'<= RHS(1)
%     S'*((diag(I(:,2)*Ah(i,:)*K))'*ones_diag')'<= RHS(2)
% end
%%
%%
for i=1:nBarriers
    [sMax;-sMin]'*P1(:,i)<=-bx'*lambda_cbf(:,i)+delta_cbf(i)+Cb*bh(i)'+Ah(i,:)*k_added
    [eye(nLandmark/2);-eye(nLandmark/2)]'*P1(:,i)>= -((diag(L*Ah(i,:)*K))'*ones_diag')'
    
    RHS = Ax'*lambda_cbf(:,i)+Cb*Ah(i,:)';
    
    [sMax;-sMin]'*P2(:,i)<= RHS(1)
    [eye(nLandmark/2);-eye(nLandmark/2)]'*P2(:,i)>=((diag(I(:,1)*Ah(i,:)*K))'*ones_diag')'
    
    [sMax;-sMin]'*P3(:,i)<= RHS(2)
    [eye(nLandmark/2);-eye(nLandmark/2)]'*P3(:,i)>=((diag(I(:,2)*Ah(i,:)*K))'*ones_diag')'

end
%% Test to understand K
% K(2,1)==0;
% K(2,3)==0;
% K(1,2)==0;
% K(1,4)==0;
%%
%constraints for stability:
% bx'*lambda_clf<= (delta_clf-z'*K*L+Cl*z'*xe)'
% Ax'*lambda_clf >= (-z'*K*I+z'*Cl)'
%%
% S'*((diag(L*z'*K))'*ones_diag')'<=-bx'*lambda_clf+delta_clf+(Cl*z'*xe)'
% RHS = Ax'*lambda_clf-Cl*z;
% -S'*((diag(I(:,1)*z'*K))'*ones_diag')'<= RHS(1)
% -S'*((diag(I(:,2)*z'*K))'*ones_diag')'<= RHS(2)
%%

y1: [sMax;-sMin]'*Q1 <= -bx'*lambda_clf+delta_clf+(Cl*z'*xe)'-z'*k_added
y2: [eye(nLandmark/2);-eye(nLandmark/2)]'*Q1 >= ((diag(L*z'*K))'*ones_diag')'

RHS = Ax'*lambda_clf-Cl*z;
%
y3: [sMax;-sMin]'*Q2<= RHS(1)
y4: [eye(nLandmark/2);-eye(nLandmark/2)]'*Q2 >= -((diag(I(:,1)*z'*K))'*ones_diag')'

y5: [sMax;-sMin]'*Q3<= RHS(2)
y6: [eye(nLandmark/2);-eye(nLandmark/2)]'*Q3 >= -((diag(I(:,2)*z'*K))'*ones_diag')'
%%

lambda_clf>=0
lambda_cbf>=0

delta_b<=0
delta_clf<=0
delta_b>=-0.1
delta_clf>=-0.1

P1>=0
P2>=0
P3>=0
Q1>=0
Q2>=0
Q3>=0

cvx_end
end

