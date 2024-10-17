function epipolarConstraint_test
x1=randn(2,1);
x2=randn(2,1);
[Rt,vt,R0,v0,vVec]=rot_randGeodFun(eye(3));
[Tt,dTt,T0,dT]=real_randGeodFun(randn(3,1));

figure(1)
check_der(@(t) evaluateFandDF(x1,x2,Rt,Tt,vVec,dT,t),'function')

wFix=randn(6,1);

figure(2)
check_der(@(t) evaluateDFandDdF(x1,x2,Rt,Tt,vVec,dT,wFix,t),'function')

function [f,df]=evaluateFandDF(x1,x2,Rt,Tt,vVec,dT,t)
[f,Jf]=epipolarConstraint(Rt(t),Tt(t),x1,x2);
df=Jf*[vVec;dT];

% w=[vVec;dT];
% xL=[x1;1];
% xR=[x2;1];
% R=Rt(t);
% T=Tt(t);
% E=hat(T)*R;
% f=xL'*E*xR;
% df=-xL'*hat(T)*R*hat(xR)*w(1:3)+xR'*R'*hat(xL)*w(4:6);


function [df,ddf]=evaluateDFandDdF(x1,x2,Rt,Tt,vVec,dT,wFix,t)
[~,Jf,Hf]=epipolarConstraint(Rt(t),Tt(t),x1,x2);
w=[vVec;dT];
df=Jf*wFix;
ddf=w'*Hf*wFix;
% xL=[x1;1];
% xR=[x2;1];
% R=Rt(t);
% T=Tt(t);
% E=hat(T)*R;
% H1=-hat(R'*hat(T)*xL)*hat(xR);
% H2=hat(xR)*R'*hat(xL);
% % ddf=w'*[H1 H2; H2' zeros(3)]*wFix;
% %ddf=w(4:6)'*H2'*wFix(1:3)+w(1:3)'*H1*wFix(1:3)+w(1:3)'*H2*wFix(4:6);
