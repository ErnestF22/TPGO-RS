function poseEstimationCost_test
resetRands(3)
X=randn(3,1);
x=randn(2,1);
[Rt,vt,R0,v0,vVec]=rot_randGeodFun(eye(3));
[Tt,dTt,T0,dT]=real_randGeodFun(randn(3,1));
wFix=rand(6,1);

% figure(1)
% check_der(@(t) evaluateFandDF(X,x,Rt,Tt,vVec,dT,t),'function')
% 
% figure(2)
% check_der(@(t) evaluateDFandDdF(X,x,Rt,Tt,vVec,dT,wFix,t),'function')

[xt,dxt,x0,dx]=real_randGeodFun(x);

figure(3)
check_der(@(t) evaluateGradFanddGradF(X,xt,R0,T0,dx,t),'function')

function [f,df]=evaluateFandDF(X,x,Rt,Tt,vVec,dT,t)
[f,gradf]=poseEstimationCost(Rt(t),Tt(t),X,x);
df=gradf'*[vVec;dT];

function [df,ddf]=evaluateDFandDdF(X,x,Rt,Tt,vVec,dT,w,t)
[~,gradf,Hf]=poseEstimationCost(Rt(t),Tt(t),X,x);
v=[vVec;dT];
df=gradf'*w;
ddf=v'*Hf*w;

function [df,ddf]=evaluateGradFanddGradF(X,xt,R,T,dx,t)
[~,gradf,~,Jx]=poseEstimationCost(R,T,X,xt(t));
df=gradf;
ddf=Jx*dx;


