function POCEpipolarCostHessianRank
load triangulate_test_dataset_datacalibrated
R=R(:,:,2);
T=T(:,2);
xL=x(:,:,2);
xR=x(:,:,1);

[c,J,H]=epipolarConstraint(R,T,xL,xR);
cAll=sum(c);
JAll=sum(J,3);
HAll=sum(H,3);

vNull=[zeros(3,1);T];

disp('JAll*vNull')
disp(JAll*vNull)
disp('vNull''*HAll*vNull')
disp(vNull'*HAll*vNull)

% [U,S,V]=svd(HAll);
% w=V(:,end);
% 
% vVec=w(1:3);
% dT=w(4:6);
% 
% 
% Rt=rot_geodFun(R,rot_hat(R,vVec));
% Tt=@(t) T+t*dT;
% 
% %check_der(@(t) fAndDf(Rt,Tt,xL,xR,vVec,dT,t),'function')
% check_der(@(t) dfAndDdf(Rt,Tt,xL,xR,vVec,dT,t),'function')




function [f,df]=fAndDf(Rt,Tt,xL,xR,vVec,dT,t)
[c,J]=epipolarConstraint(Rt(t),Tt(t),xL,xR);
f=sum(c);
df=sum(J,3)*[vVec;dT];

function [df,ddf]=dfAndDdf(Rt,Tt,xL,xR,vVec,dT,t)
[~,J,H]=epipolarConstraint(Rt(t),Tt(t),xL,xR);
w=[vVec;dT];
df=sum(J,3)*w;
ddf=w'*sum(H,3)*w;
