function bearingCostRotationVisibilityROnly_test
%resetRands(1)
y0=[1;0];
NLandmarks=5;
offset=5*[1;1];
xLandmarks=randn(2,NLandmarks)+offset*ones(1,NLandmarks);
x0=[0;0];

funs=bearingCostFunctions('cosine');

R0=eye(2);

dRVec=rand;
R=@(t) rot_exp(R0,rot_hat(R0,t*dRVec));
t=linspace(-pi,pi);

figure(1)
check_der(@(t) costAndDer(x0,R(t),dRVec,xLandmarks,y0,funs),'function',t)

function [phi,dPhi]=costAndDer(x0,R,dRVec,xLandmarks,y0,funs)
y=bearingCompute(x0,xLandmarks);
[phi,gradPhiRVec]=bearingCostRotationVisibilityROnly(R,y,y0,funs);

dPhi=gradPhiRVec*dRVec;
