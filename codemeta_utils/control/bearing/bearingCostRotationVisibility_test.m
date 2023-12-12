function bearingCostRotationVisibility_test
%resetRands(1)
y0=[1;0];
NLandmarks=5;
offset=5*[1;1];
xLandmarks=randn(2,NLandmarks)+offset*ones(1,NLandmarks);
x0=[0;0];
%rMax=12.5;

funs=bearingCostFunctions('cosine');

S=[0 -1; 1 0];
R0=eye(2);

dRVec=rand;
R=@(t) rot_exp(R0,t*dRVec*S);
[T,~,~,dT]=real_randGeodFun(x0,'speed',rand);
t=linspace(-pi,pi);
ny=@(t) bearingComputeRanges(T(t),xLandmarks);
rMax=max(max(evalfun(ny,t)));

figure(1)
check_der(@(t) costAndDer(R(t),dRVec,T(t),dT,xLandmarks,y0,funs),'function',t)

function [phi,dPhi]=costAndDer(R,dRVec,x,dx,xLandmarks,y0,funs)
[y,ny]=bearingCompute(x,xLandmarks);
[phi,gradPhiVec]=bearingCostRotationVisibility(R,y,ny,y0,funs);

dPhi=gradPhiVec'*[dRVec;dx];
