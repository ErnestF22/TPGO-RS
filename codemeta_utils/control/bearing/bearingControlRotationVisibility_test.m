function bearingControlRotationVisibility_test
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
nyMax=max(max(evalfun(ny,t)));

figure(1)
plotfun(@(t) derAlongControl(R(t),T(t),xLandmarks,y0,funs,nyMax),t)

function dPhi=derAlongControl(R,x,xLandmarks,y0,funs,nyMax)
[y,ny]=bearingCompute(x,xLandmarks);
[~,gradPhiVec]=bearingCostRotationVisibility(R,y,ny,y0,funs);
uVec=bearingControlRotationVisibility(R,y,nyMax,y0,funs);

dPhi=gradPhiVec'*uVec;
