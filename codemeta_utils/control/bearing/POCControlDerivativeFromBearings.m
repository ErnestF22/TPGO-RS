function POCControlDerivativeFromBearings
d=3;
N=2;
funsName='cosine';
alpha=1;

xLandmarks=randn(d,N);
v0=randn(d,1);
v1=randn(d,1);
v2=randn(d,1);

x=@(t) v0+t*v1+0.5*t^2*v2;
dx=@(t) v1+t*v2;
ddx=@(t) v2;
xg=x(0);
yg=bearingCompute(xg,xLandmarks);

funs=bearingCostFunctions(funsName);

check_der(@(t) funAndDer(x(t),dx(t),xLandmarks,yg,funs),'function')

function [u,du]=funAndDer(x,dx,xLandmarks,yg,funs)
[y,ny]=bearingCompute(x,xLandmarks);
dy=bearingComputeDerivative(dx,y,ny);
u=bearingControlDirect(y,yg,funs);
du=bearingControlDirectDerivative(y,dy,yg,funs);
