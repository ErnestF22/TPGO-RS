function bearingComputeSecondDerivative_test
d=3;
N=2;
X=randn(d,N);
v0=randn(d,1);
v1=randn(d,1);
v2=randn(d,1);

x=@(t) v0+t*v1+0.5*t^2*v2;
dx=@(t) v1+t*v2;
ddx=@(t) v2;

%check_der(@(t) funAndDer(x(t),dx(t),X),'function')
check_der(@(t) derAndDder(x(t),dx(t),ddx(t),X),'function')


function [y,dy]=funAndDer(x,dx,X)
[y,ny]=bearingCompute(x,X);
dy=bearingComputeDerivative(dx,y,ny);

function [dy,ddy]=derAndDder(x,dx,ddx,X)
[y,ny]=bearingCompute(x,X);
dy=bearingComputeDerivative(dx,y,ny);
ddy=bearingComputeSecondDerivative(dx,ddx,y,ny);
