function bearingComputeDerivative_test
nTask=2;
N=3;
X=randn(2,N);
X0=randn(2,1);
v=randn(2,1);
tXEval=@(t) X0+t*v;

switch nTask
    case 1
        check_der(@(t) bearings(v,tXEval(t),X),'function')
    case 2
        t=linspace(0,1,100);
        subplot(2,1,1)
        plotfun(@(t) sumBearingDerivative(v,tXEval(t),X),t)
        title('Sum of bearing derivatives')
        subplot(2,1,2)
        plotfun(@(t) reconstructVelocity(v,tXEval(t),X),t)
        hold on
        plotfun(@(t) v,t,'rx')
        hold off
        title('Reconstructed velocity versus true one')
end

function [Y,dY]=bearings(v,XEval,X)
[Y,nY]=bearingCompute(XEval,X);
dY=bearingComputeDerivative(v,Y,nY);

function vRec=reconstructVelocity(v,XEval,X)
[dY,Q]=sumBearingDerivative(v,XEval,X);
vRec=Q\dY;

function [dY,Q]=sumBearingDerivative(v,XEval,X)
[Y,nY]=bearingCompute(XEval,X);
[dY,Q]=bearingComputeDerivative(v,Y,nY);
dY=sum(dY,2);
Q=sum(Q,3);
