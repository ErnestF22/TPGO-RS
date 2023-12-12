function bearingCostGeneral_test
global funs
%resetRands();
%funs=bearingCostFunctions('cosine');
funs=bearingCostFunctions('angleSq');
NX=10;
dimX=2;
XLandmarks=rand(dimX,NX);
XGoal=randn(dimX,1);
XEval0=randn(dimX,1);
[XEval,~,~,dXEval]=real_randGeodFun(XEval0,'speed',1);
YGoal0=bearingCompute(XGoal,XLandmarks);
[YGoal,dYGoal]=sphere_randGeodFun(YGoal0,'speed',1);

figure(1)
subplot(2,1,1)
check_der(@(t) cost(XEval(t),XLandmarks,YGoal0,dXEval),'function')
subplot(2,1,2)
check_der(@(t) dcost(XEval(t),XLandmarks,YGoal0,dXEval),'function')
figure(2)
tdgrad=@(t) dgrad(XEval(t),XLandmarks,YGoal(t),dXEval,dYGoal(t));
check_der(tdgrad,'function')


function [c,dc]=cost(XEval,X,YGoal,dXEval)
global funs
[YEval,nYEval]=bearingCompute(XEval,X);
[c,gradc]=bearingCostGeneral(YEval,YGoal,nYEval,funs);
dc=gradc'*dXEval;

function [dc,ddc]=dcost(XEval,X,YGoal,dXEval)
global funs
[YEval,nYEval]=bearingCompute(XEval,X);
[~,gradc,Dgradc]=bearingCostGeneral(YEval,YGoal,nYEval,funs);
dc=gradc'*dXEval;
ddc=dXEval'*Dgradc'*dXEval;

function [gradc,dgradc]=dgrad(XEval,X,YGoal,dXEval,dYGoal)
%includes measurements
global funs
[YEval,nYEval]=bearingCompute(XEval,X);
[~,gradc,Dgradc,Dygradc]=bearingCostGeneral(YEval,YGoal,nYEval,funs);
dgradc=Dgradc*dXEval;
NY=size(YGoal,2);
for iY=1:NY
    dgradc=dgradc+Dygradc(:,:,iY)*dYGoal(:,iY);
end
