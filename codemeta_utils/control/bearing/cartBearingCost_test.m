function cartBearingCost_test
%resetRands();
NX=10;
dimX=2;
X=rand(dimX,NX);
XGoal=randn(dimX,1);
XEval0=randn(dimX,1);
dXEval=cnormalize(randn(dimX,1));
XEval=@(t) XEval0+t*dXEval;

subplot(2,1,1)
check_der(@(t) cost(XEval(t),X,XGoal,dXEval),'function')
% subplot(2,1,2)
% check_der(@(t) dcost(XEval(t),X,XGoal,dXEval),'function')


function [c,dc]=cost(XEval,X,XGoal,dXEval)
[c,gradc]=cartBearingCost(XEval,X,XGoal);
dc=gradc'*dXEval;

function [dc,ddc]=dcost(XEval,X,XGoal,dXEval)
[~,gradc,Hessc]=cartBearingCost(XEval,X,XGoal);
dc=gradc'*dXEval;
ddc=dXEval'*Hessc'*dXEval;
