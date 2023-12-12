function cartBearingCost_testLine
%resetRands();
X=zeros(2,1);
XGoal=[randn;0];
XEval0=XGoal;
dXEval=cnormalize(randn(2,1));
%dXEval=cnormalize([1;0.1]);
XEval=@(t) XEval0+t*dXEval;

check_der(@(t) cost(XEval(t),X,XGoal,dXEval),'function',linspace(0,10))

function [c,dc]=cost(XEval,X,XGoal,dXEval)
[c,gradc]=cartBearingCost(XEval,X,XGoal);
dc=gradc'*dXEval;
