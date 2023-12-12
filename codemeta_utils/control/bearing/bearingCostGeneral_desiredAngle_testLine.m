function bearingCostGeneral_desiredAngle_testLine
global funs;
costName='cosineSq';
%costName='angleSq';

funs=bearingCostFunctions(costName);

%resetRands();
X=zeros(2,1);
XGoal=[randn;0];
%XEval0=XGoal;
XEval0=cnormalize(randn(2,1));
dXEval=cnormalize(randn(2,1));
%dXEval=cnormalize([1;0.1]);
XEval=@(t) XEval0+t*dXEval;

check_der(@(t) cost(XEval(t),X,XGoal,dXEval),'function',linspace(0,10))

function [theta,dtheta]=cost(XEval,X,XGoal,dXEval)
global funs;
[YEval,YGoal,nYEval]=bearingComputeOld(XEval,X,XGoal);
[theta,gradTheta]=bearingCostGeneral_desiredAngle(YEval,YGoal,funs,nYEval);
dtheta=gradTheta'*dXEval;
