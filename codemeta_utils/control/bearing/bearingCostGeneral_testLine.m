function bearingCostGeneral_testLine
global funs;
%costName='cosine';
costName='angleSq';

funs=bearingCostFunctions(costName);

%resetRands();
X=zeros(2,1);
XGoal=[1;0];
%testName='limit';
%testName='cost';
testName='offset';
switch testName
    case 'der'
        XEval0=cnormalize(randn(2,1));
    case {'limit','cost','offset'}
        XEval0=XGoal;
end
dXEval=cnormalize(randn(2,1));
%dXEval=cnormalize([-1;0]);
XEval=@(t) XEval0+t*dXEval;

t=linspace(0,20);
fInf=funs.f(dXEval(1));
figure(1)
switch testName
    case 'der'
        subplot(2,1,1)
        check_der(@(t) cost(XEval(t),X,XGoal,dXEval),'function',t)
        subplot(2,1,2)
        check_der(@(t) gradCost(XEval(t),X,XGoal,dXEval),'function',t)
    case 'limit'
        plotfun(@(t) dcost(XEval(t),X,XGoal,dXEval),t)
        hold on
        plotfun(@(t) fInf,t,'r')
        hold off
    case 'cost'
        plotfun(@(t) cost(XEval(t),X,XGoal,dXEval),t)
        hold on
        plotfun(@(t) fInf*t,t,'r')
        hold off
    case 'offset'
        r1=sqrt(dXEval(1)^2+dXEval(2));
        g=@(t) cost(XEval(t),X,XGoal,dXEval)/fInf-t;
        plotfun(g,t)
        hold on
        plotfun(@(t) r1,t,'r')
        hold off
end

function [c,dc]=cost(XEval,X,XGoal,dXEval)
global funs;
[YEval,YGoal,nYEval]=bearingComputeOld(XEval,X,XGoal);
[c,gradc]=bearingCostGeneral(YEval,YGoal,nYEval,funs);
dc=gradc'*dXEval;

function dc=dcost(XEval,X,XGoal,dXEval)
[~,dc]=cost(XEval,X,XGoal,dXEval);

function [gradc,dgradc]=gradCost(XEval,X,XGoal,dXEval)
global funs;
[YEval,YGoal,nYEval]=bearingComputeOld(XEval,X,XGoal);
[~,gradc,Dgradc]=bearingCostGeneral(YEval,YGoal,nYEval,funs);
dgradc=Dgradc*dXEval;
