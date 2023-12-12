function bearingCostDisplayOld(X,XGoal,funs,L,NGrid)
if ~exist('NGrid','var')
    NGrid=50;
end
optsX={'MarkerSize',10};
xGrid=linspace(0,L,NGrid);

contourAndGradientGrid(@(XEval) cost(XEval,X,XGoal,funs), @(XEval) gradCost(XEval,X,XGoal,funs),xGrid);
hold on
plot(X(1,:),X(2,:),'r*',optsX{:})
plot(XGoal(1),XGoal(2),'g*',optsX{:})
hold off

function c=cost(XEval,X,XGoal,funs)
[YEval,YGoal,nYEval]=bearingComputeOld(XEval,X,XGoal);
c=bearingCostGeneral(YEval,YGoal,nYEval,funs);

function gradc=gradCost(XEval,X,XGoal,funs)
[YEval,YGoal]=bearingComputeOld(XEval,X,XGoal);
gradc=bearingCostGeneral_gradient(YEval,YGoal,funs);
