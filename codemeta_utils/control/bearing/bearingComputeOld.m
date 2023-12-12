function [YEval,YGoal,nYEval]=bearingComputeOld(XEval,X,XGoal)
NX=size(X,2);
uVec=ones(1,NX);
YGoal=cnormalize(X-XGoal*uVec);
[YEval,nYEval]=cnormalize(X-XEval*uVec);
