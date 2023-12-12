%Compute desired angle of gradient with horizontal and its derivartive 
%function [theta,gradTheta]=bearingCostGeneral_desiredAngle(YEval,YGoal,nYEval)
%Note: if only theta is computed, nYEval is not necessary
function [theta,gradTheta]=bearingCostGeneral_desiredAngle(YEval,YGoal,funs,nYEval)
flagComputeGrad=false;
if nargout>1
    flagComputeGrad=true;
end

%Note: the desired direction is in the negative of the gradient, so the
%sign of g needs to be flipped
if ~flagComputeGrad
    g=-bearingCostGeneral_gradient(YEval,YGoal,funs);
else
    [g,H]=bearingCostGeneral_gradient(YEval,YGoal,funs,nYEval);
    g=-g;
    gradTheta=H*[0 -1; 1 0]*g/(g'*g);
end
theta=atan2(g(2),g(1));
