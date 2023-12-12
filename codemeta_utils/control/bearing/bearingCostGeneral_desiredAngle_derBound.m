%Returns bound on the second derivative of the cost and the gradient of
%the cost
%function [bdtheta,g]=bearingCostGeneral_desiredAngle_derBound(YEval,YGoal,funs,dXEval,dmax,type)
%Inputs
%   type    Type of bound, can be 'ub' (upper bound) or 'lb' (lower bound).
%           If omitted, bdtheta contains both
%Outputs
%   bdtheta Bound on the second derivative of the cost. If the input type
%           is omitted, bdtheta(1) contains the ub and bdtheta(2) contains
%           the lb.
function [bdtheta,g]=bearingCostGeneral_desiredAngle_derBound(YEval,YGoal,funs,dXEval,dmax,type)
NY=size(YEval,2);
c=bearingComputeCosine(YEval,YGoal);
f=funs.f;
df=funs.df;
ddf=funs.ddf;

g=zeros(2,1);
H=zeros(2,2,NY);
a=zeros(1,NY);

for iX=1:NY
    Yei=YEval(:,iX);
    Ygi=YGoal(:,iX);
    ci=c(iX);
    fci=f(ci);
    dfci=df(ci);
    ddfci=ddf(ci);
    [gradi,H(:,:,iX)]=bearingCostGeneral_terms(Yei,Ygi,fci,dfci,ddfci,ci);
    g=g+gradi;
end

v=[0 -1; 1 0]*g/(g'*g);

for iX=1:NY
    a(iX)=v'*H(:,:,iX)*dXEval;
end

if ~exist('type','var')
    bdtheta=dmax*[sum(a(a<0)) sum(a(a>0))];
else
    switch type
        case 'lb'
            bdtheta=dmax*sum(a(a<0));
        case 'ub'
            bdtheta=dmax*sum(a(a>0));
    end
end