%Compute the total bearing cost and its gradient
function [phi,gradPhi,HessPhi]=cartBearingCostInnerProduct(XEval,X,XGoal)
flagComputeGrad=false;
flagComputeHess=false;

if nargout>1
    flagComputeGrad=true;
end
if nargout>2
    flagComputeHess=true;
end

[dimX,NX]=size(X);
uVec=ones(1,NX);
YGoal=cnormalize(X-XGoal*uVec);
[YEval,nYEval]=cnormalize(X-XEval*uVec);
c=min(1,max(-1,sum(YGoal.*YEval)));
phi=0.5*sum(nYEval.*((1-c).^2));
% phi=0.5*sum((1-c).^2);
% phi=sum(nYEval);

if flagComputeGrad
    gradPhi=zeros(size(XEval));
    if flagComputeHess
        HessPhi=zeros(size(XEval,1));
    end
    for iX=1:NX
        Yei=YEval(:,iX);
%         nYei=nYEval(iX);
        Ygi=YGoal(:,iX);
%         Jei=(eye(dimX)-(Yei*Yei'))/nYei;
%         dYgi=Jei*Ygi;
        ci=c(iX);
        gradPhi=gradPhi+(-0.5*(1-ci)^2*Yei+(1-ci)*(eye(dimX)-(Yei*Yei'))*Ygi);
%         gradPhi=gradPhi+((1-ci)/nYei*(eye(dimX)-(Yei*Yei'))*Ygi);
%         gradPhi=gradPhi-Yei;
%         if flagComputeHess
%             Hi=-(Yei*dYgi'+dYgi*Yei'+ci*Jei)/nYei;
%             HessPhi=HessPhi+(1/(1-ci^2)-ai*ci/sqrt(1-ci^2)^3)*(dYgi*dYgi')-ai/sqrt(1-ci^2)*Hi;
%         end
    end
end

    
