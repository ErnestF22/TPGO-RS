function [c,JRT,HRT]=epipolarCost(R,T,x1,x2)
Nx=size(x1,2);

flagComputeJacobian=false;
flagComputeSecondJacobian=false;
if nargout>1
    flagComputeJacobian=true;
end
if nargout>2
    flagComputeSecondJacobian=true;
end


if ~flagComputeJacobian
    e=epipolarConstraint(R,T,x1,x2);
else
    if ~flagComputeSecondJacobian
        [e,JRTe]=epipolarConstraint(R,T,x1,x2);
    else
        [e,JRTe,HRTe]=epipolarConstraint(R,T,x1,x2);
    end        
end
c=0.5*e.^2;
if flagComputeJacobian
    et=permute(e,[2 3 1]);
    JRT=repmat(et,[1,6,1]).*JRTe;
    if flagComputeSecondJacobian
        HRT=zeros(size(HRTe));
        for ix=1:Nx
            HRT=JRTe(:,:,ix)'*JRTe(:,:,ix)+e(ix)*HRTe(:,:,ix);
        end
    end
end
