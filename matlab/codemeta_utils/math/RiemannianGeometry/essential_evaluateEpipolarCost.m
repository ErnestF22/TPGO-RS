function [c,JQ,HQ]=essential_evaluateEpipolarCost(Q,x1,x2,varargin)
flagComputeJacobian=false;
flagComputeSecondJacobian=false;
flagApproximatedHessian=false;

if nargout>1
    flagComputeJacobian=true;
end
if nargout>2
    flagComputeSecondJacobian=true;
end

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'flagapproximatedhessian'
            ivarargin=ivarargin+1;
            flagApproximatedHessian=varargin{ivarargin};
        case 'approximatedhessian'
            flagApproximatedHessian=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end


Nx=size(x1,2);

if ~flagComputeJacobian
    e=essential_evaluateEpipolarConstraint(Q,x1,x2);
else
    if ~flagComputeSecondJacobian || flagApproximatedHessian
        [e,JQe]=essential_evaluateEpipolarConstraint(Q,x1,x2);
    else
        [e,JQe,HQe]=essential_evaluateEpipolarConstraint(Q,x1,x2);
        HQ=zeros(size(HQe));
    end
    JQ=zeros(size(JQe));
    for ix=1:Nx
        JQ(:,ix)=e(ix)*JQe(:,ix);
        if flagComputeSecondJacobian
            if ~flagApproximatedHessian
                HQ(:,:,ix)=JQe(:,ix)*JQe(:,ix)'+e(ix)*HQe(:,:,ix);
            else
                HQ(:,:,ix)=JQe(:,ix)*JQe(:,ix)';
            end
        end
    end
end
c=0.5*e.^2;
