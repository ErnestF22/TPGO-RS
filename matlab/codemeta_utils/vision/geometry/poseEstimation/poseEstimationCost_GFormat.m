function [c,gradc,Hessc]=poseEstimationCost_GFormat(G,X,x,varargin)
flagVec=true;
flagApproximatedHessian=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'hatgrad'
            flagVec=false;
        case 'vecgrad'
            flagVec=true;
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

optsCost={G2R(G),G2T(G),X,x,'flagApproximatedHessian',flagApproximatedHessian};
if nargout==1
    c=poseEstimationCost(optsCost{:});
else
    if nargout==2
        [c,gradc]=poseEstimationCost(optsCost{:});
        if ~flagVec
            gradc=rot3r3_hat(G,gradc);
        end
    else
        [c,gradc,Hessc]=poseEstimationCost(optsCost{:});
        if ~flagVec
            gradc=rot3r3_hat(G,gradc);
        end
    end
end