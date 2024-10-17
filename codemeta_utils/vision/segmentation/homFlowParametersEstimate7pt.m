%Given point coordinates and flow, estimate homography flow parameters in LS
%function [alpha,A]=homFlowParametersEstimate7pt(x,dx)
function [alpha,A]=homFlowParametersEstimate7pt(x,dx,varargin)
flagPrior=false;
lambda(1)=5;    %cost parameter for data term
lambda(2)=1;    %cost parameter for regularizer (constraint) term

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'prior'
            flagPrior=true;
            ivarargin=ivarargin+1;
            alphaPrior=varargin{ivarargin};
        case 'lambda'
            ivarargin=ivarargin+1;
            lambda=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

A=homFlowParametersEstimate7ptSystemMat(x);

if ~flagPrior
    alpha=A\dx(:);
else
    alpha=(lambda(1)*(A'*A)+lambda(2)*eye(8))...
        \(lambda(1)*(A'*dx(:))+lambda(2)*alphaPrior);
end
