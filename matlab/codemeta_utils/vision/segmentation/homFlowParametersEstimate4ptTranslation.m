%Estimate plane normal and camera linear velocity from flow and angular velocity
%function homFlowParametersEstimate4pt(x,dx,w)
function [n,v,A,b]=homFlowParametersEstimate4ptTranslation(x,dx,w,varargin)
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
            nPrior=varargin{ivarargin}(1:3);
            vPrior=varargin{ivarargin}(4:6);
        case 'lambda'
            ivarargin=ivarargin+1;
            lambda=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

[A,b]=homFlowParametersEstimate4ptTranslationSystemMat(x,w);

if ~flagPrior
    nv=reshape(pinv(A)*(dx(:)-b),3,3);
    s=svd(nv);
    nv=nv-s(2)*eye(3);
else
    nvPrior=nPrior*vPrior';
    nv=(lambda(1)*(A'*A)+lambda(2)*eye(9))...
        \(lambda(1)*(A'*(dx(:)-b))+lambda(2)*nvPrior(:));
    nv=reshape(nv,3,3);
end
[U,S,V]=svd(nv);
n=S(1,1)*U(:,1);
v=V(:,1);
if flagPrior
    %use prior to pick the sign
    if v'*vPrior<0
        v=-v;
        n=-n;
    end
end
