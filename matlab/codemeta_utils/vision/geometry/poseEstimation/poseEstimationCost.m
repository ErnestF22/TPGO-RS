function [c,gradc,Hessc,Jx]=poseEstimationCost(R,T,X,x,varargin)
flagComputeGradient=false;
flagComputeHessian=false;
flagComputeJacobianX=false;
flagApproximatedHessian=false;

if nargout>1
    flagComputeGradient=true;
    if nargout>2
        flagComputeHessian=true;
        if nargout>3
            flagComputeJacobianX=true;
        end
    end
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


NX=size(X,2);
if ~flagComputeGradient
    xp=projectFromRT(R,T,X,'references');
else
    if ~flagComputeHessian || flagApproximatedHessian
        [xp,Jxp]=projectFromRT(R,T,X,'references');
    else
        [xp,Jxp,Hxp]=projectFromRT(R,T,X,'references');
    end
end
e=xp-x;
c=sum(e.^2)/2/NX;

if flagComputeGradient
    gradc=zeros(6,NX);
    for iX=1:NX
        gradc(:,iX)=Jxp(:,:,iX)'*e(:,iX)/NX;
    end
end
if flagComputeHessian
    Hessc=zeros(6,6,NX);
    for iX=1:NX
        Hc=Jxp(:,:,iX)'*Jxp(:,:,iX);
        if ~flagApproximatedHessian
            for k=1:2
                Hc=Hc+e(k)*Hxp(:,:,k,iX);
            end
        end
        Hessc(:,:,iX)=Hc/NX;
    end
end
if flagComputeJacobianX
    Jx=-permute(Jxp,[2 1 3])/NX;
end
