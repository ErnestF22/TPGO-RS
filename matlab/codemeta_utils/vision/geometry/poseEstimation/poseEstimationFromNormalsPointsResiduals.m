function [e,gradRVec,gradT,DGradR,DGradRn,DGradRx]=poseEstimationFromNormalsPointsResiduals(R,T,n,d,x2D)
flagComputeGrad=false;
flagComputeDGrad=false;
flagComputeDGradMeasurements=false;
if nargout>2
    flagComputeGrad=true;
    if nargout>3
        flagComputeDGrad=true;
        if nargout>4
            flagComputeDGradMeasurements=true;
        end
    end
end

NPoints=size(x2D,2);
x=[zeros(1,NPoints); x2D];
e=(sum(n.*(R*x))+T'*n-d)';
if flagComputeGrad
    NPlanes=size(n,2);
    gradRVec=zeros(3,NPlanes);
    gradT=n;
    for iPlane=1:NPlanes
        gradRVec(:,iPlane)=hat(x(:,iPlane))*R'*n(:,iPlane);
    end
    if flagComputeDGrad
        DGradR=zeros(3,3,NPlanes);
        for iPlane=1:NPlanes
            DGradR(:,:,iPlane)=hat(x(:,iPlane))*hat(R'*n(:,iPlane));
        end
        if flagComputeDGradMeasurements
            DGradRn=zeros(3,3,NPlanes);
            DGradRx=zeros(3,2,NPlanes);
            for iPlane=1:NPlanes
                DGradRn(:,:,iPlane)=hat(x(:,iPlane))*R'*orthComplementProjector(n(:,iPlane));
                M=-hat(R'*n(:,iPlane));
                DGradRx(:,:,iPlane)=M(:,2:3);
            end
        end
    end
end