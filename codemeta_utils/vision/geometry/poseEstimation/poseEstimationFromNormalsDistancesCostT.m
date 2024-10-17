function [cT,gradT,DGradT,DGradTn1,DGradTd1,DGradTd2]=poseEstimationFromNormalsDistancesCostT(T,n1,d1,d2)
flagComputeGrad=false;
flagComputeDGrad=false;
flagComputeDGradMeasurements=false;
if nargout>1
    flagComputeGrad=true;
    if nargout>2
        flagComputeDGrad=true;
        if nargout>3
            flagComputeDGradMeasurements=true;
        end
    end
end

eT=T'*n1-(d1-d2);
cT=0.5*eT.^2;
if flagComputeGrad
    u=ones(3,1);
    gradT=(u*eT).*n1;
    if flagComputeDGrad
        NPlanes=length(eT);
        DGradT=zeros(3,3,NPlanes);
        for iPlane=1:NPlanes
            DGradT(:,:,iPlane)=n1(:,iPlane)*n1(:,iPlane)';
        end
        if flagComputeDGradMeasurements
            DGradTn1=zeros(3,3,NPlanes);
            DGradTd1=-n1;
            DGradTd2=n1;
            for iPlane=1:NPlanes
                DGradTn1(:,:,iPlane)=(n1(:,iPlane)*T'+eT(iPlane)*eye(3))*orthComplementProjector(n1(:,iPlane));
            end
        end            
    end
end
