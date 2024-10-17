%Compute errors on plane (normals and distances) correspondences
%function [eR,eT,gradRVec,gradT,DGradR,DGradRn1,DGradRn2]=poseEstimationFromNormalsDistancesResiduals(R,T,n1,d1,n2,d2,varargin)
%The planes are supposed to have equations n'*x=d. The rigid transformation
%R,T shold map points in reference 2 to reference 1. 
function [eR,eT,gradRVec,gradT,DGradR,DGradRn1,DGradRn2]=poseEstimationFromNormalsDistancesResiduals(R,T,n1,d1,n2,d2,varargin)
flagComputeGrad=false;
flagComputeDGrad=false;
flagComputeDGradMeasurements=false;
if nargout>2
    flagComputeGrad=true;
    if nargout>4
        flagComputeDGrad=true;
        if nargout>5
            flagComputeDGradMeasurements=true;
        end
    end
end

methodR='angle';
%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'methodr'
            ivarargin=ivarargin+1;
            methodR=lower(varargin{ivarargin});
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NPlanes=size(n1,2);
eR=zeros(NPlanes,1);
eT=zeros(NPlanes,1);
for iPlane=1:NPlanes
    switch methodR
        case 'angle'
            eR(iPlane)=subspace(n1(:,iPlane),R*n2(:,iPlane));
        case 'cosine'
            eR(iPlane)=1-max(-1,min(1,n1(:,iPlane)'*R*n2(:,iPlane)));
    end
    eT(iPlane)=n1(:,iPlane)'*T-(d1(iPlane)-d2(iPlane));
end
if flagComputeGrad
    if ~strcmp(methodR,'cosine')
        error('Gradient for angle residual not implemented.')
    end
    gradRVec=zeros(3,NPlanes);
    gradT=n1;
    for iPlane=1:NPlanes
        gradRVec(:,iPlane)=-hat(n2(:,iPlane))*R'*n1(:,iPlane);
    end
    if flagComputeDGrad
        DGradR=zeros(3,3,NPlanes);
        for iPlane=1:NPlanes
            DGradR(:,:,iPlane)=-hat(n2(:,iPlane))*hat(R'*n1(:,iPlane));
        end
        if flagComputeDGradMeasurements
            DGradRn1=zeros(3,3,NPlanes);
            DGradRn2=zeros(3,3,NPlanes);
            for iPlane=1:NPlanes
                DGradRn1(:,:,iPlane)=-hat(n2(:,iPlane))*R'*orthComplementProjector(n1(:,iPlane));
                DGradRn2(:,:,iPlane)=hat(R'*n1(:,iPlane))*orthComplementProjector(n2(:,iPlane));
            end
        end
    end
end
