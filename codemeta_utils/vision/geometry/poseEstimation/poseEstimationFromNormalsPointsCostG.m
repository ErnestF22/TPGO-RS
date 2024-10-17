function [c,gradGVec,DGradG,DGradGn,DGradGd,DGradGx]=poseEstimationFromNormalsPointsCostG(G,n,d,x2D)
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

[R,T]=G2RT(G);
NPoints=size(x2D,2);
x=[zeros(1,NPoints); x2D];
e=(sum(n.*(R*x))+T'*n-d)';
c=0.5*e.^2;
if flagComputeGrad
    NPlanes=size(n,2);
    gradRVec=zeros(3,NPlanes);
    gradT=zeros(3,NPlanes);
    gradeR=zeros(3,NPlanes);
    for iPlane=1:NPlanes
        gradeR(:,iPlane)=hat(x(:,iPlane))*R'*n(:,iPlane);
        gradRVec(:,iPlane)=gradeR(:,iPlane)*e(iPlane);
        gradT(:,iPlane)=e(iPlane)*n(:,iPlane);
    end
    gradGVec=[gradRVec;gradT];
    if flagComputeDGrad
        DGradR=zeros(3,3,NPlanes);
        DGradT=zeros(3,3,NPlanes);
        DGradRT=zeros(3,3,NPlanes);
        for iPlane=1:NPlanes
            DGradR(:,:,iPlane)=hat(x(:,iPlane))*hat(R'*n(:,iPlane))*e(iPlane)+gradeR(:,iPlane)*gradeR(:,iPlane)';
            DGradT(:,:,iPlane)=n(:,iPlane)*n(:,iPlane)';
            DGradRT(:,:,iPlane)=gradeR(:,iPlane)*n(:,iPlane)';
        end
        DGradG=[DGradR DGradRT; permute(DGradRT,[2 1 3]) DGradT];
        if flagComputeDGradMeasurements
            DGradGn=zeros(6,3,NPlanes);
            DGradGd=zeros(6,NPlanes);
            DGradGx=zeros(6,2,NPlanes);
            for iPlane=1:NPlanes
                graden=(R*x(:,iPlane)+T);
                DGradRn=e(iPlane)*hat(x(:,iPlane))*R'+gradeR(:,iPlane)*graden';
                DGradTn=e(iPlane)*eye(3)+n(:,iPlane)*graden';
                DGradGn(:,:,iPlane)=[DGradRn;DGradTn];
                DGradGd(:,iPlane)=[-gradeR(:,iPlane);-n(:,iPlane)];
                gradex=R'*n(:,iPlane);
                MR=-hat(R'*n(:,iPlane))*e(iPlane)+gradeR(:,iPlane)*gradex';
                MT=n(:,iPlane)*gradex';
                DGradGx(:,:,iPlane)=[MR(:,2:3);MT(:,2:3)];
            end
        end
    end
end