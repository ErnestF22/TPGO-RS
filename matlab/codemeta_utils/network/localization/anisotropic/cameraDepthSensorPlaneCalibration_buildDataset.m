function [GCamera,GDepth,GPlaneTags,NVecPlanesDepth]=cameraDepthSensorPlaneCalibration_buildDataset(NPlanes,GTag)
NTags=size(GTag,3);
NVec=[0;0;1;0];
methodAbsolutePoses='references';

%pose camera
GCamera=eye(4);
%pose depth sensor
GDepth=RT2G(rot_randn(eye(3),0.05),[0;3;0]);

%pose planes and tags on planes
GPlanes=zeros(4,4,NPlanes);
NPlaneTags=NPlanes*NTags;
GPlaneTags=zeros(4,4,NPlaneTags);
idxPlaneTags=reshape(1:NPlaneTags,NTags,NPlanes);
for iPlane=1:NPlanes
    R=rot_randn(eye(3),1);
    T=R*[0;0;5+2*rand];
    GPlanes(:,:,iPlane)=RT2G(R,T);
    for iTag=1:NTags
        GPlaneTags(:,:,idxPlaneTags(iTag,iPlane))=GPlanes(:,:,iPlane)*GTag(:,:,iTag);
    end
end

%plane equations in world's reference frame
NVecPlanes=squeeze(rigidTransformG(invg(GPlanes),NVec,'planes',methodAbsolutePoses));

%data in depth sensor's reference frame
NVecPlanesDepth=rigidTransformG(GDepth,NVecPlanes,'planes',methodAbsolutePoses,'wc');
