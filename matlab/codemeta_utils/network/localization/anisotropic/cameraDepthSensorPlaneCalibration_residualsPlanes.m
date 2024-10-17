function [eRot,eTransl]=cameraDepthSensorPlaneCalibration_residualsPlanes(GDepth,GPlaneTags,NVecPlanesDepth)
NPlanes=size(NVecPlanesDepth,2);
NPlaneTags=size(GPlaneTags,3);
NTags=NPlaneTags/NPlanes;
if NTags~=round(NTags)
    error('size(GPlaneTags,3) must be divisible by size(NVecPlanesDepth,2)')
end

NVecPlaneTagsDepth=NVecPlanesDepth(:,reshape(repmat(1:NPlanes,NTags,1),[],1));
[nPlaneTagDepth,dPlaneTagDepth]=planeNVecToNd(NVecPlaneTagsDepth);
[RPlaneTags,TPlaneTags]=G2RT(GPlaneTags);
[RDepth,TDepth]=G2RT(GDepth);
eRot=NaN(1,NPlaneTags);
eTransl=NaN(1,NPlaneTags);
for iPlaneTag=1:NPlaneTags
    r=RPlaneTags(:,3,iPlaneTag);
    eRot(iPlaneTag)=subspace(r,RDepth*nPlaneTagDepth(:,iPlaneTag));
    eTransl(iPlaneTag)=r'*(TDepth-TPlaneTags(:,iPlaneTag))+dPlaneTagDepth(iPlaneTag);
end
