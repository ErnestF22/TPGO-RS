function sfm_displayFeatureStructureReprojection(data,iFeatureList)
if ~exist('iFeatureList','var') || isempty(iFeatureList)
    NFeatures=length(data.feature);
    iFeatureList=1:NFeatures;
end
NFeatures=length(iFeatureList);
NCols=ceil(sqrt(NFeatures));
NRows=ceil(NFeatures/NCols);

for iiFeature=iFeatureList
    iFeature=iFeatureList(iiFeature);
    img=sfm_getImageById(data,iFeature);
    
    [x,X]=sfm_getFeatureStructureByFeatureId(data,iFeature);
    P=data.projection(:,:,iFeature);
    xReproj=projectFromP(P,X);
    subplot(NCols,NRows,iiFeature)
    sfm_rawDisplayFeature(img,x);
    hold on
    plotPoints(xReproj,{'Color',[0 0.9 0.9],'Marker','o'})
    hold off
end
