function POCsfm_structureAssign
load sfmdata_castle

threshold=1.5;
NFeatures=length(data.feature);
X=data.structureFiltered.location;
for iFeature=1:NFeatures
    P=data.projection(:,:,iFeature);
    xReproj=projectFromP(P,X);
    x=data.feature(iFeature).location;
    d=euclideanDistMatrix(xReproj,x);
    [dMin,idxMin]=min(d);
    flagValid=dMin<threshold;
    data.feature(iFeature).location=xReproj(:,flagValid);
    data.feature(iFeature).structureFilteredMembership=idxMin(flagValid);
end

data=sfm_featureNormalize(data);

sfm_displayFeatureStructureReprojection(data)
