function P=sfm_getProjectionFromStructureId(data,iStructure)
NFeature=data.structure.featureCount(iStructure);
P=zeros(3,4,NFeature);
for iiFeature=1:NFeature
    iFeature=data.structure.feature(iiFeature,iStructure);
    P(:,:,iiFeature)=data.projection(:,:,iFeature);
end
