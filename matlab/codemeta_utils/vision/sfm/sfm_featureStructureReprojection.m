function data=sfm_featureStructureReprojection(data)
NFeatures=length(data.feature);
for iFeature=1:NFeatures
    x=zeros(size(data.feature(iFeature).location));
    idxStructure=data.feature(iFeature).structureFilteredMembership;
    flagStructure=idxStructure>0;
    P=data.projection(:,:,iFeature);
    x(:,flagStructure)=projectFromP(P,data.structureFiltered.location(:,idxStructure(flagStructure)));
    data.feature(iFeature).locationReprojected=x;
end

data=sfm_featureNormalize(data,'memberLocation','locationReprojected',...
    'memberLocationNormalized','locationReprojectedNormalized');
