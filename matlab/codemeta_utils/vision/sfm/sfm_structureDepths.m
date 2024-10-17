function data=POCsfm_structureDepths(data)
[NFeatures,NStructure]=size(data.structure.feature);
depths=zeros(NFeatures,NStructure);
for iStructure=1:NStructure
    P=sfm_getProjectionFromStructureId(data,iStructure);
    for iView=1:size(P,3)
        depths(iView,iStructure)=...
            projectGetDepthsFromP(P(:,:,iView),data.structure.location(:,iStructure));
    end
end

data.structure.depths=depths;
