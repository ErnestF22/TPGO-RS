%Add correspondences from features to triangulated structure
%function sfm=sfm_addFeatureStructureMembership(data)
%Requires field structureFiltered. Adds field
%feature().structureFilteredMembership
function data=sfm_addFeatureStructureMembership(data)
%alias structures
feature=data.feature;
structure=data.structureFiltered;

%allocate membership structure
NFeatures=length(data.feature);
for iFeature=1:NFeatures
    NPoints=size(feature(iFeature).location,2);
    feature(iFeature).structureFilteredMembership=zeros(1,NPoints);
end

NPoints=size(structure.location,2);
for iPoint=1:NPoints
    NFeatures=structure.featureCount(iPoint);
    for iiFeature=1:NFeatures
        iFeature=structure.feature(iiFeature,iPoint);
        iIdxFeature=structure.idxFeature(iiFeature,iPoint);
        
        feature(iFeature).structureFilteredMembership(iIdxFeature)=iPoint;
    end
end

%store results back into main data structure
data.feature=feature;
data.structureFiltered=structure;
