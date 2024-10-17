%Return the images and the 3-D points available for a feature group
%function [x,X]=sfm_getFeatureStructureCorrespondenceByFeatureId(data)
%Requires field feature().structureFilteredMembership
function [x,X]=sfm_getFeatureStructureByFeatureId(data,iFeature,varargin)
locationMemberName='location';
structureMemberName='structureFiltered';
%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'normalized'
            locationMemberName='locationNormalized';
        case 'member'
            ivarargin=ivarargin+1;
            structureMemberName=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end
feature=data.feature(iFeature);
structure=data.(structureMemberName);

flagFeatureStructure=feature.structureFilteredMembership>0;
iStructureList=feature.structureFilteredMembership(flagFeatureStructure);
x=feature.(locationMemberName)(:,flagFeatureStructure);
X=structure.location(:,iStructureList);
