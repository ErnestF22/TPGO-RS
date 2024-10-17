%Get all the features corresponding to a given 3-D point
%function [x,iFeatureList]=sfm_getFeatureLocationByStructureId(data,iStructure,varargin)
function [x,iFeatureList]=sfm_getFeatureLocationByStructureId(data,iStructure,varargin)
locationMemberName='location';
structureMemberName='structure';

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

NFeature=data.(structureMemberName).featureCount(iStructure);
iFeatureList=data.(structureMemberName).feature(1:NFeature,iStructure);
x=zeros(2,NFeature);
for iiFeature=1:NFeature
    iFeature=iFeatureList(iiFeature);
    idxFeature=data.(structureMemberName).idxFeature(iiFeature,iStructure);
    x(:,iiFeature)=data.feature(iFeature).(locationMemberName)(:,idxFeature);
end

