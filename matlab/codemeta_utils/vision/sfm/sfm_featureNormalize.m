%Normalize feature using the calibration matrices
%function data=sfm_featureNormalize(data)
%Requires data.calibration containing the transformation for each feature
%set. Applies the transformation to data.feature(iImage).location and
%stores the result in data.feature(iImage).locationNormalized
function data=sfm_featureNormalize(data,varargin)
memberLocation='location';
memberLocationNormalized='locationNormalized';
feature=data.feature;
calibration=data.calibration;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'memberlocation'
            ivarargin=ivarargin+1;
            memberLocation=varargin{ivarargin};
        case 'memberlocationnormalized'
            ivarargin=ivarargin+1;
            memberLocationNormalized=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NFeature=length(feature);
for iFeature=1:NFeature
    feature(iFeature).(memberLocationNormalized)=...
        sfm_rawTransformFeature(calibration(:,:,iFeature),...
            feature(iFeature).(memberLocation),'invert');
end
data.feature=feature;