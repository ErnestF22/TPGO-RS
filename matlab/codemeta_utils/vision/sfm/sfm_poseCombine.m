%Combine the estimated relative poses into absolute poses
%function data=sfm_poseCombine(data,varargin)
%Estimate absolute poses from relative poses.
%Optional inputs
%   'optsRotations',c   Cell array of options to pass to
%           sfm_rawRotationsCombineFactorization
%
%See also sfm_rawRotationsCombineFactorization
function data=sfm_poseCombine(data,varargin)
memberNameMatch='matchFiltered';
memberNamePose='matchPoseEstimated';
optsRotations={};
methodTranslations='essential';

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'membernamematch'
            ivarargin=ivarargin+1;
            memberNameMatch=varargin{ivarargin};
        case 'membernamepose'
            ivarargin=ivarargin+1;
            memberNamePose=varargin{ivarargin};
        case 'optsrotations'
            ivarargin=ivarargin+1;
            optsRotations=[optsRotations varargin{ivarargin}];
        case 'methodtranslations'
            ivarargin=ivarargin+1;
            methodTranslations=lower(varargin{ivarargin});
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

[Rij,Tij]=G2RT(data.(memberNamePose));
E=(sfm_getMatchIdxImg(data,memberNameMatch))';

Ri=sfm_rawRotationsCombine(Rij,E,optsRotations{:});

switch methodTranslations
    case 'essential'
        [x1,x2]=sfm_getFeatureLocationFromMatchId(data,[],'normalized',...
            'member',memberNameMatch);
        Ti=sfm_rawAverageTranslationsFromRotationsAndPoints(Ri,x1,x2,E);
    case 'none'
        Ti=zeros(3,size(Ri,3));
end
data.poseEstimated=RT2G(Ri,Ti);
