%Add relative pose for matching images using the essential matrix
%function data=sfm_essentialPose(data)
function data=sfm_essentialPose(data,varargin)
memberNameEssential='matchFilteredEssential';
matchMemberName='matchFiltered';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'essentialmembername'
            ivarargin=ivarargin+1;
            memberNameEssential=varargin{ivarargin};
        case 'matchmembername'
            ivarargin=ivarargin+1;
            matchMemberName=varargin{ivarargin};
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

NMatch=length(data.(matchMemberName));
matchPose=zeros(4,4,NMatch);

for iMatch=1:NMatch
    [x1,x2]=sfm_getFeatureLocationFromMatchId(data,iMatch,'normalized','member',matchMemberName);
    E=data.(memberNameEssential)(:,:,iMatch);
    matchPose(:,:,iMatch)=epipolarEToG(E,x1,x2);
end

data.matchPoseEstimated=matchPose;
