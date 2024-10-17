%Extracts feature information from all images and puts them all together
%for joint matching
function dataMatch=sfm_featureExtractJointMatchData(data,varargin)
NImages=length(data.feature);
%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'nimages'
            ivarargin=ivarargin+1;
            NImages=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NFeatures=arrayfun(@(x) size(x.location,2),data.feature(1:NImages));

descriptor=[data.feature(1:NImages).descriptor];
locationNormalized=[data.feature(1:NImages).locationNormalized];
location=[data.feature(1:NImages).location];
membershipPrior=zeros(1,sum(NFeatures));
cnt=1;
for iImage=1:NImages
    cntNext=cnt+NFeatures(iImage);
    membershipPrior(cnt:cntNext-1)=iImage;
    cnt=cntNext;
end

dataMatch.descriptor=descriptor;
dataMatch.location=location;
dataMatch.locationNormalized=locationNormalized;
dataMatch.membershipPrior=membershipPrior;
dataMatch.NImages=NImages;
