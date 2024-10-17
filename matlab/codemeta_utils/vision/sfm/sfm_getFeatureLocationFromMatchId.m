%Get feature locations from a match id
%function [x1,x2,m12]=sfm_getFeatureLocationFromMatchId(data,iMatch,varargin)
%Get indeces of images from data.matchIdxImg(:,iMatch), get feature
%locations for those images from data.feature(idxImg()).(locationMemberName),
%get indeces m12 of the matching features from data.(matchMemberName)(iMatch).idxMatch, 
%and use it to get the matching locations.
%Inputs
%   data
%   iMatch  index (id) of the match. If empty, implies all matches.
%Optional arguments
%   'locationMemberName',locationMemberName
%   'member',matchMemberName
%Outputs
%If iMatch contains a single element:
%   x1,x2   [2 x NMatches] array of corresponding feature points
%   m12     [2 x NMatches] array of indeces in the original arrays of
%           corresponding locations
%If iMatch contains multiple elements, the same as above, but x1,x2,m12 are
%cell arrays of arrays.
function [x1,x2,m12]=sfm_getFeatureLocationFromMatchId(data,idxMatch,varargin)
locationMemberName='location';
matchMemberName='match';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'normalized'
            locationMemberName='locationNormalized';
        case 'member'
            ivarargin=ivarargin+1;
            matchMemberName=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

matchIdxImg=sfm_getMatchIdxImg(data,matchMemberName);

if ~exist('idxMatch','var') || isempty(idxMatch)
    idxMatch=1:size(matchIdxImg,2);
end

NIdxMatch=length(idxMatch);
idxImg=matchIdxImg(:,idxMatch);

if NIdxMatch==1
    [x1,x2,m12]=getFeaturesSingle(data,idxMatch,locationMemberName,matchMemberName,idxImg);
else
    x1=cell(NIdxMatch,1);
    x2=cell(NIdxMatch,1);
    m12=cell(NIdxMatch,1);
    
    for iIdxMatch=1:NIdxMatch
        [x1{iIdxMatch},x2{iIdxMatch},m12{iIdxMatch}]=getFeaturesSingle(data,idxMatch(iIdxMatch),locationMemberName,matchMemberName,idxImg(:,iIdxMatch));
    end
end

function [x1,x2,m12]=getFeaturesSingle(data,iMatch,locationMemberName,matchMemberName,idxImg)
x1=data.feature(idxImg(1)).(locationMemberName);
x2=data.feature(idxImg(2)).(locationMemberName);
m12=data.(matchMemberName)(iMatch).idxMatch;
x1=x1(:,m12(1,:));
x2=x2(:,m12(2,:));
