%Match pairs of feature groups
%function data=sfm_matchExtract(data,varargin)
%Optional inputs
%   'matchIdxImg',m     [2 x NMatches] array with pairs of image indeces to match
%                       (default, data.matchIdxImg, or, if not present, all
%                       pairs)
%   'showstats'         Display info as each pair is processed
%   'flagSymmetryValidation',f  If true, validate matches by intersecting
%       matches in both directions (1->2 and 2->1). Default: true;
%   'optsMatch',opts    Cell array of options to be passed to vl_ubcmatch
%Input data fields 
%   data.feature.descriptor     descriptors used for matching
%Output data fields
%   data.matchIdxImg    [2 x NMatches] array with pairs of image indeces matched
%   data.matches        struct array with info on each match
%
%Compute matches between pairs of images using vl_ubcmatch

function data=sfm_matchExtract(data,varargin)
flagShowStats=false;
flagSymmetryValidation=true;
optsMatch={};
if isfield(data,'matchIdxImg')
    matchIdxImg=data.matchIdxImg;
else
    NFeature=length(data.feature);
    matchIdxImg=nchoosek(1:NFeature,2)';
end

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'showstats'
            flagShowStats=true;
        case 'optsmatch'
            ivarargin=ivarargin+1;
            optsMatch=[optsMatch{:} varargin{ivarargin}];
        case 'matchidximg'
            ivarargin=ivarargin+1;
            matchIdxImg=varargin{ivarargin};
        case 'flagsymmetryvalidation'
            ivarargin=ivarargin+1;
            flagSymmetryValidation=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if ~exist('matchIdxImg','var') || isempty(matchIdxImg)
end
data.matchIdxImg=matchIdxImg;

NMatch=size(matchIdxImg,2);
match=repmat(struct('idxImg',[],'idxMatch',[],'scores',[]),1,NMatch);
for iMatch=1:NMatch
    idxImg=matchIdxImg(:,iMatch);
    match(iMatch).idxImg=idxImg;
    if flagShowStats
        fprintf('Matching feature groups %d and %d\n',idxImg);
    end
    d1=data.feature(idxImg(1)).descriptor;
    d2=data.feature(idxImg(2)).descriptor;
    [m12,s12]=vl_ubcmatch(d1,d2,optsMatch{:});
    if flagSymmetryValidation
        [m21]=vl_ubcmatch(d2,d1,optsMatch{:});
        [m,idxm12]=sfm_rawIntersectMatch(m12,flipud(m21));
        if flagShowStats
            fprintf('\tPercentage of retained match by symmetry: %.3f%% %.3f%%\n',...
                length(idxm12)/size(m12,2)*100,length(idxm12)/size(m21,2)*100);
        end
        s=s12(idxm12);
    else
        m=m12;
        s=s12;
    end
    match(iMatch).idxMatch=m;
    match(iMatch).scores=s;
end

data.match=match;
