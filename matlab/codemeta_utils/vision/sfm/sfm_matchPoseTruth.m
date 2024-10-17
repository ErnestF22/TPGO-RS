%Add relative poses corresponding to a match
%function data=sfm_addMatchPose(data,varargin)
%Requires field poseTruth. Adds field matchPoseTruth with the relative pose
%If no member with matches is given, then the root matchIdxImg is used.
%Otherwise, it is extracted from the given match field.
%Optional arguments
%   'memberMatch',name  name of field with matches info to use for
%                           extracting image indeces
function data=sfm_matchPoseTruth(data,varargin)
memberMatch=[];
memberAbsolutePoses='poseTruth';
memberRelativePoses='matchPoseTruth';

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'membermatch'
            ivarargin=ivarargin+1;
            memberMatch=varargin{ivarargin};
        case 'memberabsoluteposes'
            ivarargin=ivarargin+1;
            memberAbsolutePoses=varargin{ivarargin};
        case 'memberrelativeposes'
            ivarargin=ivarargin+1;
            memberRelativePoses=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

matchIdxImg=sfm_getMatchIdxImg(data,memberMatch);

NMatch=size(matchIdxImg,2);
matchPoseTruth=zeros(4,4,NMatch);
for iMatch=1:NMatch
    idxImg=matchIdxImg(:,iMatch);
    G1=data.(memberAbsolutePoses)(:,:,idxImg(1));
    G2=data.(memberAbsolutePoses)(:,:,idxImg(2));
    matchPoseTruth(:,:,iMatch)=computeRelativePoseFromG(G1,G2,'references');
end
data.(memberRelativePoses)=matchPoseTruth;
