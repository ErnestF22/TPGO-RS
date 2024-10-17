%Get adjacency matrix implied by matches, weighted by the number of matches
%function A=sfm_matchAdjMatrix(data,varargin)
%Optional Inputs
%   'member',m  field for matches (default: 'matchFiltered')
function A=sfm_matchAdjMatrix(data,varargin)

matchMemberName='matchFiltered';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'member'
            ivarargin=ivarargin+1;
            matchMemberName=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NImages=sfm_getImageNumber(data);

A=zeros(NImages,NImages);
idxMatch=sfm_getMatchIdxImg(data,matchMemberName);
NMatchesPerPair=cellfun(@(x) size(x,2),{data.(matchMemberName).idxMatch});

A(sub2ind([NImages NImages],idxMatch(1,:),idxMatch(2,:)))=NMatchesPerPair;
A=A+A';
