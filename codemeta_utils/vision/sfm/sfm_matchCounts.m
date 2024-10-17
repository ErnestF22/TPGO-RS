%Returns the number of correspondences in each match
function w=sfm_matchCounts(data,varargin)
memberMatch='matchFiltered';

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'membermatch'
            ivarargin=ivarargin+1;
            memberMatch=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

w=cellfun(@(x) size(x,2),{data.(memberMatch).idxMatch})';
