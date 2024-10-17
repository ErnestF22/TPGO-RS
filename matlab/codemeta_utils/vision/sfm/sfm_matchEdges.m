%Returns match list as a list of edges
%function E=sfm_matchEdges(data,varargin)
function E=sfm_matchEdges(data,varargin)
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

E=[data.(memberMatch).idxImg]';
