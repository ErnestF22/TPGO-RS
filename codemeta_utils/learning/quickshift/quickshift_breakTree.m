%Break edges in a tree by comparing lengths with a fixed threshold
function treeEdges=quickshift_breakTree(treeDistances,treeEdges,varargin)
threshold=0.01;
%type of pairwise relation contained in D
relationType='distance';

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'threshold'
            ivarargin=ivarargin+1;
            threshold=varargin{ivarargin};
        case 'relationtype'
            ivarargin=ivarargin+1;
            relationType=lower(varargin{ivarargin});
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

switch relationType
    case 'distance'
        flagRoots=treeDistances>threshold;
    case 'similarity'
        flagRoots=treeDistances<threshold;
    otherwise
        error('Relation type not recognized')
end

idx=1:length(treeDistances);

treeEdges(flagRoots)=idx(flagRoots);
