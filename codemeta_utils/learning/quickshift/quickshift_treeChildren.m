%Return cell array with lists of children of each node (only one level)
%function treeChildren=quickshift_treeChildren(treeEdges,varargin)
%Input arguments
%   treeEdges   Vector containing the tree
%Optional arguments
%   'selfreference'   For roots, include the node itself in the as a children
function treeChildren=quickshift_treeChildren(treeEdges,varargin)
flagSelfReference=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'selfreference'
            flagSelfReference=true;
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

NEdges=length(treeEdges);
treeChildren=cell(1:NEdges);
for iEdge=1:NEdges
    idxParent=treeEdges(iEdge);
    if iEdge~=idxParent || flagSelfReference
        treeChildren{idxParent}=[treeChildren{idxParent} iEdge];
    end
end
