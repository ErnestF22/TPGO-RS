%Given a distance matrix, divide the points in clusters
function [membershipCluster,info]=quickshift_cluster(D,varargin)
phi=@(x) exp(-x.^2/2);
methodBreakTree='threshold';
optsBreakTree={};
optsTree={};
optsDensity={};
flagMembershipPrior=false;
Threshold=0.2;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'phi'
            ivarargin=ivarargin+1;
            phi=varargin{ivarargin};
        case 'optsDensity'
            ivarargin=ivarargin+1;
            optsDensity=varargin{ivarargin};
        case 'gaussian'
            ivarargin=ivarargin+1;
            phi=@(x) exp(-x.^2/(2*varargin{ivarargin}^2));
        case 'optsbreaktree'
            ivarargin=ivarargin+1;
            optsBreakTree=[optsBreakTree varargin{ivarargin}];
            %optsBreakTree={methodBreakTree varargin{ivarargin}};
        case 'optstree'
            ivarargin=ivarargin+1;
            optsTree=[optsTree varargin{ivarargin}];
        case 'methodbreaktree'
            ivarargin=ivarargin+1;
            methodBreakTree=varargin{ivarargin};
        case 'membershipprior'
            ivarargin=ivarargin+1;
            membershipPrior=varargin{ivarargin};
            flagMembershipPrior=true;
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if isempty(optsBreakTree)
    optsBreakTree={methodBreakTree Threshold};
end

%compute the density at each point
treeDensity=quickshift_density(phi,D,optsDensity{:});
if nargout>1
    info.density=treeDensity;
end

%compute the quickshift tree
[treeDistances,treeEdges]=quickshift_tree(treeDensity,D,optsTree{:});

%break tree using distances
switch lower(methodBreakTree)
    case 'threshold'
        treeEdgesClusters=quickshift_breakTree(treeDistances,treeEdges,optsBreakTree{:});
    case 'descendents'
        %treeEdgesClusters=quickshift_breakTreeDescendents(treeDistances,treeEdges,optsBreakTree{:});
        treeEdgesClusters=quickshift_breakTreeDescendentsOrdered(treeDistances,treeEdges,treeDensity,optsBreakTree{:});
    case 'descendentstopdown'
        treeEdgesClusters=quickshift_breakTreeDescendents(treeDistances,treeEdges,optsBreakTree{:});
        
end
if nargout>1
    info.treeEdges=treeEdges;
    info.treeEdgesClusters=treeEdgesClusters;
    info.treeDistances=treeDistances;
end

%break tree using prior membership information
if flagMembershipPrior
    [membershipCluster,treeEdgesPrior,clustersIndicatorsPrior]=...
        quickshift_breakTreePrior(treeDistances,treeEdges,membershipPrior);
    if nargout>1
        info.clustersIndicatorsPrior=clustersIndicatorsPrior;
        info.treeEdgesPrior=treeEdgesPrior;
    end
else
    membershipCluster=quickshift_tree2membership(treeEdgesClusters);
end


