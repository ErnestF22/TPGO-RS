%function c=medoids_centersAssignmentCost(x,mu,varargin)
%Compute the medoids cost for given centers, including assignments.
%Inputs
%   xPoints     [1 x nbPoints] data points
%   mu          [1 x nbClusters] centers
%Optional outputs: same as medoids.m
function [c,output]=medoids_centersAssignmentCost(xPoints,mu,varargin)
nbClusters=size(mu,2);
idx=medoids_assign(xPoints,mu);
cCluster=zeros(1,nbClusters);

flagBias=false;
flagPriorCenters=false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'bias'
            ivarargin=ivarargin+1;
            biasMult=varargin{ivarargin};
            flagBias=true;
        case 'priorcenters'
            ivarargin=ivarargin+1;
            priorCenters=varargin{ivarargin};
            ivarargin=ivarargin+1;
            priorCentersWeights=varargin{ivarargin};
            flagPriorCenters=true;
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

for iCluster=1:nbClusters
    optsCenter={};
    dataCluster=xPoints(:,idx==iCluster);
    if flagBias
        optsCenter=[optsCenter 'bias' biasMult{iCluster}];
    end
    if flagPriorCenters
        optsCenter=[optsCenter 'priorcenters'...
            priorCenters{iCluster} priorCentersWeights{iCluster}];
    end
    cCluster(iCluster)=medoids_centerCost(dataCluster,'xEval',mu(:,iCluster),optsCenter{:});
end
c=sum(cCluster);
output.cCluster=cCluster;
output.idx=idx;
