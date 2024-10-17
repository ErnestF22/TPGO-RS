%function [mu,output]=medoids(x,nbClusters,varargin)
%Runs k-medoids
%Inputs
%   x           datapoints
%   nbClusters  number of clusters to find
%Optional inputs
%   'muInit' mu     [d x nbClusters] Initial guess for the cluster centers mu. By default,
%       they are generated at random using a uniform distribution with the
%       same range as the datapoints
%   'bias',lambda   [1 x nbClusters] cell array of [d x nbBias]
%       matrices, where each column contains a vector lambdak such that the
%       term lambdak'*muk is added to the cost. If nbBias>1, the columns
%       are effectively summed
%   'priorCenters',z,rho
%                   [1 x nbClusters] and [1 x nbClusters] cell arrays, 
%                   where each element zk and rhok is a [d x
%                   nbPriorCenters] and [1 x nbPriorCenters]
%                   containing column
%                   vectors (resp., a row vector) of "priors" (resp.,
%                   weights) for the center muk. In practice, terms of the
%                   form rho[c]*norm(z[:,c]-mu,1) are added to the
%                   k-medoids cost.
function [mu,output]=medoids(x,nbClusters,varargin)
flagMuProvided=false;
nbIterationsMax=100;
flagDebug=false;
tolStop=1e-6;
flagBias=false;
flagVerbose=false;
flagPriorCenters=false;
flagStopIncrease=true;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'muinit'
            ivarargin=ivarargin+1;
            mu=varargin{ivarargin};
            flagMuProvided=true;
        case 'bias'
            ivarargin=ivarargin+1;
            bias=varargin{ivarargin};
            flagBias=true;
        case 'priorcenters'
            ivarargin=ivarargin+1;
            priorCenters=varargin{ivarargin};
            ivarargin=ivarargin+1;
            priorCentersWeights=varargin{ivarargin};
            flagPriorCenters=true;
        case 'nbiterationsmax'
            ivarargin=ivarargin+1;
            nbIterationsMax=varargin{ivarargin};
        case 'nostopincrease'
            flagStopIncrease=false;
        case 'debug'
            flagDebug=true;
        case 'verbose'
            flagVerbose=true;
        otherwise
            error('Optional argument not recognized')
    end
    ivarargin=ivarargin+1;
end

if ~flagMuProvided
    xMax=max(x);
    xMin=min(x);
    mu=(rand(1,nbClusters)-xMin).*(xMax-xMin)+xMin;
end

if flagDebug
    dimData=size(mu,1);
    muAll=NaN([dimData,nbClusters, nbIterationsMax+1]);
    muAll(:,:,1)=mu;
    cAll=NaN(1,nbIterationsMax);
end

%setup optional arguments for computing cost
optsAssignmentSearchSingle={};
if flagBias
    optsAssignmentSearchSingle{end+1}='bias';
    optsAssignmentSearchSingle{end+1}=bias;
end
if flagPriorCenters
    optsAssignmentSearchSingle{end+1}='priorCenters';
    optsAssignmentSearchSingle{end+1}=priorCenters;
    optsAssignmentSearchSingle{end+1}=priorCentersWeights;
end

cPrev=Inf;
for it=1:nbIterationsMax
    for iCluster=1:nbClusters
        [mu,c]=medoids_centersAssignmentSearchSingle(x,mu,iCluster,optsAssignmentSearchSingle{:});
        if flagVerbose
            fprintf('it=%d, iCluster=%d, c=%d\n',it,iCluster,c)
        end
    end
    
    if flagDebug
        cAll(it)=c;
        muAll(:,:,it+1)=mu;
    end
    
    if flagStopIncrease && cPrev-c<tolStop
       break
    end
    cPrev=c;
end

output.it=it;
output.nbIterationsMax=nbIterationsMax;
if flagDebug
    output.c=cAll(1:it);
    output.mu=squeeze(muAll(:,:,1:it+1));
end
