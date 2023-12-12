%A simple-minded nearest neighbor graph from pairwise distances using sorting
%function [idxList,kList]=nearestNeighborGraphFromD(D,k,varargin)
%Inputs
%   D   [NPoints x NPoints] matrix with pair-wise distances
%   k   number of neighbors or epsilon to use
%Optional arguments
%   'k'         compute the k-NN graph (default)
%   'epsilon'   compute the epsilon-ball NN graph
%Output
%   idxList     [NPoints x maxNNeighbors] list of indeces of neighbors for
%               each point
%   kList       number of neighbors for each point
%The algorithm uses a naive approach where each row of D is sorted and the
%first k indeces are used
function [idxList,kList]=nearestNeighborGraphFromD(D,k,varargin)
NPoints=size(D,1);
mode='k';
flagSymmetric=false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case {'k','epsilon'}
            mode=lower(varargin{ivarargin});
        case 'flagsymmetric'
            flagSymmetric=varargin{ivarargin};
        case 'symmetric'
            flagSymmetric=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

%make the diagonal contains Inf, so that a point is not neighbor of itself
D=D+diag(inf(NPoints,1)); 

switch mode
    case 'k'
        [~,idx]=sort(D);
        idxList=idx(1:k,:)'; 
        kList=k*ones(NPoints,1);
    case 'epsilon'
        A=D<k;
        [idxList,kList]=nearestNeighborGraphFromAdjacency(A);
end

if flagSymmetric
    A=sparse(NPoints,NPoints);
    for iPoint=1:NPoints
        A(iPoint,idxList(iPoint,1:kList(iPoint)))=1;
    end
    A=(A+A')>0;
    [idxList,kList]=nearestNeighborGraphFromAdjacency(A);
end
