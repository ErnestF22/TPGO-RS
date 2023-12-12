%Compute nearest neighbor graph from coordinate points
%function  [idxList,kList]=nearestNeighborGraphFromPoints(x,varargin)
%Compute pairwise distance matrix for the points in x and then calls
%nearestNeighborGraphFromD(D,varargin{:})
function  [idxList,kList]=nearestNeighborGraphFromPoints(x,varargin)
D=euclideanDistMatrix(x);
[idxList,kList]=nearestNeighborGraphFromD(D,varargin{:});

