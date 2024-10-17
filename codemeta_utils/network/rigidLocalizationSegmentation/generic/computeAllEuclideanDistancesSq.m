% Computes all the distances between any two pairs of vectors
% function dist=computeAllEuclideanDistancesSq(x)
% Inputs:
%   x       [D x N] matrix with N vectors
% Outputs
%   dist    [N x N] matrix with all the squared distances
%
function dist=computeAllEuclideanDistancesSq(x)
N=size(x,2);
xNorms=ones(N,1)*sum(x.^2,1);
kernel=x'*x;
dist=xNorms+xNorms'-2*kernel;


