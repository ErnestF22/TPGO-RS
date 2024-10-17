function [u,lambda]=bearingCluster_getBearingsScalesFromA(x,A,varargin)
B=adj2incmatrix(A,'oriented');
[u,lambda]=bearingCluster_getBearingsScalesFromB(x,B,varargin{:});
