%Run quickshift clustering on the sphere
%function [membership,info]=sphere_quickshift(X,varargin)
%The distance matrix is computed using the Riemannian distance. The
%optional arguments are forwarded to quickshift_cluster
function [membership,info]=sphere_quickshift(X,varargin)
D=sphere_dist(X,X);
[membership,info]=quickshift_cluster(D,varargin{:});
