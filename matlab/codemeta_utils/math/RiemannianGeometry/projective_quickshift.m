%Run quickshift clustering on projective space
%function [membership,info]=sphere_quickshift(X,varargin)
%The distance matrix is computed using the Riemannian distance, based on
%the distance on the sphere. The optional arguments are forwarded to
%quickshift_cluster 
function [membership,info]=projective_quickshift(X,varargin)
D=projective_dist(X,X);
[membership,info]=quickshift_cluster(D,varargin{:});
