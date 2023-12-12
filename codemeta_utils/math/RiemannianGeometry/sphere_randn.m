%function Y=sphere_randn(y,v,N)
%Same as the function lie_randn() specialized for generating random
%vectors on the sphere.
%
%See also lie_randn
function Y=sphere_randn(varargin)
Y=lie_randn(sphere_funs(),varargin{:});
