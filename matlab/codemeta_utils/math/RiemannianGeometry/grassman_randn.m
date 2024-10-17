%function Y=grassman_randn(X,v,N)
%Same as the function lie_randn() specialized for generating random
%rotation matrices.
%
%See also lie_randn
function Y=grassman_randn(varargin)
Y=lie_randn(grassman_funs(),varargin{:});
