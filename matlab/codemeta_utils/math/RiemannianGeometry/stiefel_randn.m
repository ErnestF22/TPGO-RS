%function Y=stiefel_randn(y,v,N)
%Same as the function lie_randn() specialized for generating random
%orthonormal matrices.
%
%See also lie_randn
function Y=stiefel_randn(varargin)
Y=lie_randn(stiefel_funs(),varargin{:});
