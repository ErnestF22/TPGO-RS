%function S=rot_randn(R,v,N)
%Same as the function lie_randn() specialized for generating random
%rotation matrices.
%
%See also lie_randn
function S=rot_randn(varargin)
S=lie_randn(rot_funs(),varargin{:});
