function xmean=hyperbolic_mean(x,varargin)
xmean=lie_mean(x, hyperbolic_funs(), varargin{:});
