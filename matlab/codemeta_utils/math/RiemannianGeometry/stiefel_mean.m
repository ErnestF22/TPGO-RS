function YMean=stiefel_mean(Y,varargin)
YMean=lie_mean2(Y, stiefel_funs(), varargin{:});
