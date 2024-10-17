function RMean=rot_mean(R,varargin)
RMean=lie_mean(R, rot_funs(), varargin{:});