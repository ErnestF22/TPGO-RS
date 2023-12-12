function A=spblkdiag(varargin)
for iArg=1:nargin
    varargin{iArg}=sparse(varargin{iArg});
end
A=blkdiag(varargin{:});
