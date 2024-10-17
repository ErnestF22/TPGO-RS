%Passes arguments to fprintf if flag is true, otherwise it does nothing
%function fprintfFlag(flag,varargin)
function fprintfFlag(flag,varargin)
if flag
    fprintf(varargin{:})
end
