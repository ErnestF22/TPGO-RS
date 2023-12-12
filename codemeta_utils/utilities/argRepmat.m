%Ensures that different arguments have the same number of rows/columns/etc.
%function N=argRepmat(k,x1,x2,...)
%This function is useful for writing functions for which different
%arguments could either be a single value or multiple values, but the
%operation needs to be performed in parallel on all the values.
%
%Checks dimension k of x1,x2,..., i.e., builds [size(x1,k) size(x2,k) ...]
%If the maximum N is greater than one, the arguments for whick size(x,k)=1
%are repeated in the k-th dimension (all dimensions above k are assumed to
%be singletons). If size(x,k) is different from both 1 and the maximum,
%then throws an error.
%Returns N and the modified arguments
function [N,varargout]=argRepmat(k,varargin)
d=cellfun(@(x) size(x,k), varargin);
N=max(d);
if ~all(d==N | d==1)
    error('Argument dimensions need to be consistent')
end
varargout=varargin;
r=ones(1,k);
r(k)=N;
for iArg=1:length(varargin)
    if d(iArg)==1
        varargout{iArg}=repmat(varargin{iArg},r);
    end
end

        