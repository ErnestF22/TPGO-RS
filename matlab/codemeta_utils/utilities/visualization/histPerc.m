%function [N,X]=histPerc(varagin)
%Same as hist, but the y axis represents relative percentages instead of
%counts.
%
%See also hist
%
function varargout=histPerc(varargin)

[N,X]=hist(varargin{:});

[b,k]=size(N);

N=N./(ones(b,1)*sum(N,1))*100;

if nargout==0
    bar(X,N)
    ylabel('Percentage')
end

if nargout>=1
    varargout{1}=N;
    if nargout>=2;
        varargout{2}=N;
    end
end
