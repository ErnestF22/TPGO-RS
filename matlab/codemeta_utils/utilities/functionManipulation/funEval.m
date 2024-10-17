%function varargout=evalfun(f,t)
%Evaluate the function f (which must be a function handle) by computing it at
%each point in the vector t and using the given line style. If t is the
%string 'angle', then t=linspace(-pi,pi,200).
%Output
%   ft          a [length(f(t)) lenght(t)] matrix with f evaluated at each t
%   t           the points at which f has been evaluated
function [r,varargout]=funEval(fun,t)
if ~exist('t','var') || isempty(t)
    t=linspace(0,1,100);
end
if exist('style','var')==0 || isempty(t)
    style='-';
end

if ischar(t)
    switch lower(t)
        case 'angle'
            t=linspace(-pi,pi,201);
        case {'posangle','halfangle'}
            t=linspace(0,pi,201);
        otherwise
            error('String for t not recognized!')
    end
end

for it=1:length(t)
    r(it,:)=reshape(fun(t(it)),[],1);
end

r=r';
if nargout>1
    varargout{1}=t;
end
