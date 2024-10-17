%function [ft,t]=plotfun(f,t,style,varargin)
%Plot the function f (which must be a function handle) by evaluating it at
%each point in the vector t and using the given line style. If t is the
%string 'angle', then t=linspace(-pi,pi,200). If length(f(t))>1, plots each
%entry as a separate graph.
%Optional argument
%   'plotDim'   1   default behaviour
%               2   use the first entry of f(t) for the x axis and the
%                   second for the y axis (i.e., plot f as a 2d curve)
%               3   same as 2 but plot f as a 3d curve
%Optional output
%   ft          a [length(f(t)) lenght(t)] matrix with f evaluated at each t
%   t           the points at which f has been evaluated
%
%See also evalfun
function varargout=funPlot(fun,t,style,varargin)
plotDim=1;
plotType='linear';
%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'plotdim'
            ivarargin=ivarargin+1;
            plotDim=varargin{ivarargin};
        case 'logy'
            plotType=lower(varargin{ivarargin});
    end
    ivarargin=ivarargin+1;
end

if ~exist('t','var') || isempty(t)
    t=[];
end
if ~exist('style','var') || isempty(style)
    style='-';
end
if ~iscell(style)
    style={style};
end

[r,t]=evalfun(fun,t);

if plotDim>1 && size(r,1)~=plotDim
    error('plotDim>1 but dimension of the function vector does not match!')
end

switch plotDim
    case 1
        switch plotType
            case 'linear'
                plot(t,r,style{:})
            case 'logy'
                semilogy(t,r,style{:})
        end
    case 2
        plot(r(1,:),r(2,:),style{:})
    case 3
        plot3(r(1,:),r(2,:),r(3,:),style{:})
end

if nargout>0
    varargout{1}=r;
    varargout{2}=t;
end