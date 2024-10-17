% compare two functions
function [varargout] = funCompare(fa,fb,t,varargin)
% INPUTS:
%   fa, fb := Functions of variable t
%   t [optional] := array of time to compare fa and fb
%   varargin := see below
% OUTPUTS:
%   varargout := {1} numerical evaluation of fa
%                {2} numerical evaluation of fb
%                {3} time of evaluation

flagDisplayError=true;
%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'nodisplay'
            flagDisplayError=false;
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if ~exist('t','var')
    t=linspace(-1,1,100);
end
if ischar(t)
    switch t
        case 'angle'
            t=linspace(-pi,pi,101);
        otherwise
            error('String for t not recognized')
    end
end

[faeval,ta]=evalfun(fa,t);
[fbeval,tb]=evalfun(fb,t);

plot(ta,faeval)
hold on
plot(tb,fbeval,'rx')
hold off

if flagDisplayError
    disp(['Maximum absolute error: ' num2str(max(abs(faeval-fbeval)))])
    disp(['Maximum relative error: ' num2str(max(abs((faeval-fbeval)./faeval)))])
end

if nargout>=1
    varargout{1}=faeval;
    if nargout>=2
        varargout{2}=fbeval;
        if nargout>=3
            varargout{3}=t;
        end
    end
end