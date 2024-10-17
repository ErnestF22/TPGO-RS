%Numerically check the differential of a map
%function funCheckDifferential(F,x,DF,dx,t)
%Inputs
%   F   handle to map of the form F(x)
%   x   handle to curve in the domain of F
%   DF  handle to differential in the form DF(x,dx)
%   dx  handle to tangent of the curve x 

function [varargout] = funCheckDifferential(F,x,DF,dx,t, varargin)

% If it does NOT exist a variable named 't' create it 
if ~exist('t','var')
    % Create a linearly spaced vector t between -1 and 1 of size 100
    t=linspace(-1,1,100);
end

Fx=@(t) F(x(t));
dFx=@(t) DF(x(t),dx(t));

if (~isempty(varargin))
    [Ft, dFt, appder] = funCheckDer(Fx,dFx,t, varargin);
else
    [Ft, dFt, appder] = funCheckDer(Fx,dFx,t);
end

if nargout>=1
    varargout{1}=Ft;
    if nargout>=2
        varargout{2}=dFt;
        if nargout>=3
            varargout{3}=appder;
        end
    end
end