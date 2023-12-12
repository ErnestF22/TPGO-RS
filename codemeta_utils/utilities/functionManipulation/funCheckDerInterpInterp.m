%Numerically check the derivative of a function, where both are defined by interpolation
%function funCheckDerInterpInterp(tGrid,FGrid,dFGrid,t,varargin)
function funCheckDerInterpInterp(tGrid,FGrid,dFGrid,t,varargin)
if ~exist('t','var') || isempty(t)
    t=linspace(min(tGrid),max(tGrid));
end

    
F=@(t) real_interp(tGrid,FGrid,t,'spline');
dF=@(t) real_interp(tGrid,dFGrid,t,'spline');
funCheckDer(F,dF,t,varargin{:});
