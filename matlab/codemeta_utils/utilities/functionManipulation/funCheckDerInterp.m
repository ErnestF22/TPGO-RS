%Numerically check the derivative of a function defined by interpolation
%function funCheckDerInterp(tGrid,FGrid,dF,t,varargin)
function funCheckDerInterp(tGrid,FGrid,dF,t,varargin)
F=@(t) real_interp(tGrid,FGrid,t,'spline');
funCheckDer(F,dF,t,varargin{:});
