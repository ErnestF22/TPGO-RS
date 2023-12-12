function dF=funApproxDerInterp(tGrid,FGrid,t,varargin)
F=@(t) real_interp(tGrid,FGrid,t,'spline');
dF=funApproxDer(F,t,varargin{:});
