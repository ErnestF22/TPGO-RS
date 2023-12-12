%function [xt,dxt,x0,dx]=real_randGeodFun(x0,varargin)
%Creates a random normal geodesic (straight line) starting from x0
function [xt,dxt,x0,dx0,ddxt]=real_randGeodFun(x0,varargin)
v=cnormalize(randn(size(x0)));
[xt,dxt,x0,dx0,ddxt]=real_geodFun(x0,v,varargin{:});

