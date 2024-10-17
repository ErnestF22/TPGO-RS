function cnormalizeDiff_test
[x,dx,~,~,ddx]=real_randGeodFun(randn(3,1));

xNorm=@(t) cnormalize(x(t));
dxNorm=@(t) cnormalizeDiff(x(t),dx(t));
ddxNorm=@(t) cnormalizeDDiff(x(t),dx(t),ddx(t));

check_der(xNorm,dxNorm)
check_der(dxNorm,ddxNorm)
