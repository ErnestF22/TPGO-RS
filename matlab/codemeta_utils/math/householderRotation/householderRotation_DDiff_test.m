function householderRotation_DDiff_test
k=3;
[x,dx,~,~,ddx]=real_randGeodFun(randn(3,1),'speed','quadratic');
R0=@(t) householderRotation(x(t),k);
dR0Vec=@(t) householderRotation_Diff(x(t),k,dx(t));
dR0=@(t) R0(t)*hat(dR0Vec(t));
ddR0Vec=@(t) householderRotation_DDiff(x(t),k,dx(t),ddx(t));

%check_der(R0,dR0)
check_der(dR0Vec,ddR0Vec)
