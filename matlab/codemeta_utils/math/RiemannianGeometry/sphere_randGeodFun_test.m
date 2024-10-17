function sphere_randGeodFun_test
[xt,dxt]=sphere_randGeodFun(sphere_randn(eye(3,1)));

n=@(t) norm(xt(t));

figure(1)
plotfun(n,'angle')
title('Norm of x(t)')
figure(2)
check_der(xt,dxt,'angle')