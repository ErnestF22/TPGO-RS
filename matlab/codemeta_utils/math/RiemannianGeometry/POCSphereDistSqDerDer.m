function POCSphereDistSqDerDer
x=sphere_randGeodFun();
d=@(t) sphere_dist(sphere_eye(),x(t));
t=linspace(0,2*pi);
plot(t,funApproxDerDer(@(t) d(t)^2,t));
grid on