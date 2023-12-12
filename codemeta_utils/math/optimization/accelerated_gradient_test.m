function accelerated_gradient_test
f = @(x) x.^4;
gradf=@(x) 4 * x.^3;
x0=1;
alpha=0.1;
xs=accelerated_gradient(gradf,x0,alpha,20);

plotfun(f,linspace(-x0,x0));
hold on
plot(xs,f(xs),'r-o')
hold off
