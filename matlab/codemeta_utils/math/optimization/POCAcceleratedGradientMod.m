function POCAcceleratedGradientMod
f = @(x) x.^4;
gradient=@(x) 4 * x.^3;
x0=1;
alpha=0.2;
maxIt=100;
xs=accelerated_gradient(gradient,x0,alpha,maxIt);

y0=x0;
ysB=[y0 zeros(1,maxIt)];
xsB=[y0 zeros(1,maxIt)];
epsilons=zeros(1,maxIt);
for k=1:(maxIt+1)
    y=ysB(k);
    x=xsB(k);
    g=gradient(y);
    epsilon=(k-1)/(k+2);
    epsilons(k)=epsilon;
    xy=[x;y];
    A=[0 1;-epsilon 1+epsilon];
    B=[-alpha; -(1+epsilon)*alpha];
    xyplus=A*xy+B*g;
%     x_plus=y-alpha*g;
%     y_plus=-epsilon*x+(1+epsilon)*y-(1+epsilon)*alpha*g;
%     x_plus=0*xy(1)+xy(2)-alpha*g;
%     y_plus=-epsilon*xy(1)+(1+epsilon)*xy(2)-(1+epsilon)*alpha*g;
    x_plus=xyplus(1);
    y_plus=xyplus(2);
    ysB(k+1)=y_plus;
    xsB(k+1)=x_plus;
end

% plotfun(f,linspace(-x0,x0));
% hold on
% plot(xsB,f(xsB),'r-o')
% hold off

semilogy(f(xsB))
disp(f(xsB(end)))
% disp(max(abs(xs-xsB)))
% disp(epsilons)