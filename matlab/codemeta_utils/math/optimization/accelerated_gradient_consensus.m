function accelerated_gradient_consensus
resetRands()
N=50;
A=adjgallery(N,'kNeigh',2);
L=adj2laplacianmatrix(A);

f=@(x) 0.5*x'*L*x;
gradf=@(x) L*x;

x0=rand(N,1);
alpha=0.12;
maxIt=1000;

%xs=accelerated_gradient(gradf,x0,alpha,20);
[xsGrad,fxsGrad]=gradient(f,gradf,x0,2*alpha,maxIt);
[xsAccGrad,fxsAccGrad]=accelerated_gradient(f,gradf,x0,alpha,maxIt);

figure(1)
subplot(2,1,1)
semilogy(abs([fxsGrad;fxsAccGrad]'))
subplot(2,1,2)
plot([ones(N,1)'*xsGrad;ones(N,1)'*xsAccGrad]')

function [xs,fxs]=gradient(f,gradf,x0,alpha,maxIt)
xs=[x0 zeros(length(x0),maxIt)];
fxs=[f(x0) zeros(1,maxIt)];

for it=1:(maxIt+1)
    x=xs(:,it);
    g=gradf(x);
    
    x_plus=x-alpha*g;
    fx_plus=f(x_plus);
    
    xs(:,it+1)=x_plus;
    fxs(it+1)=fx_plus;
end

function [xs,fxs]=accelerated_gradient(f,gradf, y0, alpha, maxIt)
d=length(y0);
ys=[y0 zeros(d,maxIt)];
xs=[y0 zeros(d,maxIt)];
fxs=[f(y0) zeros(1,maxIt)];
for it=1:(maxIt+1)
    y=ys(:,it);
    x=xs(:,it);
    g=gradf(y);
    x_plus=y-alpha*g;
    y_plus=x_plus+((it-1)/(it+2))*(x_plus-x);
    ys(:,it+1)=y_plus;
    xs(:,it+1)=x_plus;
    fxs(it+1)=f(x_plus);
end
    
