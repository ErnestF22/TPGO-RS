close all
clear variables
randn('state',0)
n=7;
p=3;
e=stiefel_eye(zeros(n,p));
h=grassman_tangentProj(e,randn(7,3));
h=h/sqrt(grassman_metric(e,h,h));

m=logspace(-16,0,50);

for(im=1:length(m))
    disp(im)
    h1=h*m(im);
    y=grassman_exp(e,h1);
    h1est=grassman_log(e,y);
    
    errest=h1-h1est;

    relerrest(im)=sqrt(sum(errest(:).^2)/sum(h1(:).^2));
end

loglog(m,relerrest)

