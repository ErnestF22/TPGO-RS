close all
clear variables
randn('state',0)
n=7;
p=3;
e=stiefel_eye(zeros(n,p));
h=stiefel_tangentProj(e,randn(7,3));
h=h/sqrt(stiefel_metric(e,h,h));

m=logspace(-16,0,50);

for(im=1:length(m))
    disp(im)
    h1=h*m(im);
    y=stiefel_exp(e,h1);
    h1est=stiefel_log(e,y);
    errest=h1-h1est;
    
    [h1ref,allerr{im}]=stiefel_log_refine(e,y);
    errref=h1-h1ref;
    
    
    relerrest(im)=sqrt(sum(errest(:).^2)/sum(h1(:).^2));
    relerrref(im)=sqrt(sum(errref(:).^2)/sum(h1(:).^2));
end

loglog(m,relerrest,m,relerrref)
legend('non-refined','refined')

    