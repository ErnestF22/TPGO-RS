%randn('state',0)

n=7;        %dimension of data
p=3;
N=10;       % number of samples

xtruth=[eye(p);zeros(n-p,p)];         %true estimate
sigmanoise=0.1;                 %variance of added noise


Y=zeros([n p N]);
%initialization
for(iN=1:N)
    H=stiefel_tangentProj(xtruth,randn(n,p));
    H=H/sqrt(stiefel_metric(xtruth,H,H));
    Y(:,:,iN)=stiefel_exp(xtruth,sigmanoise*randn*H);
end

Ymean=stiefel_mean(Y);

display(Ymean)
