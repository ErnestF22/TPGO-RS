function Sigma=randnCovariance(d,N)
Sigma=zeros([d d N]);
for iN=1:N
    preSigma=randn(d);
    Sigma(:,:,iN)=preSigma*preSigma';
end
