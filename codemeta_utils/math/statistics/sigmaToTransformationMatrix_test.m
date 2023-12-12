function sigmaToTransformationMatrix_test
sigma=5;
d=4;
N=30000;
Sigma=sigma^2*eye(d);
SigmaHalf=sigmaToTransformationMatrix(sigma,d);
v=SigmaHalf*randn(4,N);
SigmaEst=zeros(d);
for n=1:N
    SigmaEst=SigmaEst+v(:,n)*v(:,n)';
end
SigmaEst=SigmaEst/(N-1);
disp('[Sigma SigmaEst]')
disp([Sigma SigmaEst])
