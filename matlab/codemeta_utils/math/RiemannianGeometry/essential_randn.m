function Q=essential_randn(Q,sigma,N)
if ~exist('Q','var') || isempty(Q)
    Q=essential_eye();
end
if ~exist('sigma','var') || isempty(sigma)
    sigma=100;
end
if ~exist('N','var') || isempty(N)
    N=size(Q,3);
end

if N~=size(Q,3)
    Q=repmat(Q(:,:,1),[1 1 N]);
end

if size(sigma,3)==1 && N>1
    sigma=repmat(sigma,[1 1 N]);
end
if size(sigma,3)>1 && N==1
    Q=repmat(Q,[1 1 N]);
    N=size(sigma,3);
end

SigmaHalf=sigmaToTransformationMatrix(sigma,6);

for iN=1:N
    Q1=Q(:,:,iN);
    vQ=rand(1)*essential_vee(Q1,essential_randTangentNormVector(Q1));
    vQ=essential_hat(Q1,SigmaHalf(:,:,iN)*vQ);
    vQ=essential_tangentProj(Q1,vQ);
    Q(:,:,iN)=essential_exp(Q1,vQ);
end
