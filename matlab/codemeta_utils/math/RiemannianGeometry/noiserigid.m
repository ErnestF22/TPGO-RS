%Given a rigid body transformation g, return a g1 by perturbing rotation
%and translation with Gaussian noise

%%AUTORIGHTS%%

function g1=noiserigid(g,sigmaNoiseXR,sigmaNoiseT)
flagJointNoise=false;
if(~exist('sigmaNoiseT','var'))
    flagJointNoise=true;
end

N=size(g,3);

if size(sigmaNoiseXR,3)==1
    sigmaNoiseXR=repmat(sigmaNoiseXR,[1 1 N]);
end
if ~flagJointNoise
    if size(sigmaNoiseT,3)==1
        sigmaNoiseT=repmat(sigmaNoiseT,[1 1 N]);
    end
end

if ~flagJointNoise
    SigmaHalfT=sigmaToTransformationMatrix(sigmaNoiseT);
    SigmaHalfR=sigmaToTransformationMatrix(sigmaNoiseXR);
else
    SigmaHalf=sigmaToTransformationMatrix(sigmaNoiseXR,6);
end


g1=g;

for iN=1:N
    if ~flagJointNoise
        vT=SigmaHalfT(:,:,iN)*randn(3,1);
        vR=SigmaHalfR(:,:,iN)*randn(3,1);
    else
        v=SigmaHalf(:,:,iN)*randn(6,1);
        vR=v(1:3);
        vT=v(4:6);
    end
    g1(1:3,4,iN)=g(1:3,4,iN)+vT;
    g1(1:3,1:3,iN)=g(1:3,1:3,iN)*rot(vR);
end
