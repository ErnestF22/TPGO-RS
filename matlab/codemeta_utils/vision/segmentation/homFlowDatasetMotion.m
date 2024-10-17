function [G,x]=homFlowDatasetMotion(X,G0,varargin)
NCameras=3;
sigmaRotation=0.1;
sigmaTranslation=0.3;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'nnewcameras'
            ivarargin=ivarargin+1;
            NCameras=varargin{ivarargin};
        case 'sigmarotation'
            ivarargin=ivarargin+1;
            sigmaRotation=varargin{ivarargin};
        case 'sigmatranslation'
            ivarargin=ivarargin+1;
            sigmaTranslation=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

w=generateRandomWalk(sigmaRotation*cnormalize(randn(3,1)),NCameras,sigmaRotation);
v=generateRandomWalk(sigmaTranslation*cnormalize(randn(3,1)),NCameras,sigmaTranslation);

[R0,T0]=G2RT(G0);
R=cat(3,R0,rot_exp(R0,rot_hat(R0,w)));
T=[zeros(3,1) v]+T0*ones(1,NCameras+1);
G=RT2G(R,T);

x=projectFromG(G,X,'references');

function v=generateRandomWalk(v0,N,sigma)
v=zeros(3,N);
v(:,1)=v0;
for iN=2:N
    v(:,iN)=v(:,iN-1)+sigma*max(randn(3,1),2*sigma);
end
