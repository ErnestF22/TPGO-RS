function [SEst,RbsEst,taubEst,output]=affineSfM(xb,varargin)
RFirstFrame=eye(3);
flagCollect=nargout>3;
flagShowDiagnosticInfo=false;
flagFlip=true;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'rfirstframe'
            ivarargin=ivarargin+1;
            RFirstFrame=varargin{ivarargin};
        case 'showdiagnosticinfo'
            flagShowDiagnosticInfo=true;
        case 'flagflip'
            ivarargin=ivarargin+1;
            flagFlip=varargin{ivarargin};
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

F=size(xb,3);
P=size(xb,2);
W=matStack(xb);

%2-D translation and centering of W
c=mean(W,2);
WCentered=W-c*ones(1,P);
if flagCollect
    output.W=W;
    output.c=c;
    output.WCentered=W;
end

%rank 3 factorization
[UWcentered,SWCentered,VWCentered]=svd(WCentered,'econ');
MEstAffine=UWcentered(:,1:3);
SEstAffine=SWCentered(1:3,1:3)*VWCentered(:,1:3)';
if flagShowDiagnosticInfo
    disp('Error of rank-3 affine factorization')
    disp(norm(WCentered-MEstAffine*SEstAffine,'fro'))
end
if flagCollect
    output.MEstAffine=MEstAffine;
    output.SEstAffine=SEstAffine;
end

%metric upgrade
[RReducedEst,KEst]=affineSfMMetricUpgrade(matUnstack(MEstAffine,2));
SEstMetric=KEst\SEstAffine;
if flagCollect
    output.RReducedEst=RReducedEst;
    output.KEst=KEst;
    output.SEst=SEstMetric;
end

%complete third row of rotations
REstMetric=zeros(3,3,F);
for f=1:F
    REstMetric(:,:,f)=orthCompleteBasis(RReducedEst(:,:,f)')';
end
if flagCollect
    output.REstMetric=REstMetric;
end

%flip z coordinate
if flagFlip
    KFlip=diag([1 1 -1]);
    REstMetric=multiprod(KFlip,multiprod(REstMetric,KFlip));
    SEstMetric=KFlip'*SEstMetric;
    if flagCollect
        output.KFlip=KFlip;
    end
end

%align to RFirstFrame
KAlign=rot_proj(REstMetric(:,:,1))'*RFirstFrame;
RbsEst=multiprod(REstMetric,KAlign);
SEst=KAlign\SEstMetric;
if flagCollect
    output.KAlign=KAlign;
end

%Estimated translation in body frame
%Note that the third coordinate is arbitrary
taubEst=[reshape(c,2,F); -2*ones(1,F)];
