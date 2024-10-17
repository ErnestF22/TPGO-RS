%Estimate multi-homography with an initial guess for the rotations
%function [T,NVec,RRel]=homographyEstimateIMUAssisted(x,RRel,idxX,varargin)
%This function estimates discrete homographies from multiple views and planes,
%uses an estimate of the relative rotations (RRel) to transform them into
%continuous-like homographies, and uses
%homographyContinuousEstimateNuclearNorm to compute the translation,
%rotation corrections and plane norms. The process can be iterated to use
%successively better rotation estimates.
%Inputs
%   x       [2 x NPoints x NViews] array of images
%   RRel    initial guess for the relative rotation
%   idxX    [1 x NPoints] vector indicating the images' memberships
%           If omitted or empy, all points are assumed to belong to the
%           same plane.
%Optional arguments
%   'flagCompensateRotations',f     if true, correct RRel by the rotation
%       correction estimated through the homographies (default: true)
%   'maxIt', it                     number of times to repeat the
%       estimation (default: 2). At each iteration, the rotation guess is
%       updated with the correction estimated by the homography. If
%       'flagCompensateRotation' is false, maxIt is automatically set to 1.
function [T,NVec,RRel]=homographyEstimateIMUAssisted(x,RRel,idxX,varargin)
flagCompensateRotation=true;
maxIt=2;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'flagcompensaterotation'
            ivarargin=ivarargin+1;
            flagCompensateRotation=varargin{ivarargin};
        case 'maxit'
            ivarargin=ivarargin+1;
            maxIt=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if ~exist('idxX','var') || isempty(idxX)
    idxX=ones(1,size(x,2));
end

%normalize membership indeces (no empty clusters)
idxX=mapValues(idxX);

%If compensation of the rotations is not allowed, then it does not make
%sense to iterate
if ~flagCompensateRotation
    maxIt=1;
end

NFrames=size(x,3)-1;
NPlanes=max(idxX);
NPairs=NFrames*NPlanes;
H=zeros(3,3,NPairs);
HContinuous=zeros(3,3,NPairs);
E=zeros(NPairs,2);
idxPlanes=reshape(1:NPairs,NFrames,NPlanes);
for iPlane=1:NPlanes
    xGroup=x(:,idxX==iPlane,:);
    HPlane=homographyEstimateLinear(xGroup);
    HPlane=homographyNormalize(HPlane,xGroup);
    H(:,:,idxPlanes(:,iPlane))=HPlane;
    E(idxPlanes(:,iPlane),1)=1:NFrames;
    E(idxPlanes(:,iPlane),2)=iPlane;
end

for it=1:maxIt
    for iPlane=1:NPlanes
        HPlane=H(:,:,idxPlanes(:,iPlane));
        HContinuous(:,:,idxPlanes(:,iPlane))=multiprod(multitransp(RRel),HPlane)-repmat(eye(3),1,1,size(HPlane,3));
    end
    [w,v,n]=homographyContinuousEstimateNuclearNorm(HContinuous,E);

    T=-squeeze(multiprod(RRel,permute(v,[1 3 2])));
    NVec=planeNScaledToNVec(n);

    if flagCompensateRotation
        I=rot_eye(RRel);
        RCompensate=rot_exp(I,rot_hat(I,-w));
        RRel=multiprod(RRel,RCompensate);
    end
end
