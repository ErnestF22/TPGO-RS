function [w,v,n]=homographyContinuousEstimateRefine(H,E,w0,v0,n0,varargin)
maxIt=100;
flagCollectIterations=false;
tolDist=1e-6;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'maxit'
            ivarargin=ivarargin+1;
            maxIt=varargin{ivarargin};
        case 'toldist'
            ivarargin=ivarargin+1;
            tolDist=varargin{ivarargin};
        case 'collect'
            flagCollectIterations=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NFrames=size(w0,2);
NPlanes=size(n0,2);

p=1;
if flagCollectIterations
    w=zeros([size(w0) maxIt+1]);
    v=zeros([size(v0) maxIt+1]);
    n=zeros([size(n0) maxIt+1]);
end
w(:,:,1)=w0;
v(:,:,1)=v0;
n(:,:,1)=n0;

%Prepare stucture with fixed data for estimation over frames
infoFrame=repmat(struct('H',[],'idxPlane',[],'NPlanes',[]),NFrames,1);
for iFrame=1:NFrames
    flagFrame=E(:,1)==iFrame;
    infoFrame(iFrame).NPlanes=sum(flagFrame);
    infoFrame(iFrame).H=H(:,:,flagFrame);
    infoFrame(iFrame).idxPlane=E(flagFrame,2);
end

%Prepare stucture with fixed data for estimation over planes
infoPlane=repmat(struct('H',[],'idxFrame',[]),NPlanes,1);
for iPlane=1:NPlanes
    flagPlane=E(:,2)==iPlane;
    infoPlane(iPlane).H=H(:,:,flagPlane);
    infoPlane(iPlane).idxFrame=E(flagPlane,1);
end

wNext=zeros(size(w0));
vNext=zeros(size(v0));
nNext=zeros(size(n0));

for it=1:maxIt
    wPrev=w(:,:,p);
    vPrev=v(:,:,p);
    nPrev=n(:,:,p);
    
    %estimate rotational velocities
    for iFrame=1:NFrames
        wHat=mean(infoFrame(iFrame).H-reshape(vPrev(:,iFrame)*vec(nPrev(:,infoFrame(iFrame).idxPlane))',3,3,[]),3);
        wNext(:,iFrame)=-vee3(wHat);
    end

    %estimate translational velocities
    for iFrame=1:NFrames
        A_f=infoFrame(iFrame).H+hat3(repmat(wNext(:,iFrame),1,infoFrame(iFrame).NPlanes));
        b_f=nPrev(:,infoFrame(iFrame).idxPlane);
        vNext(:,iFrame)=lowRankLS(A_f,b_f);
    end

    %estimate scaled plane normals
    for iPlane=1:NPlanes
        A_i=multitransp(infoPlane(iPlane).H+hat3(wNext(:,infoPlane(iPlane).idxFrame)));
        b_i=vNext(:,infoPlane(iPlane).idxFrame);
        nNext(:,iPlane)=lowRankLS(A_i,b_i);
    end
    
    if flagCollectIterations
        p=p+1;
    end
    w(:,:,p)=wNext;
    v(:,:,p)=vNext;
    n(:,:,p)=nNext;
    
    if normL2Inf(wPrev-wNext)<tolDist ...
        && normL2Inf(vPrev-vNext)<tolDist ...
        && normL2Inf(nPrev-nNext)<tolDist
        break
    end
end

if flagCollectIterations
    w=w(:,:,1:p);
    v=v(:,:,1:p);
    n=n(:,:,1:p);
end

function n=normL2Inf(A)
n=max(sqrt(sum(A.^2)));

%solves min_a \sum_i ||A_i-ab_i'||^2
function a=lowRankLS(A,b)
b=permute(b,[1 3 2]);
bNormSq=sum(b.^2);
%vectorized version of a=mean({A_i*b_i/(b_i'*b_i)}_i)
a=mean(multidiv(multiprod(A,b),bNormSq),3);
