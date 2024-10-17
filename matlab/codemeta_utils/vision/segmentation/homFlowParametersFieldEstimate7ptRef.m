function [alphaHat,alpha,lambdaSchedule,alphaLocal]=homFlowParametersFieldEstimate7ptRef(x,dx,idxList,varargin)
lambda=zeros(1,3);
lambda(1)=5;         %cost parameter for data term
lambda(2)=1e-5;      %cost parameter for soft constraint between alpha and alphaHat
lambda(3)=1;         %cost parameter for L1 norm
flagCollect=false;
maxItAlternation=10;
dAlpha=8;            %dimension of the parameter vectors (is a constant)

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'maxitalternation'
            ivarargin=ivarargin+1;
            maxItAlternation=varargin{ivarargin};
        case 'collect'
            flagCollect=true;
        case 'lambda'
            ivarargin=ivarargin+1;
            lambda=varargin{ivarargin};
        case 'lambda2'
            ivarargin=ivarargin+1;
            lambda2=varargin{ivarargin};
            lambda=assignExpand(lambda,lambda2,2);
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NX=size(x,2);
maxItLambda=size(lambda,1);

totalIt=maxItAlternation*maxItLambda;
lambdaSchedule=zeros(totalIt,3);

%preallocate alpha and alphaHat
if ~flagCollect
    alpha=zeros(dAlpha,NX);
    alphaHat=zeros(size(alpha));
else
    alpha=zeros(dAlpha,NX,totalIt+1);
    alphaHat=zeros(size(alpha)-[0 0 1]);
end

%initialize alpha from local estimates
for iX=1:NX
    alpha(:,iX,1)=homFlowParametersEstimate7pt(x(:,idxList(iX,:)),dx(:,idxList(iX,:)));
end
alphaLocal=alpha(:,:,1);

%p is a counter for tracking different slices of the output when collecting
%results 
p=1;
%pLambda is a counter for tracking the lambda schedule
pLambda=1;
for itLambda=1:maxItLambda
    lambdaCurrent=lambda(itLambda,:);
    for itAlternation=1:maxItAlternation
        %store lambda for this iteration
        lambdaSchedule(pLambda,:)=lambdaCurrent;
        pLambda=pLambda+1;
        
        %estimate alphaHat from alpha at each dimension
        for id=1:dAlpha
            alphaHat(id,:,p)=medianFieldRegularized(alpha(id,:,p),idxList,lambdaCurrent(2:3));
        end
        pPrev=p;
        if flagCollect
            p=p+1;
        end
        %estimate alpha from alphaHat at each point
        for iX=1:NX
            alpha(:,iX,p)=homFlowParametersEstimate7pt(x(:,idxList(iX,:)),dx(:,idxList(iX,:)),...
                'prior',alphaHat(:,iX,pPrev),'lambda',lambdaCurrent(1:2));
        end
    end
end

%Assigns lNew to l(:,k), repeating the other entries if necessary
function l=assignExpand(l,lNew,k)
lNew=shiftdim(lNew);
d=size(l,1);
dNew=size(lNew,1);
if dNew>1 && ~(d==1 || d==dNew)
    error('Assignment of column not possible due to size mismatch')
end
if dNew>d
    l=repmat(l,dNew,1);
end
l(:,k)=lNew;

