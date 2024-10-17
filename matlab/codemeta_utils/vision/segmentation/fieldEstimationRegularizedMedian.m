%Estimate a field by using local estimates plus a median regularizer
%function [alphaHat,alpha,lambdaSchedule,alphaLocal]=...
%   fieldEstimationRegularizedMedian(fieldData,idxList,kList,funEstimate,varargin)
%Inputs
%   fieldData   [d x NPoints] fitting data at each location
%   idxList     [NPoints x max(kList)] neighbor's indeces for each location
%   kList       number of neighbors for each location
%   funEstimate Handle to function to do the local fitting. Must have the
%               following prototypes:
%       alpha=funEstimate(localFieldData)
%       Estimates the parameters alpha from the given data by minimizing
%       some funtion g(alpha)
%       alpha=funEstimate(localFieldData,'prior',alphaHat,'lambda',lambda)
%       Estimates the parameters alpha from the given data plus a least
%       squares regularization, i.e., minimizes 
%           lambda(1)*g(alpha)+lambda(2)*norm(alpha-alphaHat)^2.
%Optional arguments
%   'lambda',lambda     [NLambdas x 3] vectors of function parameters to use in
%                       sequence
%   'lambda2',lambda2   Replaces the current lambda(:,2) with the given
%                       lambda2 (expanding the other entries if necessary)
%   'maxItAlternation',maxIt
%                       Maximum number of coordinate descent for each
%                       lambda
%Outputs
%    alphaHat       Field from the median regularization
%    alpha          Field from the fitting
%    lambdaSchedule Lambdas used at each iteration
%    alphaLocal     Initial field estimate from local fittings
function [alphaHat,alpha,lambdaSchedule,alphaInit]=...
    fieldEstimationRegularizedMedian(fieldData,idxList,kList,funEstimate,varargin)
lambda=zeros(1,3);
lambda(1)=5;         %cost parameter for data term
lambda(2)=1e-5;      %cost parameter for soft constraint between alpha and alphaHat
lambda(3)=1;         %cost parameter for L1 norm
flagCollect=false;
maxItAlternation=10;
dAlpha=8;            %dimension of the parameter vectors (is a constant)
flagGivenInit=false;

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
        case 'dimalpha'
            ivarargin=ivarargin+1;
            dAlpha=varargin{ivarargin};
        case 'alphainit'
            ivarargin=ivarargin+1;
            alphaInit=varargin{ivarargin};
            flagGivenInit=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NX=size(fieldData,2);
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

if ~flagGivenInit
    alphaInit=fieldEstimationLocal(fieldData,idxList,kList,funEstimate,'dimAlpha',dAlpha);
end
alpha(:,:,1)=alphaInit;

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
            alphaHat(id,:,p)=medianFieldRegularized(alpha(id,:,p),idxList,kList,lambdaCurrent(2:3));
        end
        pPrev=p;
        if flagCollect
            p=p+1;
        end
        %estimate alpha from alphaHat at each point
        for iX=1:NX
            alpha(:,iX,p)=funEstimate(fieldData(:,idxList(iX,1:kList(iX))),...
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
