%Generic RANSAC function
%function [model,output]=ransac(data,funModelEstimation,modelOrder,funResiduals,threshold,varargin)
%Inputs
%   data                Original input data for the estimation
%   funModelEstimation  Handle to function for estimating a model. It
%       should have prototype
%           model=funModelEstimation(subData);
%   modelOrder          Order of the model (number of datapoints needed)
%   funResiduals        Function to compute residuals give a model. It
%       should have prototype
%           e=funResiduals(data,model)
%   threshold           Threshold for deciding outliers
%Optional inputs
%   'NTrials',NTrials   Number of trials (default 100)
%   'funSelection',f    Handle to function for selecting a subset of data
%       from data and computing the data dimension. It should have prototype
%           subData=f(data,idx)
%       When idx is omitted, it should return the number of datapoints
%   'collectResiduals'  Add a field allResiduals to the output structure
%       with all the residuals from all the trials
%   'inlierEstimate'    Estimate the final model from all the inliers
%Outputs
%   model   The estimated model
%   output  Struct containing the following fields
%       flagInlier      flag vector indicating the inliers
%       NInlier         number of inliers (i.e., sum(output.flagInlier))
%       residuals       vector of residuals for the selected model
function [model,output]=ransac(data,funModelEstimation,modelOrder,funResiduals,threshold,varargin)
NTrials=100;
funSelection=@funSelectionColumns;
flagCollectResiduals=false;
flagInlierEstimate=false;
flagWaitBar=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'ntrials'
            ivarargin=ivarargin+1;
            NTrials=varargin{ivarargin};
        case 'collectresiduals'
            flagCollectResiduals=true;
        case 'inlierestimate'
            flagInlierEstimate=true;
        case 'waitbar'
            flagWaitBar=true;
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NData=funSelection(data);

bestModel=[];
bestNInliers=-1;
bestFlagInliers=zeros(1,NData);
bestResiduals=[];
bestTrial=0;
if flagCollectResiduals
    allResiduals=zeros(NTrials,NData);
end
if flagWaitBar
    w=getTextWaitBar(NTrials);
    w(0);
end
for it=1:NTrials
    idx=randperm(NData,modelOrder);
    subData=funSelection(data,idx);
    model=funModelEstimation(subData);
    r=abs(funResiduals(model,data));
    
    flagInliers=r<threshold;
    NInliers=sum(flagInliers);
    if NInliers>bestNInliers
        bestNInliers=NInliers;
        bestFlagInliers=flagInliers;
        bestModel=model;
        bestResiduals=r;
        bestTrial=it;
    end
    if flagCollectResiduals
        allResiduals(it,:)=r;
    end
    if flagWaitBar
        w(it)
    end
end

if flagInlierEstimate
    subData=funSelection(data,bestFlagInliers);
    model=funModelEstimation(subData);
else
    model=bestModel;
end
if nargout>1
    output.NInliers=bestNInliers;
    output.flagInliers=bestFlagInliers;
    output.residuals=bestResiduals;
    output.model=bestModel;
    output.threshold=threshold;
    output.trial=bestTrial;
    if flagCollectResiduals
        output.allResiduals=allResiduals;
    end
end

function subData=funSelectionColumns(data,idx)
if nargin==1
    subData=size(data,2);
else
    subData=data(:,idx,:);
end
