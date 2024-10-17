%Generic RANSAC function
%function [modelsBest,output]=ransacQuickShift(data,funModelEstimation,modelOrder,funResiduals,varargin)
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
%   'optsSampleModels',opts cell array with options for ransacSampleModels
%   'collectResiduals'  Add a field allResiduals to the output structure
%       with all the residuals from all the trials
%   'inlierEstimate'    Estimate the final model from all the inliers
%   'waitbar'   Display waitbar
%Outputs
%   model   The estimated model
%   output  Struct containing the following fields
%       flagInlier      flag vector indicating the inliers
%       NInlier         number of inliers (i.e., sum(output.flagInlier))
%       residuals       vector of residuals for the selected model
function [modelsBest,output]=ransacQuickShift(data,funModelEstimation,modelOrder,funResiduals,varargin)
methodNormalization='linear';
methodScale='fixed';

thresholdCorrelationDistance=0.5;
thresholdImportance=0;
scaleFixed=0.1;

funSelection=@ransacFunSelectionColumns;

optsModelsSample={};
optsResiduals={};
optsIkose={[],'factorSigmaInliers',1};

flagInlierEstimate=true;
flagUniformOutput=true;
flagWaitBar=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'optssamplemodels'
            ivarargin=ivarargin+1;
            optsModelsSample=[optsModelsSample varargin{ivarargin}];
        case 'optsresiduals'
            ivarargin=ivarargin+1;
            optsResiduals=[optsResiduals varargin{ivarargin}];
        case 'funselection'
            ivarargin=ivarargin+1;
            funSelection=varargin{ivarargin};
        case 'methodnormalization'
            ivarargin=ivarargin+1;
            methodNormalization=lower(varargin{ivarargin});
        case 'methodscale'
            ivarargin=ivarargin+1;
            methodScale=lower(varargin{ivarargin});
        case 'inlierestimate'
            flagInlierEstimate=true;
        case 'flaginlierestimate'
            ivarargin=ivarargin+1;
            flagInlierEstimate=varargin{ivarargin};
        case 'thresholdcorrelationdistance'
            ivarargin=ivarargin+1;
            thresholdCorrelationDistance=varargin{ivarargin};
        case 'thresholdimportance'
            ivarargin=ivarargin+1;
            thresholdImportance=varargin{ivarargin};
        case 'flaguniformoutput'
            ivarargin=ivarargin+1;
            flagUniformOutput=varargin{ivarargin};
        case 'scale'
            ivarargin=ivarargin+1;
            scaleFixed=varargin{ivarargin};
            methodScale='fixed';
        case 'waitbar'
            flagWaitBar=true;
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

optsModelsSample=[optsModelsSample {'funSelection' funSelection}];
models=ransacModelsSample(data,funModelEstimation,modelOrder,optsModelsSample{:});
residuals=ransacResiduals(data,models,funResiduals,optsResiduals{:});

[NSamples,NData]=size(residuals);
scales=zeros(NSamples,1);
residualsNorm=zeros(NSamples,NData);
importances=zeros(NSamples,1);
if flagWaitBar
    w=getTextWaitBar(MSamples);
    w(0);
end
for it=1:NSamples
    switch methodScale
        case 'fixed'
            scales(it)=scaleFixed;
        case 'ikose'
            scales(it)=ikose(residuals(it,:),optsIkose{:});
        otherwise
            error('Method for chosing the scale not recognized')
    end
    switch methodNormalization
        case 'gaussian'
            residualsNorm(it,:)=exp(-residuals(it,:).^2/(2*(2.5*scales(it))^2));
        case 'linear'
            residualsNorm(it,:)=max(0,2.5*scales(it)-residuals(it,:));
        otherwise
            error('Method for normalizing residuals not recognized')
    end
    importances(it)=sum(residualsNorm(it,:));
    if flagWaitBar
        w(it)
    end
end

D=pdist2(residualsNorm,residualsNorm,'correlation');
[treeVectorDistance,treeVectorMembership]=quickshift_tree(importances,D);
flagRoots=treeVectorDistance>thresholdCorrelationDistance ...
    & importances'>thresholdImportance;
modelsBest=models(flagRoots)';
idxBest=find(flagRoots);

%if requested, re-fit each estimate using all the inliers
if flagInlierEstimate
    if nargout>1
        if flagUniformOutput
            output.modelsBestCluster=cell2mat(modelsBest);
        else
            output.modelsBestCluster=modelsBest;
        end            
    end
    for iModel=1:length(modelsBest)
        idxModel=idxBest(iModel);
        flagInliers=residuals(idxModel,:)<2.5*scales(idxModel);
        modelsBest{iModel}=funModelEstimation(funSelection(data,flagInliers));
    end
end

if flagUniformOutput
    modelsBest=cell2mat(modelsBest);
end

if nargout>1
    output.idxBest=idxBest;
    output.models=models;
    output.residuals=residuals;
    output.scales=scales;
    output.residualsNorm=residualsNorm;
    output.importances=importances;
    output.treeVectorDistance=treeVectorDistance;
    output.treeVectorMembership=treeVectorMembership;
end

