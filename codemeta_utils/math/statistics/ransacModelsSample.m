%Sample models for RANSAC
%function models=ransacModelsSample(data,funModelEstimation,modelOrder,varargin)
%Inputs
%   data                Original input data for the estimation
%   funModelEstimation  Handle to function for estimating a model. It
%       should have prototype
%           model=funModelEstimation(subData);
%   modelOrder          Order of the model (number of datapoints needed)
%   threshold           Threshold for deciding outliers
%Optional inputs
%   'NTrials',NTrials   Number of trials (default 100)
%   'funSelection',f    Handle to function for selecting a subset of data
%       from data and computing the data dimension. It should have prototype
%           subData=f(data,idx)
%       When idx is omitted, it should return the number of datapoints
%   'waitbar'   Display waitbar
%   'flagUniformOutput',flag    if true, models are concatenated in a
%       single matrix row-wise (default: false). E.g. if each model is a 
%       [d x 1] vector, then models becomes a [d x NSamples] matrix
%Outputs
%   models  [NSamples x 1] cell array of data
function models=ransacModelsSample(data,funModelEstimation,modelOrder,varargin)
NSamples=100;
funSelection=@ransacFunSelectionColumns;
flagWaitBar=false;
flagUniformOutput=false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'nsamples'
            ivarargin=ivarargin+1;
            NSamples=varargin{ivarargin};
        case 'funselection'
            ivarargin=ivarargin+1;
            funSelection=varargin{ivarargin};
        case 'waitbar'
            flagWaitBar=true;
        case 'flaguniformoutput'
            ivarargin=ivarargin+1;
            flagUniformOutput=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NData=funSelection(data);
models=cell(NSamples,1);

if flagWaitBar
    w=getTextWaitBar(NSamples);
    w(0);
end
for it=1:NSamples
    idx=randperm(NData,modelOrder);
    subData=funSelection(data,idx);
    models{it}=funModelEstimation(subData);
    if flagWaitBar
        w(it)
    end
end

if flagUniformOutput
    models=models';
    models=[models{:}];
end

