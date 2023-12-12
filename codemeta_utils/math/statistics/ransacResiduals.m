%Compute residuals from the data and models
%function residuals=ransacResiduals(data,models,funResiduals,varargin)
%Inputs
%   data                Original input data for the estimation
%   models              [1 x NModels] cell array containing the models to
%       evaluate
%   funResiduals        Function to compute residuals give a model. It
%       should have prototype
%           e=funResiduals(data,model)
%Optional inputs
%   'waitbar'   Display waitbar
%Outputs
%   residuals   [NSamples x NData] matrix with the residuals for each model
function residuals=ransacResiduals(data,models,funResiduals,varargin)
flagWaitBar=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'waitbar'
            flagWaitBar=true;
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NSamples=length(models);
%get the residuals for the first model and their length
residuals1=funResiduals(models{1},data);
NData=length(residuals1);
if flagWaitBar
    w=getTextWaitBar(NSamples);
    w(0);
end

%allocate the entire matrix for the resiudals and compute all of them
residuals=zeros(NSamples,NData);
residuals(1,:)=residuals1;
if flagWaitBar
    w(1);
end
for it=2:NSamples
    residuals(it,:)=abs(funResiduals(models{it},data));
    if flagWaitBar
        w(it)
    end
end

