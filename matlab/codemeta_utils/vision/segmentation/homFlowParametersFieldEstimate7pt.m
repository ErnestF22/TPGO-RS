function [alphaHat,alpha,lambdaSchedule,alphaLocal]=homFlowParametersFieldEstimate7pt(x,dx,idxList,kList,varargin)
fieldData=[x;dx];
[alphaHat,alpha,lambdaSchedule,alphaLocal]=...
    fieldEstimationRegularizedMedian(fieldData,idxList,kList,@funEstimate,varargin{:});

function alpha=funEstimate(fieldData,varargin)
x=fieldData(1:2,:);
dx=fieldData(3:4,:);
alpha=homFlowParametersEstimate7pt(x,dx,varargin{:});

