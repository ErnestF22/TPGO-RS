function [alphaHat,alpha,lambdaSchedule,alphaInit]=homFlowParametersFieldEstimate4ptTranslation(x,dx,w,idxList,kList,varargin)
fieldData=[x;dx];
alphaInit=fieldEstimationLocal(fieldData,idxList,kList,@funEstimate,'dimAlpha',6);
s=pickSign(alphaInit(4:6,:));
alphaInit=alphaInit.*(ones(6,1)*s);

[alphaHat,alpha,lambdaSchedule]=...
    fieldEstimationRegularizedMedian(fieldData,idxList,kList,@funEstimate,varargin{:},...
    'dimAlpha',6,'alphaInit',alphaInit);

    function alpha=funEstimate(fieldData,varargin)
    x=fieldData(1:2,:);
    dx=fieldData(3:4,:);
    [n,v]=homFlowParametersEstimate4pt(x,dx,w,varargin{:});
    alpha=[n;v];
    end
end

function s=pickSign(v)
vMedian=median(v,2);
s=sign(vMedian'*v);
end
