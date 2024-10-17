%Filter matches using essential matrix
%function [e,flagInlier,threshold]=sfm_rawMatchFilterWithEssential(x1,x2,E12,threshold)
%If the threshold is omitted, it is automatically computed by taking the
%maximum between the one given by sfm_rawThresholdEstimate and 0.005
%Returns the residuals, a flag vector indicating the outliers and the
%threshold used.
function [e,flagInlier,threshold]=sfm_rawMatchFilterWithEssential(x1,x2,E12,varargin)
flagGivenThreshold=false;
methodResiduals='sampsonabs';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'threshold'
            ivarargin=ivarargin+1;
            threshold=varargin{ivarargin};
            flagGivenThreshold=true;
        case 'methodresiduals'
            ivarargin=ivarargin+1;
            methodResiduals=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

e=epipolarResiduals(E12,x1,x2,methodResiduals);

if ~flagGivenThreshold || isempty(threshold)
    threshold=max(sfm_rawThresholdEstimate(e),0.005);
end

flagInlier=e<threshold;
