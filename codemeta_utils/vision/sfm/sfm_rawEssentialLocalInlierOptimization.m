%Refine estimate of E by iterating between linear estimation and outlier detection
%function [E,residuals,flagInlier]=sfm_rawEssentialLocalInlierOptimization(E,x1,x2,threshold)
function [E,residuals,flagInlier]=sfm_rawEssentialLocalInlierOptimization(E,x1,x2,flagInlier,threshold,varargin)
maxIt=10;
flagDisplayStats=true;

ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'flagdisplaystats'
            ivarargin=ivarargin+1;
            flagDisplayStats=varargin{ivarargin};
        case 'maxit'
            ivarargin=ivarargin+1;
            maxIt=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end


fprintfFlag(flagDisplayStats,'# inliers: %d\n',sum(flagInlier));

for it=1:maxIt
    prevFlagInlier=flagInlier;
    [E,residuals,flagInlier]=epipolarEstimateESampleAndOutliers(x1,x2,flagInlier,'threshold',threshold);
    
    fprintfFlag(flagDisplayStats,'# inliers: %d\n',sum(flagInlier));
    
    if all(prevFlagInlier==flagInlier)
        break
    end
end
