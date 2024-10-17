%Extract essential matrix from image points using RANSAC
%function [E,residuals,flagInlier,allResiduals,allESamples,allE] ...
%       = sfm_rawEssentialRansac(x1,x2,NSamples,threshold,varargin)
%Inputs
%   x1,x2       [2 x Nx] arrays of image points
%   NSamples    Number of samples of points to draw
%   threshold   Threshold for inliners (using SampsonAbs residuals)
%Optional Inputs
%   'progressBar'       Show text progressbar
%   'keepBestSample'    If the method used for estimating E from a sample
%                       of points produces multiple estimates, keep only
%                       the best one in allESamples (see below)
%   'methodE',m         Method to use for the estimation
%                       (see epipolarEstimateESampleAndOutliers)
%   'NRetriesMax',n     Number of data point samples to try to get a least
%                       one essential matrix sample. If reached, give up
%                       and skip to the next essential matrix sample.
%Outputs
%   E               The essential matrix selected by RANSAC
%   residuals       All the residuals for all points for the selected solution
%   flagInlier      [1 x Nx] bool mask indicating which points are inliers
%   allResiduals    All the residuals from all points and from all the samples.
%                   Useful to plot for an initial choice of the threshold
%   allESamples     All the essential matrices obtained from the point samples
%   allE            Array of best essential matrices after each sampling
%   
function [E,residuals,flagInlier,allResiduals,allESamples,allE] = sfm_rawEssentialRansac(x1,x2,NSamples,threshold,varargin)
%default values
flagGetAllResiduals=nargout>3;
flagGetAllESamples=nargout>4;
flagGetAllE=nargout>5;
flagDebug=false;

flagProgressBar=false;
flagKeepBestSample=false;
NRetriesMax=10;

methodE='5pt';
sampleSize=[];

%optional arguments
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'progressbar'
            flagProgressBar=true;
        case 'keepbestsample'
            flagKeepBestSample=true;
        case 'debug'
            flagProgressBar=false;
            flagDebug=true;
        case 'methode'
            ivarargin=ivarargin+1;
            methodE=lower(varargin{ivarargin});
        case 'samplesize'
            ivarargin=ivarargin+1;
            sampleSize=lower(varargin{ivarargin});
        case 'nretriesmax'
            ivarargin=ivarargin+1;
            NRetriesMax=lower(varargin{ivarargin});
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end


%fill in empty settings and other initializations
if isempty(sampleSize)
    switch methodE
        case '5pt'
            sampleSize = 5;
        case '5pt+validation'
            sampleSize = 8;
        case '8pt'
            sampleSize = 8;
        otherwise
            error('methodE not recognized')
    end
end

NPoints=size(x1,2);
if NPoints<sampleSize
    warning('Not enough matches to compute the essential matrix')
    E=diag([1,1,0]);
    residuals=Inf(1,NPoints);
    flagInlier=false(1,NPoints);
    return
end
[E,residuals,flagInlier]=epipolarEstimateESampleAndOutliers(x1,x2,true(1,NPoints),...
    'threshold',threshold,'methodE','8pt');
NInliers=sum(flagInlier);
if flagDebug
    fprintf('(Sample: 0, Inliers: %d) ',NInliers)
end
%Note that this is the result that will be returned in case every sample fails


%prepare structures for data collection and logging
if flagGetAllResiduals
    allResiduals=cell(NSamples,1);
end
if flagGetAllESamples
    allESamples=cell(NSamples,1);
end
if flagGetAllE
    allE=zeros(3,3,NSamples+1);
end    
if flagProgressBar
    w=getTextWaitBar(NSamples);
end

%main loop for generating essential matrix samples
for i=1:NSamples
    if flagGetAllE
        allE(:,:,i)=E;
    end

    %try to get at least one essential matrix sample by trying up to
    %NRetriesMax image pair samples
    for nRetries=1:NRetriesMax
        indices = randperm(NPoints);
        sampleInd = indices(1:sampleSize);
        flagSample=false(1,NPoints);
        flagSample(sampleInd)=true;
        
        [ESample,curResiduals,curFlagInlier]=epipolarEstimateESampleAndOutliers(x1,x2,flagSample,...
            'threshold',threshold,'methodE',methodE);
        if size(ESample,3)>0
            break
        end
    end
    if nRetries==NRetriesMax && size(ESample,3)==0
        %maximum number of retires reached, skip to the next sample
        continue
    end
    
    curNInliers = sum(curFlagInlier,2);
    
    if flagGetAllResiduals
        allResiduals{i}=curResiduals;
    end
    if ~flagKeepBestSample && flagGetAllESamples
        allESamples{i}=ESample;
    end
    
    %if multiple samples, choose only one with the most inliers
    if size(ESample,3)>1
        [curNInliers,idxMaxCurNInliers]=max(curNInliers);
        ESample=ESample(:,:,idxMaxCurNInliers);
        curResiduals=curResiduals(idxMaxCurNInliers,:);
        curFlagInlier=curFlagInlier(idxMaxCurNInliers,:);
    end
    if flagKeepBestSample && flagGetAllESamples
        allESamples{i}=ESample;
    end    
    
    %check if sample has better inliers than the best so far
    if curNInliers > NInliers
        NInliers  = curNInliers;
        residuals = curResiduals;
        flagInlier = curFlagInlier;
        E = ESample;
        if flagDebug
            fprintf('(Sample: %d, Inliers: %d) ',i,curNInliers)
        end
    end
    if flagProgressBar
        w(i)
    end
end
if flagDebug
    fprintf('\n')
end
if flagGetAllE
    allE(:,:,NSamples+1)=E;
end

%consolidate data collection outputs
if flagGetAllResiduals
    allResiduals=cat(1,allResiduals{:});
end
if flagGetAllESamples
    allESamples=cat(3,allESamples{:});
end
