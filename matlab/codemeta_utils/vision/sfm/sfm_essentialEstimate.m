%Estimate essential matrices from the pairwise matches using  the 8 point 
% algorithm. Also, filters matches based on estimated Essential matrices.
%function data=my_sfm_addMatchEssentialEstimated(data)
%Adds field matchEssential with the essential
%matrix and field matchFiltered with the inlier matches.

function data=sfm_essentialEstimate(data,varargin)

NIter = 500;
threshold = 1e-3;

flagRefine=false;
flagShowStats=false;
memberNameEssential='matchEssentialEstimated';

ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'showstats'
            flagShowStats=true;
        case 'niter'
            ivarargin=ivarargin+1;
            NIter=varargin{ivarargin};
        case 'threshold'
            ivarargin=ivarargin+1;
            threshold=varargin{ivarargin};
        case 'flagrefine'
            ivarargin=ivarargin+1;
            flagRefine=varargin{ivarargin};
        case 'refine'
            flagRefine=true;
        case 'membernameessential'
            ivarargin=ivarargin+1;
            memberNameEssential=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end


NMatch=length(data.match);
matchEssentialEstimated=zeros(3,3,NMatch);
match=data.match;
%matchFiltered=repmat(struct('idxImg',[],'idxMatch',[],'scores',[],'residuals',[]),1,NMatch);

for iMatch=1:NMatch
    
    if flagShowStats
        fprintf('Computing Essential matrix between images %d and %d\n',data.match(iMatch).idxImg);
    end
    
    [x1,x2]=sfm_getFeatureLocationFromMatchId(data,iMatch,'normalized');
    
    % Compute essential matrix and inliers
    [E,residuals,flagInlier] = ...
        sfm_rawEssentialRansac(x1,x2,NIter,threshold);
    if flagShowStats
        fprintf('\t Num. of inliers: %d\n',sum(flagInlier));
    end
    
    if flagRefine
        [E,residuals,flagInlier]=epipolarEstimateESampleAndOutliers(x1,x2,flagInlier,...
            'threshold',threshold,'methodE','8pt');
%         [E,residuals,flagInlier]=...
%             sfm_rawEssentialLocalInlierOptimization(E,x1,x2,flagInlier,threshold);
%         E=epipolarEstimateESampson(x1(:,flagInlier),x2(:,flagInlier),E);
    end
    
    matchEssentialEstimated(:,:,iMatch)=E;
    match(iMatch).flagInlier=flagInlier;
    match(iMatch).residuals=residuals;
  
end

data.(memberNameEssential)=matchEssentialEstimated;
data.match=match;


end
