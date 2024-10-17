%Re-estimate essential matrices from inliers alone
%function data=sfm_essentialRefine(data,varargin)
function data=sfm_essentialRefine(data,varargin)
flagShowStats=false;
memberEssential='matchFilteredEssential';
memberMatch='matchFiltered';
methodResiduals='sampsonabs';
method='minSampson';

ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'showstats'
            flagShowStats=true;
        case 'memberessential'
            ivarargin=ivarargin+1;
            memberEssential=varargin{ivarargin};
        case 'membermatch'
            ivarargin=ivarargin+1;
            memberMatch=varargin{ivarargin};
        case 'method'
            ivarargin=ivarargin+1;
            method=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NMatch=length(data.(memberMatch));
for iMatch=1:NMatch
    
    if flagShowStats
        fprintf('Refining Essential matrix between images %d and %d\n',data.match(iMatch).idxImg);
    end
    
    [x1,x2]=sfm_getFeatureLocationFromMatchId(data,iMatch,'normalized','member',memberMatch);
    flagInliers=data.(memberMatch)(iMatch).flagInlier;
    x1=x1(:,flagInliers);
    x2=x2(:,flagInliers);
    
    E=data.(memberEssential)(:,:,iMatch);
    if flagShowStats
        fprintf('\tEstimating from %d inliers\n',sum(flagInliers));
        residuals=epipolarResiduals(E,x1,x2,methodResiduals);
        fprintf('\tMSE residuals before refininig: %.4e\n',sqrt(mean(residuals.^2)));
    end
    
    switch lower(method)
        case '8pt'
            [E,residuals]=epipolarEstimateESampleAndOutliers(x1,x2,[],'methodE','8pt');
        case 'minsampson'
            E=epipolarEstimateESampson(x1,x2,E);
            residuals=epipolarResiduals(E,x1,x2,methodResiduals);
    end
    
    if flagShowStats
        fprintf('\tMSE residuals after refininig: %.4e\n',sqrt(mean(residuals.^2)));
    end
    
    data.(memberEssential)(:,:,iMatch)=E;
end
