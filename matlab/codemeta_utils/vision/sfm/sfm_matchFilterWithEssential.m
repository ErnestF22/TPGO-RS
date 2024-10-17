%Filter matches using essential matrices
%function data=sfm_matchFilterWithEssential(data,varargin)
%Requires fields feature().locationNormalized, match, matchEssentialTruth.
%Adds field matchFiltered and match.residuals, match.flagInlier
%Optional inputs
%   'threshold',th                      threshold for filtering
%   'memberNameEssential',memberName    member containing essential
%       matrices (default: matchEssentialTruth)
function data=sfm_matchFilterWithEssential(data,varargin)
threshold=[];
flagUseFeaturesNumber=false;
thresholdFeaturesNumber=30;
matchMemberName='match';
flagShowStats=false;

memberNameEssential='matchEssentialTruth';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'showstats'
            flagShowStats=true;
        case 'threshold'
            ivarargin=ivarargin+1;
            threshold=varargin{ivarargin};
        case 'membernameessential'
            ivarargin=ivarargin+1;
            memberNameEssential=varargin{ivarargin};
        case 'flagusefeaturesnumber'
            ivarargin=ivarargin+1;
            flagUseFeaturesNumber=varargin{ivarargin};
        case 'thresholdfeaturesnumber'
            ivarargin=ivarargin+1;
            thresholdFeaturesNumber=varargin{ivarargin};
            flagUseFeaturesNumber=true;
        case 'membernamematch'
            ivarargin=ivarargin+1;
            matchMemberName=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

match=data.match;
NMatch=length(match);
matchFiltered=match;
matchUpdated=match;
matchEssential=data.(memberNameEssential);
flagKeepMatch=true(1,NMatch);
for iMatch=1:NMatch
    [x1,x2]=sfm_getFeatureLocationFromMatchId(data,iMatch,'normalized','member',matchMemberName);
    if flagUseFeaturesNumber && size(x1,2)<thresholdFeaturesNumber
        flagKeepMatch(iMatch)=false;
        fprintfFlag(flagShowStats,'Match id %d rejected (%d points)\n',iMatch, size(x1,2));
        matchUpdated(iMatch).flagValid=false;
        matchUpdated(iMatch).flagInlier=[];
        matchUpdated(iMatch).residuals=[];
        matchUpdated(iMatch).threshold=[];
    else
        E12=matchEssential(:,:,iMatch);
        [e,flagInlier,thresholdInlier]=sfm_rawMatchFilterWithEssential(x1,x2,E12,'threshold',threshold,...
            'methodResiduals','sampsonabs');
        matchFiltered(iMatch)=structCopy(matchFiltered(iMatch),sfm_rawFilterDataWithFlag(match(iMatch),flagInlier));
        matchFiltered(iMatch).threshold=thresholdInlier;
        if flagShowStats
            fprintf('Filtering match id %d\n',iMatch);
            fprintf('\tThreshold: %.4e\n',thresholdInlier);
            fprintf('\tNumber of inliers/total: %d/%d\n', sum(flagInlier), length(flagInlier));
        end
        if flagUseFeaturesNumber && sum(flagInlier)<thresholdFeaturesNumber
            flagKeepMatch(iMatch)=false;
            matchUpdated(iMatch).flagValid=false;
            if flagShowStats
                fprintf('\tMatch rejected\n');
            end
        else
            matchUpdated(iMatch).flagValid=true;
        end
            
        
        matchUpdated(iMatch).flagInlier=flagInlier;
        matchUpdated(iMatch).residuals=e;
        matchUpdated(iMatch).threshold=thresholdInlier;
    end
end
data.matchFilteredEssential=matchEssential(:,:,flagKeepMatch);
data.matchFiltered=matchFiltered(flagKeepMatch);
if flagShowStats
    fprintf('Matches retained/total: %d/%d\n', sum(flagKeepMatch), length(flagKeepMatch));
end
data.match=matchUpdated;
