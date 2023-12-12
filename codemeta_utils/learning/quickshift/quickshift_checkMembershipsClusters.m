%Function to check self-consistency properties in the the output of quickshift_matching
function quickshift_checkMembershipsClusters(membershipCluster,membershipPrior,clusterIndicators)
%check that clusters with non-zero indicators are the same as
%the unique clusters referenced by the membership
clusters=unique(membershipCluster);
if exist('clusterIndicators','var')
    fprintf('Indeces of referred clusters...')
    clusterDiff=setdiff(clusters,find(sum(clusterIndicators,2)));
    if isempty(clusterDiff)
        disp('OK')
    else
        disp('Inconsistent')
    end
end

%check that each cluster does not contain points from the same prior
fprintf('Check that clusters contain points from different priors...')
flagPass=true;
for kCluster=clusters
    flagPoints=membershipCluster==kCluster;
    pointsPrior=membershipPrior(flagPoints);
    flagPass=and(flagPass,length(unique(pointsPrior))==length(pointsPrior));
end

if flagPass
    disp('OK')
else
    disp('Conflicts')
end
