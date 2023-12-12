%Break a tree based on prior membership information
%function [membershipCluster,treeEdges,clustersIndicators]=quickshift_breakTreePrior(treeDistances,treeEdges,membershipPrior)
function [membershipCluster,treeEdges,clustersIndicators]=quickshift_breakTreePrior(treeDistances,treeEdges,membershipPrior)
NPoints=length(treeDistances);

%vector with cluster number for each point
membershipCluster=1:NPoints;
%matrix with indicator of memberships for each cluster
clustersIndicators=sparse(1:NPoints,membershipPrior,ones(1,NPoints));

idxNotRoot=find(not(treeEdges==1:NPoints));
[~,idxIdxSortDist]=sort(treeDistances(idxNotRoot));
idxSortDist=idxNotRoot(idxIdxSortDist);

for iIdx=1:length(idxSortDist)
    %index of the candidate edge
    idx=idxSortDist(iIdx);
    
    %indeces of the two clusters that should be merged
    k1=membershipCluster(idx);
    k2=membershipCluster(treeEdges(idx));
    
    %indicator vector for the two clusters after merging
    clusterIndicatorMerged=clustersIndicators(k1,:)+clustersIndicators(k2,:);
    
    if any(clusterIndicatorMerged>1)
        %edge would create conflict, remove it
        treeEdges(idx)=idx;
    else
        %merge the clusters:
        %update cluster k1
        clustersIndicators(k1,:)=clusterIndicatorMerged;
        
        %clear cluster k2
        clustersIndicators(k2,find(clustersIndicators(k2,:)))=0;
        
        %move points from k2 to k1
        membershipCluster(membershipCluster==k2)=k1;
    end
        
end

