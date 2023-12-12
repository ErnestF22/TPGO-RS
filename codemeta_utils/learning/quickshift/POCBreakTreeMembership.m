function POCbreakTreeMembership

optsCluster={'methodBreakTree','descendents',...
    'optsBreakTree',{'ratio',1,'distancesDefault',0.1}};

[X,membershipPrior]=quickshift_test_datasets('correspondences');

%compute pairwise distances
D=sqrt(euclideanDistMatrix(X,X));

%Gaussian kernel
phi=@(x) exp(-x.^2/2);

%Do the clustering
[~,info]=quickshift_cluster(D,'phi',phi,optsCluster{:});

treeEdges=info.treeEdges;
treeDistances=info.treeDistances;

NPoints=size(X,2);

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

plotGroups(X,membershipPrior)
hold on

%Plot the 2-D tree on the 2-D points
quickshift_plotTree(X,info.treeEdges,'color',0.8*[1 1 1])
quickshift_plotTree(X,info.treeEdgesClusters,'color',[1 0 0])
quickshift_plotTree(X,treeEdges)
hold off

quickshift_checkMembershipsClusters(membershipCluster,membershipPrior,clustersIndicators)
