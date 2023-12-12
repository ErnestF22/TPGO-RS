function quickshift_plotMatch(X,membershipMatches,info)
plotGroups(X,membershipMatches)
hold on
%Plot the 2-D tree on the 2-D points
quickshift_plotTree(X,info.treeEdges,'color',0.8*[1 1 1])
quickshift_plotTree(X,info.treeEdgesClusters,'color',[1 0.5 0.5])
%quickshift_plotTree(X,info.treeDistances)
quickshift_plotScales(X,info.scales,'edgecolor',0.8*[1 1 1])
hold off
axis equal
