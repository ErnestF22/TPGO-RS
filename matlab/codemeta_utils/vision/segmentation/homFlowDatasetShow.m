function homFlowDatasetShow(X,G,NVec,idxX)
plotGroups(X,idxX)
hold on
draw3dcameraFromG(G,'references','scale',0.5)
plotPoints(G2T(G),'k-')
hold off
axis equal
