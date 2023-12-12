function POCCycle
flagDisplay=false;

sigmaNoise=0.0;
idxAnchorNode=1;
tol=1e-12;

x=[ 0 0 0;
    2 0 0;
    0 2 0;
    1 2 2;
    2 1 3]';
[d,N]=size(x);
A=adjgallery(N,'kneigh',1);

%generate measurements
u=bearingCluster_getBearingsScalesFromA(x,A,'noisy',sigmaNoise);

%find nullspace
[B,E]=adj2incmatrix(A,'oriented');
BReduced=bearingCluster_reducedMatrix(B,idxAnchorNode);
U=bearingCluster_augmentedBearingMatrix(u);
M=bearingCluster_measurementMatrixFromB(BReduced,u);
[L,s]=bearingCluster_nullSpaceBasis(M,tol);
if isempty(L)
    error('No nullspace detected, try to increase tolerance')
end

%fix one scale
Lp=bearingCluster_nullSpaceBasis_fixSingleScale(L,1);

%cluster
N=bearingCluster_clusterVectors(Lp);
display(N)

[membership,X]=bearingCluster_clusteringFromN(N);
% Q=rot_proj(X);
% L=L*X';
display(membership)

if flagDisplay
    %display
    K=length(unique(membership));
    disp(membership)

    E=adj2edges(A,'oriented');
    NEdges=size(E,1);
    c=jet(K);
    plotPoints(x)
    hold on
    for iEdge=1:NEdges    
        plotPoints(x(:,E(iEdge,:)),'-','Color',c(membership(iEdge),:))
    end
    hold off
    grid on
    axis equal

    disp(L)
end