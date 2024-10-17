function POCClustering
flagFixOneScale=false;
flagAnchorNode=true;

sigmaNoise=0.0;
tol=1e-12;

NSubnet=2;
A=blkdiag(adjgallery(3,'kneigh',2),adjgallery(4,'kneigh',2));
A(4,3)=1;
A(3,4)=1;
if NSubnet>2
    A=blkdiag(A,adjgallery(5,'kneigh',2));
    A(7,8)=1;
    A(8,7)=1;
end

NNodes=size(A,1);
x=gshow(A,'adj');
x=x';

idxAnchorNode=1;
idxNonAnchorNodes=setdiff(1:NNodes,idxAnchorNode);

directions={'oriented','directed'};
for iDirection=1:2
    disp(['# Direction: ' directions{iDirection}])
    [B{iDirection},E{iDirection}]=adj2incmatrix(A,directions{iDirection});
    [u{iDirection},lambda]=bearingCluster_getBearingsScalesFromB(x,B{iDirection},'noisy',sigmaNoise);
    
    if flagAnchorNode
        BReduced{iDirection}=bearingCluster_reducedMatrix(B{iDirection},idxAnchorNode);
    end
    U{iDirection}=bearingCluster_augmentedBearingMatrix(u{iDirection});
    M{iDirection}=bearingCluster_measurementMatrixFromB(BReduced{iDirection},u{iDirection});

    NEdges=size(E{iDirection},1);

    if flagFixOneScale
        idxLambdaFix=1;
        idxLambdaNotFix=[1:(idxLambdaFix-1) (idxLambdaFix+1):NEdges];
        M{iDirection}=M{iDirection}(:,idxLambdaNotFix);
    else
        idxLambdaNotFix=1:NEdges;
    end

    [L{iDirection},s{iDirection}]=bearingCluster_nullSpaceBasis(M{iDirection},tol);
    disp(' k =')
    disp(size(L{iDirection},2))

    N{iDirection}=bearingCluster_clusterVectors(L{iDirection});
    [membership,X]=bearingCluster_clusteringFromN(N{iDirection});

    disp(' [E; N; membership] = ')
    disp([E{iDirection}(idxLambdaNotFix,:)'; N{iDirection}; membership])
    
    disp('X = ')
    disp(X)
    
%     R=N{iDirection}(:,[1 4 5]);
%     %disp(R'*R)
% 
%     disp(' [R''*VCluster] = ')
%     disp(R'*N{iDirection})
end


save([mfilename '_data'])

