function POCquickshift_distancesIntraCluster
global NTrials
resetRands();

allN=[10 100 1e3];
allK=[3 10 100];
NTrials=100;

NN=length(allN);
NK=length(allK);
for iN=1:NN
    for iK=1:NK
        subplot(NN,NK,sub2ind([NN NK],iN,iK))
        generatePlot(allK(iK),allN(iN))
        title(sprintf('K=%d N=%d',allK(iK),allN(iN)))
    end
end




function generatePlot(K,N)
global NTrials
phi=@(x) exp(-x.^2/2);
allDist=zeros(NTrials,N);
for iTrial=1:NTrials
    X=randn(K,N);
    D=sqrt(euclideanDistMatrix(X,X));
    P=quickshift_density(phi,D);
    [vDist,vMember]=quickshift_tree(P,D);
    allDist(iTrial,:)=vDist;
end
cumDistPerc(allDist(:))
