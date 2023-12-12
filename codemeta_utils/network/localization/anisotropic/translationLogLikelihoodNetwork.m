%Cost for localization of translation using Gaussian
%function l=translationLogLikelihoodNetwork(Ri,Ti,Tij,Gammaij,E)
function l=translationLogLikelihoodNetwork(Ri,Ti,Tij,Gammaij,E)
NEdges=size(Tij,2);
l=0;
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    
    R1=Ri(:,:,iNode);
    T1=Ti(:,iNode);
    T2=Ti(:,jNode);
    T12=Tij(:,iEdge);
    Gamma12=Gammaij(:,:,iEdge);
    
    l=l+translationLogLikelihood(R1,T1,T2,T12,Gamma12);
end
