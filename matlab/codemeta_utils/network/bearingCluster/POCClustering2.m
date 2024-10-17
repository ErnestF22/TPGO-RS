function POCClustering2
%%
%resetRands(3)
%sigmaNoise=0.1; tol=0.1;
sigmaNoise=0.0; tol=1e-12;
%methodMeasurements='incidence';
methodMeasurements='cycleBasis';

d=2;
%datasetName='complex-loop-5-2';
%datasetName='loop-4';
%datasetName='trapezoid';
%datasetName='degenerate-square';
%datasetName='butterfly';

%[A,x]=bearingCluster_generateTest(datasetName,'d',d);
%[B,E]=adj2incmatrix(A,'oriented');

load('testData/testData6.mat','x','E')

[u,lambda]=bearingCluster_getBearingsScalesFromE(x,E,'noisy',sigmaNoise);

switch lower(methodMeasurements)
    case 'incidence'
        BReduced=bearingCluster_reducedMatrix(B,1);
        M=bearingCluster_measurementMatrixFromB(BReduced,u);
    case 'cyclebasis'
        CUnsigned=grCycleBasisMultiComp(E);
        C=grOrientCycleBasis(CUnsigned,E)';
        M=bearingCluster_measurementMatrixFromC(C,u);
end
[L,s]=bearingCluster_nullSpaceBasis(M,tol);

N=bearingCluster_clusterVectors(L);


[membership,info]=projective_quickshift(N,'threshold',0.5);

c=unique(membership);
cNext=max(c)+1;
Nc=length(c);
for ic=1:Nc
    idxEc=find(membership==c(ic));
    Ec=E(idxEc,:);
    membershipEc=grCompEdges(Ec);
    NcEc=max(membershipEc);
    for icEc=2:NcEc
        %relabel the different connected components in the cluster
        membership(idxEc(membershipEc==icEc))=cNext;
        cNext=cNext+1;
    end
end


fprintf('Framework contains:\n')
fprintf('\t%d shakes\n',size(L,2));
fprintf('\t%d rigid components\n',length(unique(membership)))
bearingClustering_plot(x,E,membership)


%%

save([mfilename '_data'])

