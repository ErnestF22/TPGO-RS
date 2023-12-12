function POCMultilabelClustering
%datasetName='loop-rigid-4';
datasetName='butterfly';
d=3;
[A,x]=bearingCluster_generateTest(datasetName,'d',d);
[B,E]=adj2incmatrix(A,'oriented');
[u,lambda]=bearingCluster_getBearingsScalesFromE(x,E);
BReduced=bearingCluster_reducedMatrix(B,1);

membership={};
legendText={};

legendText{end+1}='Bearing';
membership{end+1}=cluster(BReduced,u,'bearing');
legendText{end+1}='Incidence, d=2';
membership{end+1}=cluster(B,u,'incidence',2);
legendText{end+1}='Incidence, d=3';
membership{end+1}=cluster(B,u,'incidence',3);

NTests=length(membership);
figure(1)
for iTest=1:NTests
    subplot(1,NTests,iTest)
    bearingClustering_plot(x,E,membership{iTest})
    title(legendText{iTest})
end


function membership=cluster(BReduced,u,method,varargin)
switch lower(method)
    case 'bearing'
        M=bearingCluster_measurementMatrixFromB(BReduced,u);
    case 'incidence'        
        d=2;
        if ~isempty(varargin)
            d=varargin{1};
        end
        M=kron(BReduced,eye(d));
end
[L,s]=bearingCluster_nullSpaceBasis(M,1e-12);
N=bearingCluster_clusterVectors(L);
[membership,info]=projective_quickshift(N,'threshold',0.5);
        