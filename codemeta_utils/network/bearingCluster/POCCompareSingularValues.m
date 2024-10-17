function POCCompareSingularValues

sigmaNoise=0; tol=0.1;

%testName='complex-loop-5-2';
%testName='tree-3';
%testName='loop-5';
testName='loop-rigid-4';

[A,x]=bearingCluster_generateTest(testName);
[B,E]=adj2incmatrix(A,'oriented');
[u,lambda]=bearingCluster_getBearingsScalesFromB(x,B,'noisy',sigmaNoise);

methods={'incidence','cycleBasisNone','cycleBasisOrthogonal','cycleBasisCholesky'};
NMethods=length(methods);
for iMethod=1:NMethods
    switch lower(methods{iMethod})
        case 'incidence'
            BReduced=bearingCluster_reducedMatrix(B,1);
            M=bearingCluster_measurementMatrixFromB(BReduced,u);
        case {'cyclebasisnone','cyclebasisorthogonal','cyclebasischolesky'}
            E=adj2edges(A,'oriented');
            CUnsigned=grCycleBasis(E);
            C=grOrientCycleBasis(CUnsigned,E)';
            M=bearingCluster_measurementMatrixFromC(C,u,'methodNormalization',methods{iMethod}(11:end));
        otherwise
            error('Method name not recognized')
    end
    [L,s]=bearingCluster_nullSpaceBasis(M,tol);
    s=[s; zeros(max(0,size(M,2)-length(s)),1)];
    sMethods(:,iMethod)=s/max(s);
end

disp(sMethods)
plot(sMethods,'o-')
legend(methods)

