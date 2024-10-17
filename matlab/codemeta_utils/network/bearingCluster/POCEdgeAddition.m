function POCEdgeAddition

%%
dataset=2;
switch dataset
    case {1,2}
        switch dataset
            case 1
                [A,x]=bearingCluster_generateTest('butterfly');
            case 2
                [A,x]=bearingCluster_generateTest('butterfly-bend');
        end
        E.original=adj2edges(A,'oriented');
        E.rigid=[E.original;1 5];
        E.flexible=[E.original; 1 7];
    case 3
        [A,x]=bearingCluster_generateTest('loop-5+complete');
        E.original=adj2edges(A,'oriented');
        E.added=[E.original; 1 3];
end

typesAddition=fields(E);

disp('# Nullity')
for iType=1:length(typesAddition)
    type=typesAddition{iType};
    
    u.(type)=bearingCluster_getBearingsScalesFromE(x,E.(type));
    C.(type)=grOrientedCycleBasis(E.(type))';
    M.(type)=bearingCluster_measurementMatrixFromC(C.(type),u.(type),'methodNormalization','none');
    L.(type)=null(M.(type));
    [Lp.(type),l.(type)]=cnormalize(L.(type)');
    Lp.(type)=Lp.(type)';
    
    NAdd=size(E.(type),1)-size(E.original,1);
    ML.(type)=M.(type)*blkdiag(L.original,diag(ones(NAdd,1)));
    
    fprintf('\t%s\n',type)
    disp([size(null(M.(type))) size(null(ML.(type)))])
end

%%
