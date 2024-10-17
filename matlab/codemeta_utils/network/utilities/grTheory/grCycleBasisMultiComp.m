function C=grCycleBasisMultiComp(E)
membership=grCompEdges(E);
NE=size(E,1);
NComponents=max(membership);
C=[];
for iComponent=1:NComponents
    idxComponent=find(membership==iComponent);
    EComponent=mapValues(E(idxComponent,:));
    CComponent=grCycleBasis(EComponent);
    CComponentFull=zeros(NE,size(CComponent,2));
    CComponentFull(idxComponent,:)=CComponent;
    C=[C CComponentFull];
end
