function POCGraphDerivative
NNodes=7;
A=adjgallery(NNodes,'banded',2);
B=adj2incmatrix(A,'undirected');
NEdges=size(B,1);

dB=[];

%dB will have NEdges cols
%for each (i,j)\in E, add all edges (*,i) and (j,*)
for iEdge=1:NEdges
    flagNeigh=sum(and(B,ones(NEdges,1)*B(iEdge,:)),2);
    flagNeigh(iEdge)=0;
    NNeigh=sum(flagNeigh);
    dBi=zeros(NNeigh,NEdges);
end