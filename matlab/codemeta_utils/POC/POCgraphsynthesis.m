%Test for generating NN graph with minimum degree constraints
function POCgraphsynthesis
NNodes=50;
v=randn(2,NNodes);
d=euclideanDistMatrix(v);
E=reshape(permute(cat(3,ones(NNodes,1)*(1:NNodes),(ones(NNodes,1)*(1:NNodes))'),[3 1 2]),2,[])';
E=E(E(:,1)<E(:,2),:);
NEdges=size(E,1);
w=zeros(NEdges,1);
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    w(iEdge)=d(iNode,jNode);
end

B=edges2incmatrix(E);

cvx_begin
    variable e(NEdges,1)
    minimize(w'*e)
    subject to
        e>=0
        e<=1
        B'*e>=3
cvx_end


disp(e')

plotPoints(v)
hold on
plotEdges(v,E(e>0.5,:))
hold off
