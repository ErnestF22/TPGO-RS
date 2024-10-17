function G=directedNetworkStabilityMatrix(E)
%deduce number of nodes
NNodes=max(E(:));

B=edges2incmatrix(fliplr(E),NNodes,'directed');
G=(B<0)'*B;
