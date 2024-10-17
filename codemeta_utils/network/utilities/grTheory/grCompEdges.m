%Partition a list of edges into connected components
%function ncE=grCompEdges(E)
%First the graph corresponding to E is partitioned into connected
%components. Then, the label for each edge is assigned based on the label
%of its first endpoint (which should be equal to the second too).
%The labels for the components are guaranteed to be from 1 to the number of
%components
function ncE=grCompEdges(E)
ncV=grComp(E);
ncE=mapValues(ncV(E(:,1)));
