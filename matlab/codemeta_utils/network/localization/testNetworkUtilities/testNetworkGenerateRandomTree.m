function ETree=testNetworkGenerateRandomTree(t_node)

E=testNetworkGetEdges(t_node);
ETree=graphGenerateRandomTree(E);
