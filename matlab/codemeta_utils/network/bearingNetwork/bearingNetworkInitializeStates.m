function t_node=bearingNetworkInitializeStates(t_node)
Titruth=t_node.Titruth;
t_node.Ti=Titruth+randn(size(Titruth));
