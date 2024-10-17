function POCComputePreconditioner
t_node=bearingNetworkBuildTestNetwork();

funs=bearingCostFunctions('cosine');

[Tdiag,H]=bearingNetworkPreconditioner(t_node,funs,'hop0');
Tblock=bearingNetworkPreconditioner(t_node,funs,'hop1');
Tdblock=bearingNetworkPreconditioner(t_node,funs,'hop1p1');
Tideal=bearingNetworkPreconditioner(t_node,funs,'ideal');

disp(sort(eig(H)'))
disp(sort(eig(Tdiag*H)'))
disp(sort(eig(Tblock*H)'))
disp(sort(eig(Tdblock*H)'))
disp(sort(real(eig(Tideal*H)')))

