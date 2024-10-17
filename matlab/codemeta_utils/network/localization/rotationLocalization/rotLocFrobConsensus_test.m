function rotLocFrobConsensus_test
resetRands()
t_node=testNetworkBuildTestNetwork('references');

RRel=t_node.Rijtruth;
E=t_node.E;

[R,Q]=rotLocFrobConsensus(E,RRel,'collect','showProgress',...
    'maxIt',1000);%,'methodFixLeader','triangular','epsilonFactor',0.005);

c=rotLocFrobCost(E,Q,RRel);
figure(1)
semilogy(c)
disp(c(end))
save([mfilename '_data'])