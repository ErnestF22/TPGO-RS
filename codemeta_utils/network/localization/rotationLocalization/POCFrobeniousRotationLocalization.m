function POCFrobeniousRotationLocalization
resetRands()
t_node=testNetworkBuildTestNetwork('references');

RRel=t_node.Rijtruth;
E=t_node.E;

R=rotLocFrobConsensus(E,RRel,'collect','showProgress',...
    'maxIt',1000);%,'methodFixLeader','triangular','epsilonFactor',0.005);

c=rotLocFrobCost(E,R,RRel);
figure(1)
semilogy(c)
disp(c(end))
save([mfilename '_data'])