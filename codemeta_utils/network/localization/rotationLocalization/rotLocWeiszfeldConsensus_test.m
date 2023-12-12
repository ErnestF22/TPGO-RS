function rotLocWeiszfeldConsensus_test
resetRands()
t_node=testNetworkBuildTestNetwork('references');

RRel=t_node.Rijtruth;
RInit=rot_randn(t_node.Ritruth,0.01);
E=t_node.E;
R=rotLocWeiszfeldConsensus(E,RRel,'collect','showProgress',...
    'maxIt',30,'RInit',RInit,'tolUpdate',1e-6);

c=rotLocWeiszfeldCost(E,R,RRel,'showProgress');
figure(1)
semilogy(c)
disp(c(end))
save([mfilename '_data'])