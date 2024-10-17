function rotLocRiemannianConsensus_test
funs=rot3_almostGlobal_functions('type','tron','b',3);

resetRands()
t_node=testNetworkBuildTestNetwork('references');

RRel=t_node.Rijtruth;
RInit=rot_randn(t_node.Ritruth,0.01);
E=t_node.E;
epsilon=1/(2*edges2maxDegree(E)*funs.mumax);
R=rotLocRiemannianConsensus(E,RRel,'collect','showProgress',...
    'maxIt',30,'funs',funs,'epsilon',epsilon,'RInit',RInit,'tolUpdate',1e-3);

c=rotLocRiemannianCost(E,R,RRel,funs,'showProgress');
figure(1)
semilogy(c)
disp(c(end))
save([mfilename '_data'])