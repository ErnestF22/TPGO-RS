function POCLocalize_essential
resetRands();
t_node=testNetworkBuildTestNetwork('methodInit','noisytruth',...
    'varNoisyTruth',[0.5 0.5]);

optsLieMinimize={'disablePhase12','displayIt'};
totalMaxIt=1000;
E=testNetworkGetEdges(t_node);
N=testNetworkGetNumberOfNodes(t_node);

% [Ri0,Ti0]=testNetworkGetRotTransl(t_node);
% GiInit=RT2G(Ri0,Ti0);
GiInit=RT2G(rot_randn([],[],N), 8*randn(3,N));

[~,Qij]=testNetworkGetRelativeEssential(t_node);
%Qij=t_node.QEijtruth;

fEssential=@(G) essentialCostNetwork(G2R(G),G2T(G),Qij,E);

figure(1)
[GiEst,errorsGiEst]=lie_minimizeGradNewton(rot3r3_funs(),fEssential,GiInit,...
    'epsilon',1,'MaxIt',totalMaxIt,optsLieMinimize{:},'getErrors','phase3IncreaseCounterMax',30); %,'showCost'
[GiEst,errorsGiEstNoEpsilonIncrease]=lie_minimizeGradNewton(rot3r3_funs(),fEssential,GiInit,...
    'epsilon',1,'MaxIt',totalMaxIt,optsLieMinimize{:},'getErrors','phase3IncreaseCounterMax',inf); %,'showCost'
[GiEst,errorsGiEstLessEpsilonIncrease]=lie_minimizeGradNewton(rot3r3_funs(),fEssential,GiInit,...
    'epsilon',1,'MaxIt',totalMaxIt,optsLieMinimize{:},'getErrors','phase3IncreaseCounterMax',60); %,'showCost'


GiTruth=t_node.gitruth;
optsDisplay={'references','adjmatrix',t_node.A};
figure(1)
subplot(1,2,1)
testNetworkDisplay(GiInit,optsDisplay{:},'Estimated')
hold on
testNetworkDisplay(GiTruth,optsDisplay{:})
hold off
title('Init')
subplot(1,2,2)
testNetworkDisplay(GiEst,optsDisplay{:},'Estimated')
hold on
testNetworkDisplay(GiTruth,optsDisplay{:})
hold off
title('Optimized')

figure(2)
semilogy(errorsGiEst.t,errorsGiEst.cost);
hold on
semilogy(errorsGiEstNoEpsilonIncrease.t,errorsGiEstNoEpsilonIncrease.cost,'r');
semilogy(errorsGiEstLessEpsilonIncrease.t,errorsGiEstLessEpsilonIncrease.cost,'g');
hold off
legend('With epsilon increases','Without epsilon increases','Less epsilon increases')

save([mfilename '_data'])
