function bearingNetworkControlDynamic_test
resetRands(4)

t_node=bearingNetworkBuildTestNetwork(7);
%t_node.dTi=0.01*randn(size(t_node.Ti));

t_node.Ti=5*randn(size(t_node.Titruth));

t_cost.funsBearings=bearingCostFunctions('angleSq');
t_cost.funsRanges=bearingCostFunctions('squared');
t_cost.flagUseRanges=true;

t_cost.flagDynamic=true;
t_cost.flagUseFeatures=true;

optsControl={};
optsControlDynamic={'alpha',[2 2 2],...
    'optsFeatures',{'remove','threshold',3},...
    'optsBearings',{'remove','threshold',3}};
%argSolver={'odeSolver',@odeEuler,'optsSolver',{'MaxStep',.01}};
argSolver={};

figure(1)
[t,x,t_node]=bearingNetworkEvolve(t_node,'tFinal',100,'t_cost',t_cost,...
    'lambda',0.0001,...
    'optsControl',optsControl,argSolver{:},...
    'optsControlDynamic',optsControlDynamic,...
    'showOdeProgress');

output=bearingNetworkEvolveStats(t_node,t,x,'t_cost',t_cost,'cost','angles','bearingRanges');

figure(2)
bearingNetworkPlot(t_node,'flagPlotFeatureNodes',t_cost.flagUseFeatures)
hold on
bearingNetworkPlotTrajectories(x)
hold off

figure(3)
semilogy(t,output.phi)
title('Cost')

% figure(4)
% plot(t,output.m)
% title('Centroid')
% 
% figure(5)
% plot(t,output.a)
% title('Bearing angle distances')

figure(6)
plot(t,output.r)
title('Relative ranges corresponding to bearings')

figure(7)
plot(t,[squeeze(x(3,:,:)); squeeze(x(4,:,:))],'-')
title('Velocities')

save([mfilename '_data'])

