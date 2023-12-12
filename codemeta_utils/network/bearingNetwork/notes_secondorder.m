function notes_secondorder
resetRands(1)

TLong=500;
TShort=3;
k=[2 2 2];
kZero=[0 0 0];
kNoFeatures=[k(1:2) 0];
allConditions={{'b',k,TLong},{'bd',k,TLong},...
    {'b',kNoFeatures,TShort},{'bd',kNoFeatures,TLong},...
    {'b',kZero,TShort},{'bd',kZero,TLong}};
lambda=0;
d=2;
t_node=bearingNetworkBuildTestNetwork(7);
%t_node.dTi=0.01*randn(size(t_node.Ti));

t_node.Ti=5*randn(size(t_node.Titruth));

t_cost.funsBearings=bearingCostFunctions('angleSq');
t_cost.funsRanges=bearingCostFunctions('squared');

t_cost.flagDynamic=true;
t_cost.alpha=[1 1];

t_cost0=t_cost;
t_node0=t_node;
NConditions=length(allConditions);
for iCondition=1:NConditions
    
    condition=allConditions{iCondition};
    baseFileName=fileNameClean([mfilename '_D' num2str(d) '_L' num2str(lambda)...
        '_' condition{1} '_K' num2str(condition{2},'%d')]);

    disp(baseFileName)

    t_cost=t_cost0;
    t_node=t_node0;
    switch condition{1}
        case 'b'
            t_cost.flagUseRanges=false;
        case 'bd'
            t_cost.flagUseRanges=true;
        otherwise
            error('Test condition for cost not recognized');
    end

    t_cost.flagUseFeatures=condition{2}(3)>0;
    optsControl={};
    optsControlDynamic={'alpha',condition{2},...
        'optsFeatures',{'remove','threshold',3},...
        'optsBearings',{'remove','threshold',3}};
    %argSolver={'odeSolver',@odeEuler,'optsSolver',{'MaxStep',.01}};
    argSolver={};

    
    figure(1)
    [t,x,t_node]=bearingNetworkEvolve(t_node,'tFinal',condition{3},'t_cost',t_cost,...
        'lambda',lambda,...
        'optsControl',optsControl,argSolver{:},...
        'optsControlDynamic',optsControlDynamic,...
        'showOdeProgress');

    output=bearingNetworkEvolveStats(t_node,t,x,'t_cost',t_cost,'cost','angles','bearingRanges','rangedifference');

    figure(2)
    bearingNetworkPlot(t_node,'flagPlotFeatureNodes',t_cost.flagUseFeatures)
    hold on
    bearingNetworkPlotTrajectories(x)
    hold off

%     figure(3)
%     semilogy(t,output.phi)
%     title('Cost')

    % figure(4)
    % plot(t,output.m)
    % title('Centroid')
    % 
    % figure(5)
    % plot(t,output.a)
    % title('Bearing angle distances')

%     figure(6)
%     plot(t,output.r)
%     title('Relative ranges corresponding to bearings')
% 
%     figure(7)
%     plot(t,[squeeze(x(3,:,:)); squeeze(x(4,:,:))],'-')
%     title('Velocities')

    save([baseFileName])
end