function cdc16_collisionAvoidance
resetRands(1)
NNodes=7;
TFinal=1500;
optsControl={};
t_cost.funsRanges=bearingCostFunctions('squared');
t_node=bearingNetworkBuildTestNetwork(NNodes,2,'Edges',bnFranchi_graphEdges(NNodes));
t_node.Ti=-t_node.Titruth+0.1*randn(2,NNodes);
Eca=adj2edges(ones(NNodes));
rca=1.5;
t_node.Eca=Eca;
t_node.rca=rca;
t_cost.funsBearings=bearingCostFunctions('cosine');
t_cost.alpha=1;

%argSolver={'odeSolver',@odeEuler,'optsSolver',{'MaxStep',.01}};
%argSolver={'odeSolver',@ode15s};
argSolver={'odeSolver',@ode113};

t_node0=t_node;
t_cost0=t_cost;
allFlags={false,true};
baseName=[mfilename '_data'];
for iFlag=1:length(allFlags)
    t_node=t_node0;
    t_cost=t_cost0;
    flagCA=allFlags{iFlag};
    t_cost.flagCollisionAvoidance=flagCA;
    figure(1)
    [t,x,t_node,ode]=bearingNetworkEvolve(t_node,'tFinal',TFinal,'t_cost',t_cost,...
        'optsControl',optsControl,argSolver{:},...
        'showOdeProgress');
    bearingNetworkPlotTrajectories(x)
    hold on
    bearingNetworkPlot(t_node)
    hold off
    results{iFlag}.t_node=t_node;
    results{iFlag}.x=x;
    results{iFlag}.t=t;
    
    figure(2)
    d=computeAllDistances(x,Eca);
    plot(t,min(d,[],2))
    ax=axis;
    ax(3:4)=[0 5];
    axis(ax)
    results{iFlag}.d=d;
    save([baseName '_CA' num2str(flagCA)])
end
save(baseName)



function d=computeAllDistances(x,Eca)
NIt=size(x,3);
NEdges=size(Eca,1);
d=zeros(NIt,NEdges);
for it=1:NIt
    d(it,:)=cnorm(x(:,Eca(:,1),it)-x(:,Eca(:,2),it));
end



