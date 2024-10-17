function output=bearingNetworkPreconditioner_test_singleTrial(NNodes,TDelay)

t_node=bearingNetworkBuildTestNetwork(NNodes);
t_node.Ti=5*randn(size(t_node.Ti));

funsBearings=bearingCostFunctions('cosine');

precondTypes={'eye','hop0','hop1n','hop1np1','idealpinv'};
NTypes=length(precondTypes);
phi=struct();
t=struct();
for iType=1:NTypes
    type=precondTypes{iType};
    
    [T.(type),H]=bearingNetworkPreconditioner(t_node,funsBearings,type);
end

neumannOrders=[1 2 10 50];
for iOrder=1:length(neumannOrders)
    order=neumannOrders(iOrder);
    type=['neumann' num2str(order)];
    T.(type)=bearingNetworkPreconditioner(t_node,funsBearings,'neumann',order);
end

neumannOrders=[1 2 10 50];
for iOrder=1:length(neumannOrders)
    order=neumannOrders(iOrder);
    type=['neumannOpt' num2str(order)];
    T.(type)=bearingNetworkPreconditioner(t_node,funsBearings,'neumannOpt',order);
end

types=fields(T);
NTypes=length(types);
for iType=1:NTypes
    type=types{iType};

    disp(type)
    optsControl={'preconditioner',T.(type),'preconditionerDelay',TDelay};
    [t.(type),x.(type)]=bearingNetworkEvolve(t_node,'tFinal',50,...
        'funsbearings',funsBearings,'optsControl',optsControl,...
        'odeSolver',@odeEuler,'optsSolver',{'MaxStep',.1,'OutputFcn',@odeWaitBar},...
        'getCost');
    outputEvolve=bearingNetworkEvolveStats(t_node,t.(type),x.(type),'cost','funsbearings',funsBearings);


    phi.(type)=outputEvolve.phi;
end

output.t=t;
output.x=x;
output.phi=phi;
output.T=T;
output.xInit=t_node.Ti;
