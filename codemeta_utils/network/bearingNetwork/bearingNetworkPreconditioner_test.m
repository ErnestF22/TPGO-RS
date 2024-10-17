function bearingNetworkPreconditioner_test
%resetRands()

% %make the network
% A=adjgallery(4,'full');
% t_node=testNetworkCreateStruct(A);
% xtruth=[1 0; 0 0; -1 0.1; -1 -0.1]';
% t_node=bearingNetworkAddGroundTruth(t_node,xtruth);
% 
% xInit=xtruth;
% thetaInit=pi/4;
% xInit(:,1)=[cos(thetaInit); sin(thetaInit)];
% %xInit(:,[3 4])=[2.1 -0.1; 2 0.1]'; 
% t_node.Ti=xInit;

t_node=bearingNetworkBuildTestNetwork(11);
t_node.Ti=5*randn(size(t_node.Ti));

funs=bearingCostFunctions('cosine');

precondTypes={'eye','hop1n','hop1np1'};%'hop0','idealpinv'}; %,'ideal','idealn','hop1'
NTypes=length(precondTypes);
phi=struct();
t=struct();

for iType=1:NTypes
    type=precondTypes{iType};
    
    [T.(type),H]=bearingNetworkPreconditioner(t_node,funs,type);
end

neumannOrders=[1 2];% 10 50 5000];
for iOrder=1:length(neumannOrders)
    order=neumannOrders(iOrder);
    type=['neumann' num2str(order)];
    T.(type)=bearingNetworkPreconditioner(t_node,funs,'neumann',order);
end

neumannOrders=[2];% 10 50 5000];
for iOrder=1:length(neumannOrders)
    order=neumannOrders(iOrder);
    type=['neumannOpt' num2str(order)];
    T.(type)=bearingNetworkPreconditioner(t_node,funs,'neumannOpt',order);
end

types=fields(T);
NTypes=length(types);
for iType=1:NTypes
    type=types{iType};

    figure(1)
    disp(type)
    disp(sort(abs(eig(T.(type)*H)))')
    optsControl={'preconditioner',T.(type),'preconditionerDelay',5};
    [t.(type),x.(type),output]=bearingNetworkEvolve(t_node,'tFinal',50,...
        'optsControl',optsControl,...
        'odeSolver',@odeEuler,'optsSolver',{'MaxStep',.1},...
        'showOdeProgress','getCost');

    phi.(type)=output.phi;
end

c=rbg(NTypes);
figure(2)
for iType=1:NTypes
    type=types{iType};
    semilogy(t.(type),phi.(type),'color',c(iType,:));
    hold on
end
hold off
legend(fields(phi),'location','southEast')

save([mfilename '_data'])
