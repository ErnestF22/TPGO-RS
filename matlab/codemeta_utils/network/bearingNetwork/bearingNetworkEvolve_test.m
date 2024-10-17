function bearingNetworkEvolve_test
%try to find local minima for the cost with projection
resetRands()

flagPreconditioner=false;
flagLeaders=0;

%make the network
% A=adjgallery(4,'full');
% t_node=testNetworkCreateStruct(A);
% 
% xtruth=[1 0; 0 0; -1 0.1; -1 -0.1]';
% t_node=bearingNetworkAddGroundTruth(t_node,xtruth);

t_node=bearingNetworkBuildTestNetwork(7);
% xInit=xtruth;
% thetaInit=pi/4;
% xInit(:,1)=[cos(thetaInit); sin(thetaInit)];
% %xInit(:,[3 4])=[2.1 -0.1; 2 0.1]'; 
% t_node.Ti=xInit;

t_node.Ti=5*randn(size(t_node.Titruth));

t_cost.funsBearings=bearingCostFunctions('angleSq');
optsControl={};
if flagLeaders
    s=0.001;
    switch flagLeaders
        case 1
            idxLeader=4;
            dxLeader=@(t) t*s*[0;1];
            optsControl=[optsControl {'leader',idxLeader,dxLeader}];
        case 2
            idxLeader=[3 4];
            dxLeader=@(t) t*s*[0;1]*[1 1];
            optsControl=[optsControl {'leader',idxLeader,dxLeader}];
            t_node.Ti(:,idxLeader(2))=t_node.Ti(:,idxLeader(1))+(t_node.Titruth(:,idxLeader(2))-t_node.Titruth(:,idxLeader(1)));
    end
    t_node.Ti=t_node.Ti-(t_node.Ti(:,idxLeader(1))-t_node.Titruth(:,idxLeader(1)))*ones(1,t_node.NNodes);
end
%optsControl={'project',1};

if flagPreconditioner
    T=bearingNetworkPreconditioner(t_node,funsBearings,'hop1n');
    % T=bearingNetworkPreconditioner(t_node,funsBearings,'neumannOpt',50);
    optsControl=[optsControl {'preconditioner',T,'preconditionerDelay',5}];
end

figure(1)
[t,x,t_node]=bearingNetworkEvolve(t_node,'tFinal',10,'t_cost',t_cost,...
    'optsControl',optsControl,'odeSolver',@odeEuler,'optsSolver',{'MaxStep',.1},...
    'showOdeProgress');
output=bearingNetworkEvolveStats(t_node,t,x,'t_cost',t_cost,'cost','angles','residuals');

figure(2)
bearingNetworkPlot(t_node)
hold on
plotPoints(t_node.Ti0,'rx')
plot(squeeze(x(1,:,:))',squeeze(x(2,:,:))','-.','color',[1 0.75 0])
hold off

figure(3)
semilogy(t,output.phi)
title('Cost')

figure(4)
plot(t,output.m)
title('Centroid')

figure(5)
plot(t,output.a)
title('Bearing angle distances')

