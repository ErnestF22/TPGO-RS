function POCRelativePositionFormationControl
resetRands()

NNodes=5;
d=2;
TFinal=20;
ECost=[...
     1     2
     2     1
     1     3
     3     1
     2     3
     3     2
     2     4
     4     2
     3     4
     4     3
     3     5
     5     3
     4     5
     5     4
     ];
t_node=relativePositionNetworkTestNetwork(NNodes);
TiTruth=t_node.TiTruth;
%EControl=ECost([1:10 12 14],:);
EControl=[ECost(ECost(:,1)>ECost(:,2),:); 1 3; 3 5];
TijTruth=relativePositionNetworkCompute(TiTruth,ECost);
TijControlTruth=relativePositionNetworkCompute(TiTruth,EControl);

Ti0=5*randn(d,NNodes);

t_nodeSim=t_node;
t_nodeSim.Ti=Ti0;
t_nodeSim.TijTruth=TijControlTruth;
t_nodeSim.E=EControl;

t_nodeStats=t_nodeSim;
t_nodeStats.E=ECost;
t_nodeStats.TijTruth=TijTruth;
figure(1)
[t,x,t_nodeSim,handleOde]=relativePositionNetworkEvolve(t_nodeSim,'TFinal',TFinal);

TiSim=x;
t_nodeStats.Ti=t_nodeSim.Ti;

figure(2)
bearingNetworkPlotTrajectories(TiSim);
hold on
relativePositionNetworkPlot(t_nodeStats,'optsEdgesTi',{'k','LineWidth',5},'optsEdgesTiTruth',{'b'})
relativePositionNetworkPlot(t_nodeSim,'optsEdgesTi',{'r','LineWidth',3},'flagEdgesTiTruth',false)
hold off
axis equal

output=relativePositionNetworkEvolveStats(t_nodeStats,t,x,'cost','gradient','control',handleOde);

figure(3)
semilogy(t,output.phi)
title('Cost')


figure(4)
dphi=shiftdim(sum(sum(output.gradPhi.*output.control)));
semilogy(t,-dphi)
title('Time derivative of cost')

figure(4)
funCheckDerInterpInterp(t,output.phi,dphi)


end
