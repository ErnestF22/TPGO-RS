function POCbearingNetworkMinimize
figDir='../papers/control/notes/figures';
figDim=[300,200];
flagSaveFigure=0;

resetRands()
%costName='cosine';
costName='angleSq';

funs=bearingCostFunctions(costName);

t_node=bearingNetworkBuildTestNetwork();
E=t_node.E;
ytruth=t_node.Yijtruth;

NIt=3000;
%x0=t_node.Ti;
x0=6*randn(2,t_node.NNodes);
x0=x0-(mean(x0,2)+[0;1])*ones(1,t_node.NNodes);
x=zeros([size(x0) NIt]);
x(:,:,1)=x0;
m=zeros(2, NIt);
m(:,1)=mean(x0,2);
c=zeros(NIt,1);
c(1)=compute(x0,E,ytruth,funs);
ny=zeros(NIt,size(E,1));
[~,ny(1,:)]=bearingNetworkComputeBearings(x0,E);
epsilon=0.2;
for it=2:NIt
    [c(it),g,ny(it,:)]=compute(x(:,:,it-1),E,ytruth,funs);
    x(:,:,it)=x(:,:,it-1)-epsilon*g;
    m(:,it)=mean(x(:,:,it),2);
end
xFinal=x(:,:,end);
t_node.Ti=xFinal;

figure(1)
semilogy(c)
savefigure(fullfile(figDir,'bearingNetwork_cost'),'epsc',figDim,flagSaveFigure)

figure(2)
plot(x0(1,:),x0(2,:),'rx')
hold on
plot(squeeze(x(1,:,:))',squeeze(x(2,:,:))','color',[1 0.75 0])
bearingNetworkPlot(t_node)
hold off
axis equal
savefigure(fullfile(figDir,'bearingNetwork_trajectories'),'epsc',figDim,flagSaveFigure)

figure(3)
plot(m')
savefigure(fullfile(figDir,'bearingNetwork_centroid'),'epsc',figDim,flagSaveFigure)

figure(4)
plot(ny)
title('Inter-neighbor distances')

function [c,g,ny]=compute(x,E,ytruth,funs)
[y,ny]=bearingNetworkComputeBearings(x,E);
[c,g]=bearingNetworkCost(E,y,ytruth,ny,funs);
