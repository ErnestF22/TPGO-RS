function POCBearingDynamicControlNetwork
resetRands();
figDir='../../../papers/control/notes/figures';
figDim=[300,150];
flagSaveFigure=0;
lambdas=[0];% 0.1 1];
TFinal=30;

costNameBearings='cosine';
costNameRanges='squared';

funsBearings=bearingCostFunctions(costNameBearings);
funsRanges=bearingCostFunctions(costNameRanges);

t_node=bearingNetworkBuildTestNetwork();
NNodes=t_node.NNodes;

EBearings=t_node.E;
ygBearings=t_node.Yijtruth;

ERanges=t_node.Er;
ygRanges=t_node.Yrijtruth;
nygRanges=t_node.nYrijtruth;

xg=t_node.Titruth;
offset=[0;0];
x0=6*randn(2,t_node.NNodes);
x0=x0-(mean(x0,2)+offset)*ones(1,t_node.NNodes);

z0=[zeros(size(x0));x0];

for ilambda=1:length(lambdas)
    disp(['lambda=' num2str(lambdas(ilambda))])
    [t,z]=ode45(@(t,z) closedLoop(z,EBearings,NNodes,ygBearings,funsBearings,lambdas(ilambda)), [0 TFinal], z0);

    Nt=length(t);
    zResh=reshape(z,[],4,NNodes);
    zFinal=squeeze(zResh(end,:,:));%reshape(z(end,:),4,NNodes);
    xFinal=zFinal(3:4,:);
    t_node.Ti=xFinal;

    %compute cost
    c=zeros(Nt,1);
    for it=1:Nt
        c(it)=cost(squeeze(zResh(it,:,:)),EBearings,ygBearings,funsBearings);
    end

    baseFileName=['bearingNetworkDynamic_' num2str(ilambda)];
    
    figure(1)
    subplot(2,1,1)
    plot(t,reshape(zResh(:,1,:),Nt,[]))
    hold on
    plot(t,reshape(zResh(:,2,:),Nt,[]))
    hold off
    ylabel('Velocities')
    xlabel('Time')
    subplot(2,1,2)
    plot(t,reshape(zResh(:,3:4,:),Nt,[]))
    hold on
    plot(t,reshape(zResh(:,3:4,:),Nt,[]))
    hold off
    ylabel('Positions')
    xlabel('Time')
    savefigure(fullfile(figDir,[baseFileName '_plots']),'epsc',figDim,flagSaveFigure)

    figure(2)
    bearingNetworkPlot(t_node)
    hold on
    plot(x0(1,:),x0(2,:),'r*')
    plot([x0(1,:); xFinal(1,:)],[x0(2,:); xFinal(2,:)],'k:')
    plot(squeeze(zResh(:,3,:)), squeeze(zResh(:,4,:)),'Color',[1 0.75 0])
    hold off
    savefigure(fullfile(figDir,[baseFileName '_trajectories']),'epsc',figDim,flagSaveFigure)

    figure(3)
    plot(t,c)
    ylabel('Energy')
    xlabel('Time')
    savefigure(fullfile(figDir,[baseFileName '_energy']),'epsc',figDim,flagSaveFigure)
end

function c=cost(z,E,YTruth,funs)
x=z(3:4,:);
nYij= bearingNetworkComputeRanges(x,E);
Yij= bearingNetworkComputeBearings(x,E);
cPot=bearingNetworkCost(E,Yij,YTruth,nYij,funs);
dx=z(1:2,:);
cKin=sum(dx(:).^2);
c=cPot+cKin;

%state z:
%   first two rows:     dx
%   second two rows:    x
function dz=closedLoop(z,E,NNodes,YTruth,funs,lambda)
z=reshape(z,4,[]);
[Y,dY]=measurement(z,E);

u=control(Y,YTruth,dY,E,NNodes,funs);
dz=model(z,u,lambda);
dz=dz(:);

function [Y,dY]=measurement(z,E)
nY=bearingNetworkComputeRanges(z(3:4,:),E);
Y=bearingNetworkComputeBearings(z(3:4,:),E);
dY=bearingNetworkDerivative(z(1:2,:),Y,nY,E);

function u1=control1(y,yg,E,funs)
u1=-bearingNetworkCostGradient(E,y,yg,funs);

function u2=control2(Y,dY,NNodes,E)
B=bearingBuildB(Y,NNodes,E);
u2=-reshape(B'*dY(:),2,[]);

function u=control(y,yg,dy,E,NNodes,funs)
alpha=3;
u1=control1(y,yg,E,funs);
u2=control2(y,dy,NNodes,E);
u=alpha*u1+u2;

function dz=model(z,u,lambda)
m=1;
dz=zeros(size(z));
dz(1:2,:)=1/m*u-lambda*z(1:2,:);
dz(3:4,:)=z(1:2,:);

function dYij=bearingNetworkDerivative(v,Yij,nYij,E)
[d,NNodes]=size(v);
R=bearingBuildR(nYij,d);
B=bearingBuildB(Yij,NNodes,E);
dYij=reshape(R\B*v(:),2,[]);

function R=bearingBuildR(nYij,d)
R=kron(diag(nYij),eye(d));

function B=bearingBuildB(Yij,NNodes,E)
d=size(Yij,1);
NEdges=size(E,1);
idxNodes=reshape(1:d*NNodes,d,NNodes);
idxEdges=reshape(1:d*NEdges,d,NEdges);
B=zeros(d*NEdges,d*NNodes);
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    
    PYij=eye(2)-Yij(:,iEdge)*Yij(:,iEdge)';
    B(idxEdges(:,iEdge),idxNodes(:,iNode))=-PYij;
    B(idxEdges(:,iEdge),idxNodes(:,jNode))=PYij;
end
