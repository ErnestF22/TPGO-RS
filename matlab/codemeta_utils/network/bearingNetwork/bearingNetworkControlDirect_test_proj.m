function bearingNetworkControlDirect_test_proj
figDir='../../../papers/control/notes/figures';
figDim=[300,200];
flagSaveFigure=0;

resetRands(1)
costNameBearings='cosine';

funsBearings=bearingCostFunctions(costNameBearings);

t_node=bearingNetworkBuildTestNetwork();
NNodes=t_node.NNodes;

EBearings=t_node.E;
ygBearings=t_node.Yijtruth;
nygBearings=t_node.nYijtruth;

xg=t_node.Titruth;
offset=[0;0];
x0=6*randn(2,t_node.NNodes);
x0=x0-(mean(x0,2)+offset)*ones(1,t_node.NNodes);

%change initial conditions to match final scale
nEdgeBearingProject=8;
[~,ny0Bearings]=bearingNetworkComputeBearings(x0,EBearings);
x0=x0/ny0Bearings(nEdgeBearingProject)*nygBearings(nEdgeBearingProject);
TFinal=5000;

for flagUseRanges=false %[false true]
    baseFileName=fullfile(figDir,'bearingNetwork_proj_');
    disp('# Bearings only')
    dx=@(t,x) control(x,EBearings,ygBearings,funsBearings,nEdgeBearingProject);
    cEval=@(x) cost(x,EBearings,ygBearings,funsBearings,nEdgeBearingProject);
    baseFileName=[baseFileName 'B_'];

    figure(1)
    optsOde = odeset('OutputFcn',@odeplot);
    [t,x]=ode45(dx,[0 TFinal],x0(:),optsOde);

    x=reshape(x',2,t_node.NNodes,[]);
    xFinal=x(:,:,end);
    t_node.Ti=xFinal;

    figure(1)
    plot(x0(1,:),x0(2,:),'rx')
    hold on
    plot(squeeze(x(1,:,:))',squeeze(x(2,:,:))','-.','color',[1 0.75 0])
    plot([xg(1,:); xFinal(1,:)],[xg(2,:); xFinal(2,:)],'k:')
    bearingNetworkPlot(t_node)
    hold off
    axis equal
    savefigure([baseFileName 'trajectories'],'epsc',figDim,flagSaveFigure)


    Nit=length(t);
    c=zeros(Nit,1);
    nyProject=zeros(Nit,1);
    for it=1:Nit
        [c(it),nyProject(it)]=cEval(x(:,:,it));
    end

    figure(2)
    semilogy(c)
    savefigure([baseFileName 'cost'],'epsc',figDim,flagSaveFigure)

    figure(3)
    plot(t,nyProject)
    
    figure(4)
    plot(t,squeeze(sum(x,2))/NNodes)
    ax=axis();
    ax(3:4)=[-1 1];
    axis(ax)
    savefigure([baseFileName 'centroid'],'epsc',figDim,flagSaveFigure)

    save([baseFileName '_data'])
end


function dx=control(x,EBearings,ygBearings,funsBearings,nEdgeBearingProject)

x=reshape(x,2,[]);
yBearings=bearingNetworkComputeBearings(x,EBearings);

dx=bearingNetworkControlDirect(EBearings,yBearings,ygBearings,funsBearings,'project',nEdgeBearingProject);
dx=dx(:);

function [c,nyProject]=cost(x,EBearings,ygBearings,funsBearings,nEdgeBearingProject)

x=reshape(x,2,[]);
[yBearings,nyBearings]=bearingNetworkComputeBearings(x,EBearings);
nyProject=nyBearings(nEdgeBearingProject);
c=bearingNetworkCost(EBearings,yBearings,ygBearings,nyBearings,funsBearings);
