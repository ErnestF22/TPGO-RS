function bearingNetworkControlDirect_test
figDir='../../../papers/control/notes/figures';
figDim=[300,200];
flagSaveFigure=2;

resetRands(1)
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

TFinal=3000;

for flagUseRanges=[false true]
    baseFileName=fullfile(figDir,'bearingNetwork_');
    if ~flagUseRanges
        disp('# Bearings only')
        dx=@(t,x) control(x,EBearings,ygBearings,funsBearings);
        cEval=@(x) cost(x,EBearings,ygBearings,funsBearings);
        baseFileName=[baseFileName 'B_'];
    else
        disp('# Bearings + ranges')
        dx=@(t,x) control(x,EBearings,ygBearings,funsBearings,ERanges,ygRanges,nygRanges,funsRanges);
        cEval=@(x) cost(x,EBearings,ygBearings,funsBearings,ERanges,ygRanges,nygRanges,funsRanges);
        baseFileName=[baseFileName 'BR_'];
    end

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
    for it=1:Nit
        c(it)=cEval(x(:,:,it));
    end

    figure(2)
    semilogy(c)
    savefigure([baseFileName 'cost'],'epsc',figDim,flagSaveFigure)

    figure(3)
    plot(t,squeeze(sum(x,2))/NNodes)
    ax=axis();
    ax(3:4)=[-1 1];
    axis(ax)
    savefigure([baseFileName 'centroid'],'epsc',figDim,flagSaveFigure)

    save([baseFileName '_data'])
end


function dx=control(x,EBearings,ygBearings,funsBearings,ERanges,ygRanges,nygRanges,funsRanges)
flagUseRanges=false;
if exist('ERanges','var')
    flagUseRanges=true;
end

x=reshape(x,2,[]);
yBearings=bearingNetworkComputeBearings(x,EBearings);

if ~flagUseRanges
    dx=bearingNetworkControlDirect(EBearings,yBearings,ygBearings,funsBearings);
else
    [yRanges,nyRanges]=bearingNetworkComputeBearings(x,ERanges);
    alpha=[1 1];
    dx=bearingNetworkControlDirect(EBearings,yBearings,ygBearings,funsBearings,...
        'ranges',ERanges,yRanges,ygRanges,nyRanges,nygRanges,funsRanges,'alpha',alpha);
end
dx=dx(:);

function c=cost(x,EBearings,ygBearings,funsBearings,ERanges,ygRanges,nygRanges,funsRanges)
flagUseRanges=false;
if exist('ERanges','var')
    flagUseRanges=true;
end

x=reshape(x,2,[]);
[yBearings,nyBearings]=bearingNetworkComputeBearings(x,EBearings);
if ~flagUseRanges
    c=bearingNetworkCost(EBearings,yBearings,ygBearings,nyBearings,funsBearings);
else
    [yRanges,nyRanges]=bearingNetworkComputeBearings(x,ERanges);
    alpha=[1 1];

    c=bearingNetworkCostCombined(...
        EBearings,ERanges,...
        yBearings,yRanges,ygBearings,ygRanges,...
        nyBearings,nyRanges,nygRanges,...
        funsBearings,funsRanges,alpha);
end
