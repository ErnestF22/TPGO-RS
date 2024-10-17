function presentation_singleIntegrators_motion_proj
figDir='../../../papers/control/notes/figures';
figDim=[300,200];
flagSaveFigure=0;
TFinal=7500;

resetRands(1)
allCostNamesBearings={'cosine'};%,'angleSq'};

t_node=bearingNetworkBuildTestNetwork();
NNodes=t_node.NNodes;

EBearings=t_node.E;
ygBearings=t_node.Yijtruth;

ERanges=t_node.Er;

idxLeader=5;
nEdgeBearingProject=8;
dxLeader=@(t) [-5/1000;0.01*sin(t/2/1000*pi)];

xg=t_node.Titruth;
x0=6*randn(2,t_node.NNodes);
%x0=x0+(xg(:,idxLeader)-x0(:,idxLeader))*ones(1,t_node.NNodes);
%x0=x0-mean(x0,2)*ones(1,t_node.NNodes);
x0=xg+0.01*randn(2,t_node.NNodes);


NCostBearings=length(allCostNamesBearings);
for iCostBearings=1:NCostBearings
    costNameBearings=allCostNamesBearings{iCostBearings};
    disp(['# Cost ' costNameBearings])
    funsBearings=bearingCostFunctions(costNameBearings);
    flagUseRanges=false;
    baseFileName=['bearingNetworkMotion_' costNameBearings '_'];
    disp('## Pure bearing formation')
    %dx=@(t,x) control(t,x,EBearings,ygBearings,funsBearings);
    dx=@(t,x) closedLoop(x,xg,EBearings,funsBearings,'leader',idxLeader,dxLeader(t));
    resEval=@(t,x) residuals(x,xg,EBearings);
    baseFileName=[baseFileName 'b_'];

    figure(1)
    optsOde = odeset('OutputFcn',@odeplot);
    [t,x]=ode45(dx,[0 TFinal],x0(:),optsOde);

    x=reshape(x',2,t_node.NNodes,[]);
    xFinal=x(:,:,end);
    t_node.Ti=xFinal;


    %compute statistics
    Nit=length(t);
    NEdgesBearings=size(EBearings,1);
    phi=zeros(Nit,1);
    c=zeros(Nit,NEdgesBearings);
    d=zeros(Nit,NEdgesBearings);

    w=getTextWaitBar(Nit);
    w(0)
    for it=1:Nit
        c(it,:)=resEval(t(it),x(:,:,it));
        d(it,:)=bearingNetworkComputeRanges(x(:,:,it),EBearings);
        w(it)
    end
    m=squeeze(sum(x,2))/NNodes;

    save([baseFileName 'data'])

    figBaseFileName=fullfile(figDir,baseFileName);
    figure(1)
    plot(x0(1,:),x0(2,:),'rx')
    hold on
    plot(squeeze(x(1,:,:))',squeeze(x(2,:,:))','-.','color',[1 0.75 0])
    plot([xg(1,:); xFinal(1,:)],[xg(2,:); xFinal(2,:)],'k:')
    bearingNetworkPlot(t_node)
    hold off
    axis equal
    savefigure([figBaseFileName 'trajectories'],'epsc',figDim,flagSaveFigure)

    figure(2)
    semilogy(t,funsBearings.f(c))

    figure(3)
    plot(t,m)
    savefigure([figBaseFileName 'centroid'],'epsc',figDim,flagSaveFigure)
end

function dx=closedLoop(x,xg,EBearings,funsBearings,varargin)
flagUseLeader=false;
nEdgeBearingProject=1;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'leader'
            flagUseLeader=true;
            ivarargin=ivarargin+1;
            idxLeader=varargin{ivarargin};
            ivarargin=ivarargin+1;
            dxLeader=varargin{ivarargin};
            NNodes=size(xg,2);
            idxNodes=reshape(1:2*NNodes,size(xg));
        case 'project'
            ivarargin=ivarargin+1;
            nEdgeBearingProject=varargin{ivarargin};
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

ygBearings=bearingNetworkComputeBearings(xg,EBearings);
dx=control(x,EBearings,ygBearings,funsBearings,nEdgeBearingProject);

if flagUseLeader
    dx(idxNodes(:,idxLeader))=dxLeader;
end

function dx=control(x,EBearings,ygBearings,funsBearings,nEdgeBearingProject)

x=reshape(x,2,[]);
yBearings=bearingNetworkComputeBearings(x,EBearings);

dx=bearingNetworkControlDirect(EBearings,yBearings,ygBearings,funsBearings,...
    'alpha',5,'project',nEdgeBearingProject);
dx=dx(:);

function [c,q]=residuals(x,xg,EBearings,ERanges)
flagUseRanges=false;
if exist('ERanges','var')
    flagUseRanges=true;
end
x=reshape(x,2,[]);
yBearings=bearingNetworkComputeBearings(x,EBearings);
ygBearings=bearingNetworkComputeBearings(xg,EBearings);

c=bearingNetworkComputeBearingsCosines(yBearings,ygBearings);

if flagUseRanges
    [yRanges,nyRanges]=bearingNetworkComputeBearings(x,ERanges);
    [ygRanges,nygRanges]=bearingNetworkComputeBearings(xg,ERanges);
    q=bearingNetworkComputeRangeResiduals(yRanges,ygRanges,nyRanges,nygRanges);
end

% function c=cost(x,EBearings,ygBearings,funsBearings,ERanges,ygRanges,nygRanges,funsRanges)
% flagUseRanges=false;
% if exist('ERanges','var')
%     flagUseRanges=true;
% end
% 
% x=reshape(x,2,[]);
% [yBearings,nyBearings]=bearingNetworkComputeBearings(x,EBearings);
% if ~flagUseRanges
%     c=bearingNetworkCost(EBearings,yBearings,ygBearings,nyBearings,funsBearings);
% else
%     [yRanges,nyRanges]=bearingNetworkComputeBearings(x,ERanges);
%     alpha=[1 1];
% 
%     c=bearingNetworkCostCombined(...
%         EBearings,ERanges,...
%         yBearings,yRanges,ygBearings,ygRanges,...
%         nyBearings,nyRanges,nygRanges,...
%         funsBearings,funsRanges,alpha);
% end
