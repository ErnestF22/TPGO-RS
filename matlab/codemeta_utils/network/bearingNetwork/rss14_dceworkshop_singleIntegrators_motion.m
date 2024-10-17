function rss14_dceworkshop_singleIntegrators_motion
figDir='../../../papers/control/notes/figures';
figDim=[300,200];
flagSaveFigure=0;
TFinal=7500;

resetRands(1)
allCostNamesBearings={'cosine','angleSq'};
costNameRanges='squared';

funsRanges=bearingCostFunctions(costNameRanges);

t_node=bearingNetworkBuildTestNetwork();
NNodes=t_node.NNodes;

EBearings=t_node.E;
ygBearings=t_node.Yijtruth;

ERanges=t_node.Er;
ygRanges=t_node.Yrijtruth;
nygRanges=t_node.nYrijtruth;

idxLeader=4;
dxLeader=@(t) [-5/1000;0.01*sin(t/2/1000*pi)];

xg=t_node.Titruth;
x0=6*randn(2,t_node.NNodes);
x0=x0+(xg(:,idxLeader)-x0(:,idxLeader))*ones(1,t_node.NNodes);
%x0=x0-mean(x0,2)*ones(1,t_node.NNodes);
%x0=xg+0.01*randn(2,t_node.NNodes);


NCostBearings=length(allCostNamesBearings);
for iCostBearings=1:NCostBearings
    costNameBearings=allCostNamesBearings{iCostBearings};
    disp(['# Cost ' costNameBearings])
    funsBearings=bearingCostFunctions(costNameBearings);
    for flagUseRanges=[false true]
        baseFileName=[mfilename '_' costNameBearings '_'];
        if ~flagUseRanges
            disp('## Pure bearing formation')
            %dx=@(t,x) control(t,x,EBearings,ygBearings,funsBearings);
            dx=@(t,x) closedLoop(x,xg,EBearings,funsBearings,'leader',idxLeader,dxLeader(t));
            resEval=@(t,x) residuals(x,xg,EBearings);
            baseFileName=[baseFileName 'b_'];
        else
            disp('## Bearing+distance formation')
            %dx=@(t,x) control(t,x,EBearings,ygBearings,funsBearings,ERanges,ygRanges,nygRanges,funsRanges);
            dx=@(t,x) closedLoop(x,xg,EBearings,funsBearings,'ranges',ERanges,funsRanges,...
                'leader',idxLeader,dxLeader(t));
            resEval=@(t,x) residuals(x,xg,EBearings,ERanges);
            baseFileName=[baseFileName 'bd_'];
        end

        figure(1)
        optsOde = odeset('OutputFcn',@odeplot);
        [t,x]=ode15s(dx,[0 TFinal],x0(:),optsOde);

        x=reshape(x',2,t_node.NNodes,[]);
        xFinal=x(:,:,end);
        t_node.Ti=xFinal;


        %compute statistics
        Nit=length(t);
        NEdgesBearings=size(EBearings,1);
        phi=zeros(Nit,1);
        c=zeros(Nit,NEdgesBearings);
        d=zeros(Nit,NEdgesBearings);
        if flagUseRanges
            NEdgesRanges=size(ERanges,1);
            q=zeros(Nit,NEdgesRanges);
        end
            
        w=getTextWaitBar(Nit);
        w(0)
        for it=1:Nit
            if ~flagUseRanges
                c(it,:)=resEval(t(it),x(:,:,it));
            else
                [c(it,:),q(it,:)]=resEval(t(it),x(:,:,it));
            end
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
end

function dx=closedLoop(x,xg,EBearings,funsBearings,varargin)
flagUseRanges=false;
flagUseLeader=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'ranges'
            flagUseRanges=true;
            ivarargin=ivarargin+1;
            ERanges=varargin{ivarargin};
            ivarargin=ivarargin+1;
            funsRanges=varargin{ivarargin};
        case 'leader'
            flagUseLeader=true;
            ivarargin=ivarargin+1;
            idxLeader=varargin{ivarargin};
            ivarargin=ivarargin+1;
            dxLeader=varargin{ivarargin};
            NNodes=size(xg,2);
            idxNodes=reshape(1:2*NNodes,size(xg));
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

ygBearings=bearingNetworkComputeBearings(xg,EBearings);
if ~flagUseRanges
    dx=control(x,EBearings,ygBearings,funsBearings);
else
    [ygRanges,nygRanges]=bearingNetworkComputeBearings(xg,ERanges);
    dx=control(x,EBearings,ygBearings,funsBearings,ERanges,ygRanges,nygRanges,funsRanges);
end

if flagUseLeader
    dx(idxNodes(:,idxLeader))=dxLeader;
end

function dx=control(x,EBearings,ygBearings,funsBearings,ERanges,ygRanges,nygRanges,funsRanges)
flagUseRanges=false;
if exist('ERanges','var')
    flagUseRanges=true;
end

x=reshape(x,2,[]);
yBearings=bearingNetworkComputeBearings(x,EBearings);

if ~flagUseRanges
    dx=bearingNetworkControlDirect(EBearings,yBearings,ygBearings,funsBearings,...
        'alpha',5);
else
    [yRanges,nyRanges]=bearingNetworkComputeBearings(x,ERanges);
    alpha=[5 5];
    dx=bearingNetworkControlDirect(EBearings,yBearings,ygBearings,funsBearings,...
        'ranges',ERanges,yRanges,ygRanges,nyRanges,nygRanges,funsRanges,'alpha',alpha);
end
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
