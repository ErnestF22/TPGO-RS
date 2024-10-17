function journal14_singleIntegrators_motion
figDir='../../../papers/control/notes/figures';
figDim=[300,200];
flagSaveFigure=0;
TFinal=7500;

resetRands(1)
allCostNamesBearings={'cosine'};%,'angleSq'};
costNameRanges='squared';

t_node=bearingNetworkBuildTestNetwork();
NNodes=t_node.NNodes;

EBearings=t_node.E;
ygBearings=t_node.Yijtruth;


%idxLeader=4;
idxLeader=[3 4];
dxLeader=@(t) [-5/1000;0.01*sin(t/2/1000*pi)];

xg=t_node.Titruth;
x0=6*randn(2,t_node.NNodes);
%x0=x0+(xg(:,idxLeader)-x0(:,idxLeader))*ones(1,t_node.NNodes);
%x0=x0-mean(x0,2)*ones(1,t_node.NNodes);
x0=xg+0.01*randn(2,t_node.NNodes);
if length(idxLeader)>1
    x0(:,idxLeader(2))=x0(:,idxLeader(1))+(t_node.Titruth(:,idxLeader(2))-t_node.Titruth(:,idxLeader(1)));    
end



NCostBearings=length(allCostNamesBearings);
for iCostBearings=1:NCostBearings
    costNameBearings=allCostNamesBearings{iCostBearings};
    disp(['# Cost ' costNameBearings])
    
    t_cost.funsBearings=bearingCostFunctions(costNameBearings);
    t_cost.alpha=5;
    
    for flagPreconditioner=false%[false true]
        baseFileName=[mfilename '_' costNameBearings '_'];
        if ~flagPreconditioner
            disp('## Bearing formation')
            optsControl={};
            baseFileName=[baseFileName 'b_'];
        else
            disp('## Bearing+preconditioner formation')
            T=bearingNetworkPreconditioner(t_node,funsBearings,'hop1n');
            optsControl={'preconditioner',T,'preconditionerDelay',5};
            baseFileName=[baseFileName 'bp_'];
        end
        
        optsControl=[optsControl {'leader',idxLeader,@(t) dxLeader(t)*ones(size(idxLeader))}];
        figure(1)
        t_node.Ti=x0;
        [t,x,t_node]=bearingNetworkEvolve(t_node,'tFinal',TFinal,'t_cost',t_cost,...
            'showArgs',...
            'optsControl',optsControl,'odeSolver',@ode15s,...%@odeEuler,'optsSolver',{'MaxStep',1.1},...
            'showOdeProgress');
        output=bearingNetworkEvolveStats(t_node,t,x,'cost');

        save([baseFileName 'data'])
        
        %produce figures
        figBaseFileName=fullfile(figDir,baseFileName);
        figure(2)
        bearingNetworkPlot(t_node)
        hold on
        bearingNetworkPlotTrajectories(x)        
        hold off
        savefigure([figBaseFileName 'trajectories'],'epsc',figDim,flagSaveFigure)

        figure(3)
        semilogy(t,output.phi)
        savefigure([figBaseFileName 'cost'],'epsc',figDim,flagSaveFigure)

        figure(4)
        plot(t,output.m)
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
