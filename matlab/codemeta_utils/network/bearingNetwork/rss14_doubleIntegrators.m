function rss14_doubleIntegrators
figDir='../../../papers/control/notes/figures';
figDim=[300,200];
flagSaveFigure=0;

resetRands(1)
allCostNamesBearings={'cosine','angleSq'};
costNameRanges='squared';

alpha=[1 1];    %bearing vs ranges
alphaks=[3 6];  %kinetic

lambdas=[0 0.01 0.1 1];

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
%x0=6*randn(2,t_node.NNodes);
x0=xg+3*randn(2,t_node.NNodes);

x0=x0-(mean(x0,2)+offset)*ones(1,t_node.NNodes);

dx0=zeros(size(x0));
z0=[dx0;x0];
flagUseRanges=true;
TFinal=1500;
NCostBearings=length(allCostNamesBearings);
NLambdas=length(lambdas);
NAlphaks=length(alphaks);
for iCostBearings=1:NCostBearings
    costNameBearings=allCostNamesBearings{iCostBearings};
    disp(['# Cost ' costNameBearings])
    funsBearings=bearingCostFunctions(costNameBearings);
    for iAlphak=1:NAlphaks
        alphak=alphaks(iAlphak);
        disp(['## Alphak' num2str(alphak)])
        for iLambda=1:NLambdas
            baseFileName=['bearingNetworkDynamic_' costNameBearings '_'];
            if ~flagUseRanges
                dz=@(t,z)  closedLoop(z,lambdas(iLambda),EBearings,ygBearings,funsBearings,alphak);
                cEval=@(z) cost(z,EBearings,ygBearings,funsBearings,alpha);
                resEval=@(x) residuals(x,EBearings,ygBearings,ERanges,ygRanges,nygRanges);
            else
                dz=@(t,z)  closedLoop(z,lambdas(iLambda),EBearings,ygBearings,funsBearings,alphak,...
                    ERanges,ygRanges,nygRanges,funsRanges,...
                    alpha);
                cEval=@(z) cost(z,EBearings,ygBearings,funsBearings,...
                    ERanges,ygRanges,nygRanges,funsRanges,...
                    alpha);
                resEval=@(x) residuals(x,EBearings,ygBearings,ERanges,ygRanges,nygRanges);
            end            
            baseFileName=[baseFileName num2str(iLambda) '_' num2str(alphak) '_'];

            figure(1)
            optsOde = odeset('OutputFcn',@odeplot);
            [t,z]=ode15s(dz,[0 TFinal],z0(:),optsOde);

            z=reshape(z',4,t_node.NNodes,[]);
            x=z(3:4,:,:);
            xFinal=x(:,:,end);
            t_node.Ti=xFinal;


            %compute statistics
            Nit=length(t);
            NEdgesBearings=size(EBearings,1);
            NEdgesRanges=size(ERanges,1);
            phi=zeros(Nit,1);
            c=zeros(Nit,NEdgesBearings);
            d=zeros(Nit,NEdgesBearings);
            q=zeros(Nit,NEdgesRanges);

            w=getTextWaitBar(Nit);
            w(0)
            for it=1:Nit
                phi(it)=cEval(z(:,:,it));
                [c(it,:),q(it,:)]=resEval(x(:,:,it));
                d(it,:)=bearingNetworkComputeRanges(x(:,:,it),EBearings);
                w(it)
            end
            m=squeeze(sum(x,2))/NNodes;

            save([baseFileName 'data'])

            figBaseFileName=fullfile(figDir,baseFileName);
            figure(1)
            plot(x0(1,:),x0(2,:),'rx')
            hold on
            plot(squeeze(x(1,:,:))',squeeze(x(2,:,:))','-','color',[1 0.75 0])
            plot([xg(1,:); xFinal(1,:)],[xg(2,:); xFinal(2,:)],'k:')
            bearingNetworkPlot(t_node)
            hold off
            axis equal
            savefigure([figBaseFileName 'trajectories'],'epsc',figDim,flagSaveFigure)

            figure(2)
            semilogy(phi)
            savefigure([figBaseFileName 'cost'],'epsc',figDim,flagSaveFigure)

            figure(3)
            plot(t,m)
            ax=axis();
            ax(3:4)=[-1 1];
            axis(ax)
            savefigure([figBaseFileName 'centroid'],'epsc',figDim,flagSaveFigure)
        end
    end
end
%state z:
%   first two rows:     dx
%   second two rows:    x
function dz=closedLoop(z,lambda,EBearings,ygBearings,funsBearings,alphak,...
    ERanges,ygRanges,nygRanges,funsRanges,alpha)
flagUseRanges=false;
if exist('ERanges','var')
    flagUseRanges=true;
end

z=reshape(z,4,[]);
[yBearings,dyBearings]=measurementBearings(z,EBearings);
if ~flagUseRanges
    u=control(EBearings,yBearings,dyBearings,ygBearings,funsBearings,alphak);
else
    [yRanges,nyRanges,dqRanges]=measurementRanges(z,ERanges,ygRanges);
    u=control(EBearings,yBearings,dyBearings,ygBearings,funsBearings,alphak,...
        ERanges,yRanges,dqRanges,ygRanges,nyRanges,nygRanges,funsRanges,alpha);
end

dz=model(z,u,lambda);
%dz=[zeros(2,size(z,2)); u];
dz=dz(:);

function [y,dy]=measurementBearings(z,E)
x=z(3:4,:);
dx=z(1:2,:);
[y,ny]=bearingNetworkComputeBearings(x,E);
dy=bearingNetworkComputeBearingsDerivative(dx,E,y,ny);

function [y,ny,dq]=measurementRanges(z,E,yg)
x=z(3:4,:);
dx=z(1:2,:);
[y,ny]=bearingNetworkComputeBearings(x,E);
dq=bearingNetworkComputeRangeResidualsDerivatives(dx,E,yg);

function u=control(EBearings,yBearings,dyBearings,ygBearings,funsBearings,alphak,...
    ERanges,yRanges,dqRanges,ygRanges,nyRanges,nygRanges,funsRanges,alpha)
flagUseRanges=false;
if exist('ERanges','var')
    flagUseRanges=true;
end

if ~flagUseRanges
    u1=bearingNetworkControlDirect(EBearings,yBearings,ygBearings,funsBearings);
    u2=bearingNetworkControlDynamic(EBearings,dyBearings);
else
    u1=bearingNetworkControlDirect(EBearings,yBearings,ygBearings,funsBearings,...
        'ranges',ERanges,yRanges,ygRanges,nyRanges,nygRanges,funsRanges,'alpha',alpha);
    u2=bearingNetworkControlDynamic(EBearings,dyBearings,...
        'ranges',ERanges,dqRanges,ygRanges,'alpha',alpha);
end    
u=u1+alphak*u2;

function dz=model(z,u,lambda)
m=1;
dz=zeros(size(z));
dz(1:2,:)=1/m*u-lambda*z(1:2,:);
dz(3:4,:)=z(1:2,:);

function c=cost(z,EBearings,ygBearings,funsBearings,ERanges,ygRanges,nygRanges,funsRanges,alpha)
flagUseRanges=false;
if exist('ERanges','var')
    flagUseRanges=true;
end

z=reshape(z,4,[]);
x=z(3:4,:);

[yBearings,nyBearings]=bearingNetworkComputeBearings(x,EBearings);
if ~flagUseRanges
    cPot=bearingNetworkCost(EBearings,yBearings,ygBearings,nyBearings,funsBearings);
else
    [yRanges,nyRanges]=bearingNetworkComputeBearings(x,ERanges);

    cPot=bearingNetworkCostCombined(...
        EBearings,ERanges,...
        yBearings,yRanges,ygBearings,ygRanges,...
        nyBearings,nyRanges,nygRanges,...
        funsBearings,funsRanges,alpha);
end
dx=z(1:2,:);

cKin=0.5*sum(dx(:).^2);
c=cPot+cKin;

function [c,q]=residuals(x,EBearings,ygBearings,ERanges,ygRanges,nygRanges)
flagUseRanges=false;
if exist('ERanges','var')
    flagUseRanges=true;
end
x=reshape(x,2,[]);
yBearings=bearingNetworkComputeBearings(x,EBearings);
c=bearingNetworkComputeBearingsCosines(yBearings,ygBearings);

if flagUseRanges
    [yRanges,nyRanges]=bearingNetworkComputeBearings(x,ERanges);
    q=bearingNetworkComputeRangeResiduals(yRanges,ygRanges,nyRanges,nygRanges);
end