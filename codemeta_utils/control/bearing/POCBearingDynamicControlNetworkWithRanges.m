function POCBearingDynamicControlNetworkWithRanges
resetRands()
figDir='../../../papers/control/notes/figures';
figDim=[300,150];
flagSaveFigure=0;
flagShowCost=false;
%methodIntegration='euler'; dt=0.002;
methodIntegration='ode45';


lambdas=1;%[0.01 0.1 1];
TFinal=10;

alpha=[1 5];    %bearing vs ranges
alphak=1;       %kinetic

%costNameBearings='cosine';
costNameBearings='angleSq';
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
dx0=zeros(size(x0));
% x0=xg;
% dx0=5*rand(size(x0));
% dx0=dx0-mean(dx0,2)*ones(1,t_node.NNodes);

z0=[dx0;x0];

for ilambda=1:length(lambdas)
    disp(['lambda=' num2str(lambdas(ilambda))])
    dz=@(t,z)  closedLoop(z,EBearings,ygBearings,funsBearings,...
        ERanges,ygRanges,nygRanges,funsRanges,...
        lambdas(ilambda),alpha,alphak);
    cEval=@(z) cost(z,EBearings,ygBearings,funsBearings,...
        ERanges,ygRanges,nygRanges,funsRanges,...
        alpha);

    figure(1)
    subplot(1,1,1)
    switch methodIntegration
        case 'euler'
            t=0:dt:TFinal;
            Nt=length(t);
            z=zeros([Nt numel(z0)]);
            tCurrent=0;
            zCurrent=z0(:);
            z(1,:)=zCurrent;
            set(gcf,'NextPlot','ReplaceChildren')
            for it=2:Nt
                dCurrent=dz(tCurrent,zCurrent);
                tCurrent=tCurrent+dt;
                zCurrent=zCurrent+dt*dCurrent(:);
                t(it)=tCurrent;
                z(it,:)=zCurrent;
                plot(t(1:it),z(1:it,:),'EraseMode','none')
                ax=axis;
                ax(2)=TFinal;
                axis(ax);
                drawnow expose;
            end
        case 'ode45'
            optsOde = odeset('OutputFcn',@odeplot);
            [t,z]=ode45(dz,[0 TFinal],z0,optsOde);
    end
%     t=linspace(0,0.1,100);
%     Nt=length(t);
%     z=zeros(Nt,4*t_node.NNodes);
%     z(1,:)=z0(:);
%     for it=2:Nt
%         d=dz(t(it-1),z(it-1,:));
%         z(it,:)=z(it-1,:)+0.01*(t(it)-t(it-1)*d');
%     end
    
    Nt=length(t);
    zResh=reshape(z,[],4,NNodes);
    zFinal=squeeze(zResh(end,:,:));%reshape(z(end,:),4,NNodes);
    xFinal=zFinal(3:4,:);
    t_node.Ti=xFinal;

    if flagShowCost
        %compute cost
        c=zeros(Nt,1);
        dc=zeros(Nt,1);
        w=getTextWaitBar(Nt);
        w(0)
        for it=1:Nt
            [c(it),gradc]=cEval(squeeze(zResh(it,:,:)));
            dc(it)=gradc(:)'*dz(t(it),reshape(zResh(it,:,:),[],1));
            w(it)
        end
    end

    baseFileName=['bearingNetworkDynamic_T' num2str(TFinal) '_l' num2str(ilambda)];
    save([baseFileName '_data'])
    
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
    plot([xg(1,:); xFinal(1,:)],[xg(2,:); xFinal(2,:)],'k:')
    plot([xg(1,:); x0(1,:)],[xg(2,:); x0(2,:)],'k:')
    plot(squeeze(zResh(:,3,:)), squeeze(zResh(:,4,:)))%,'-','Color',[1 0.75 0])
    hold off
    savefigure(fullfile(figDir,[baseFileName '_trajectories']),'epsc',figDim,flagSaveFigure)

    if flagShowCost
        figure(3)
    %     semilogy(t,c)
        subplot(2,1,1)
        plot(t,c)
        subplot(2,1,2)
        plot(t,dc,'r')
        ylabel('Energy')
        xlabel('Time')
    end
    savefigure(fullfile(figDir,[baseFileName '_energy']),'epsc',figDim,flagSaveFigure)
end

function [c,dc]=costAndDerPot(z,dz,EBearings,ygBearings,funsBearings,ERanges,ygRanges,nygRanges,funsRanges,alpha)
x=z(3:4,:);
x2=z(1:2,:);
dx=dz(3:4,:);
dx2=dz(1:2,:);
[yBearings,nyBearings]=bearingNetworkComputeBearings(x,EBearings);
[yRanges,nyRanges]=bearingNetworkComputeBearings(x,ERanges);

[c,gradPot]=bearingNetworkCostCombined(...
    EBearings,ERanges,...
    yBearings,yRanges,ygBearings,ygRanges,...
    nyBearings,nyRanges,nygRanges,...
    funsBearings,funsRanges,alpha);

c=c+0.5*x2(:)'*x2(:);

dc=[x2(:);gradPot(:)]'*[dx2(:);dx(:)];

function [c,dc]=costAndDer(cEval,z,dz)
dz=reshape(dz,4,[]);
[c,gradc]=cEval(z);
dx=dz(3:4,:);
dx2=dz(1:2,:);
dc=gradc(:)'*[dx2(:);dx(:)];

function [c,gradc]=cost(z,EBearings,ygBearings,funsBearings,ERanges,ygRanges,nygRanges,funsRanges,alpha)
x=z(3:4,:);

[yBearings,nyBearings]=bearingNetworkComputeBearings(x,EBearings);
[yRanges,nyRanges]=bearingNetworkComputeBearings(x,ERanges);

[cPot,gradPot]=bearingNetworkCostCombined(...
    EBearings,ERanges,...
    yBearings,yRanges,ygBearings,ygRanges,...
    nyBearings,nyRanges,nygRanges,...
    funsBearings,funsRanges,alpha);

% u1=bearingNetworkControlDirect(EBearings,yBearings,ygBearings,funsBearings,...
%     ERanges,yRanges,ygRanges,nyRanges,nygRanges,funsRanges,alpha);

dx=z(1:2,:);
% dyBearings=bearingNetworkDerivative(dx,yBearings,nyBearings,EBearings);
% dqRanges=bearingNetworkComputeRangeResidualsDerivatives(dx,yRanges,ygRanges,ERanges);

cKin=0.5*sum(dx(:).^2);
c=cPot+cKin;

% u2=bearingNetworkControlDynamic(EBearings,dyBearings,...
%     ERanges,dqRanges,ygRanges,alpha);
% 
% u1b=bearingNetworkControlDirect(EBearings,yBearings,ygBearings,funsBearings,...
%     ERanges,yRanges,ygRanges,nyRanges,nygRanges,funsRanges,alpha);
% u2b=bearingNetworkControlDynamic(EBearings,dyBearings,...
%     ERanges,dqRanges,ygRanges,alpha);
% 
% u=alphak(1)*u1b+alphak(2)*u2b;
% dc=dx(:)'*gradPot(:)+dx(:)'*u(:);
gradc=[dx;gradPot];
%gradc=[zeros(size(dx(:)));alphak(1)*gradPot(:)];

%state z:
%   first two rows:     dx
%   second two rows:    x
function dz=closedLoop(z,EBearings,ygBearings,funsBearings,...
    ERanges,ygRanges,nygRanges,funsRanges,lambda,alpha,alphak)
z=reshape(z,4,[]);
[yBearings,dyBearings]=measurementBearings(z,EBearings);
[yRanges,nyRanges,dqRanges]=measurementRanges(z,ERanges,ygRanges);
u=control(EBearings,yBearings,dyBearings,ygBearings,funsBearings,...
    ERanges,yRanges,dqRanges,ygRanges,nyRanges,nygRanges,funsRanges,alpha,alphak);

dz=model(z,u,lambda);
%dz=[zeros(2,size(z,2)); u];
dz=dz(:);

function [y,dy]=measurementBearings(z,E)
x=z(3:4,:);
dx=z(1:2,:);
[y,ny]=bearingNetworkComputeBearings(x,E);
dy=bearingNetworkDerivative(dx,y,ny,E);

function [y,ny,dq]=measurementRanges(z,E,yg)
x=z(3:4,:);
dx=z(1:2,:);
[y,ny]=bearingNetworkComputeBearings(x,E);
dq=bearingNetworkComputeRangeResidualsDerivatives(dx,y,yg,E);

function u=control(EBearings,yBearings,dyBearings,ygBearings,funsBearings,...
    ERanges,yRanges,dqRanges,ygRanges,nyRanges,nygRanges,funsRanges,alpha,alphak)
u1=bearingNetworkControlDirect(EBearings,yBearings,ygBearings,funsBearings,...
    ERanges,yRanges,ygRanges,nyRanges,nygRanges,funsRanges,alpha);
u2=bearingNetworkControlDynamic(EBearings,dyBearings,...
    ERanges,dqRanges,ygRanges,alpha);

u=u1+alphak*u2;

function dz=model(z,u,lambda)
m=1;
dz=zeros(size(z));
dz(1:2,:)=1/m*u-lambda*z(1:2,:);
dz(3:4,:)=z(1:2,:);

function dy=bearingNetworkDerivative(v,y,ny,E)
[d,NNodes]=size(v);
R=bearingBuildR(ny,d);
B=bearingBuildB(y,NNodes,E);
dy=reshape(R\B*v(:),2,[]);

function R=bearingBuildR(ny,d)
R=kron(diag(ny),eye(d));

function B=bearingBuildB(y,NNodes,E)
d=size(y,1);
NEdges=size(E,1);
idxNodes=reshape(1:d*NNodes,d,NNodes);
idxEdges=reshape(1:d*NEdges,d,NEdges);
B=zeros(d*NEdges,d*NNodes);
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    
    PYij=eye(2)-y(:,iEdge)*y(:,iEdge)';
    B(idxEdges(:,iEdge),idxNodes(:,iNode))=-PYij;
    B(idxEdges(:,iEdge),idxNodes(:,jNode))=PYij;
end
