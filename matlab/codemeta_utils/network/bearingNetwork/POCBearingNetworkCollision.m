%Test simple case of collisions between agents
function POCBearingNetworkCollision
resetRands(1)
allCostNamesBearings={{'squared','field'},{'angleSq','gradient'},{'cosine','gradient'}};
costNameRanges='squared';

funsRanges=bearingCostFunctions(costNameRanges);

t_node=bearingNetworkBuildTestNetwork();
NNodes=t_node.NNodes;

EBearings=t_node.E;
ygBearings=t_node.Yijtruth;

ERanges=t_node.Er;
ygRanges=t_node.Yrijtruth;
nygRanges=t_node.nYrijtruth;


NSubset=[4  6 7];

NNodes=length(NSubset);
idxValidE=ismember(EBearings(:,1),NSubset) & ismember(EBearings(:,2),NSubset);
EBearings=mapValues(EBearings(idxValidE,:),[NSubset;1:NNodes]');
xg=t_node.Titruth(:,NSubset);

x0=[0 0.1 0.2; -2.4 -0.9 0];
xg=[-2 -2 0; 2 0 0];
ygBearings=bearingNetworkComputeBearings(xg,EBearings);

t_node=[];
TFinal=1.5;
NCostBearings=length(allCostNamesBearings);
for iCostBearings=1%1:NCostBearings
    costNameBearings=allCostNamesBearings{iCostBearings}{1};
    methodField=allCostNamesBearings{iCostBearings}{2};
    disp(['# Cost ' costNameBearings])
    funsBearings=bearingCostFunctions(costNameBearings);
    for flagUseRanges=[false true]
        baseFileName=[mfilename '_' costNameBearings '_'];
        if ~flagUseRanges
            disp('## Pure bearing formation')
            dx=@(t,x) control(x,EBearings,ygBearings,funsBearings,methodField);
            cEval=@(x) cost(x,EBearings,ygBearings,funsBearings);
            resEval=@(x) residuals(x,EBearings,ygBearings);
            baseFileName=[baseFileName 'b_'];
        else
            disp('## Bearing+distance formation')
            dx=@(t,x) control(x,EBearings,ygBearings,funsBearings,methodFieldERanges,ygRanges,nygRanges,funsRanges);
            cEval=@(x) cost(x,EBearings,ygBearings,funsBearings,ERanges,ygRanges,nygRanges,funsRanges);
            resEval=@(x) residuals(x,EBearings,ygBearings,ERanges,ygRanges,nygRanges);
            baseFileName=[baseFileName 'bd_'];
        end

        figure(1)
        optsOde = odeset('OutputFcn',@odeplot,'MaxStep',1e-2,'RelTol',1e-6);
        [t,x]=ode45(dx,[0 TFinal],x0(:),optsOde);

        x=reshape(x',2,NNodes,[]);
        xFinal=x(:,:,end);
        %t_node.Ti=xFinal;


        %compute statistics
        Nit=length(t);
        NEdgesBearings=size(EBearings,1);
        phi=zeros(Nit,1);
        c=zeros(Nit,NEdgesBearings);
        d=zeros(Nit,NEdgesBearings);
        u=zeros([Nit size(x0)]);
        if flagUseRanges
            NEdgesRanges=size(ERanges,1);
            q=zeros(Nit,NEdgesRanges);
        end
            
        w=getTextWaitBar(Nit);
        w(0)
        for it=1:Nit
            phi(it)=cEval(x(:,:,it));
            if ~flagUseRanges
                c(it,:)=resEval(x(:,:,it));
            else
                [c(it,:),q(it,:)]=resEval(x(:,:,it));
            end
            d(it,:)=bearingNetworkComputeRanges(x(:,:,it),EBearings);
            u(it,:)=dx(0,x(:,:,it));
            w(it)
        end
        m=squeeze(sum(x,2))/NNodes;
        
        save([baseFileName 'data'])

        figure(1)
        plot(x0(1,:),x0(2,:),'rx')
        hold on
        plot(squeeze(x(1,:,:))',squeeze(x(2,:,:))','-.','color',[1 0.75 0])
        plot([xg(1,:); x0(1,:)],[xg(2,:); x0(2,:)],'k:')
        plot(xg(1,:),xg(2,:),'o');
        for iNode=1:NNodes
            text(xg(1,iNode),xg(2,iNode),num2str(iNode));
        end
        %bearingNetworkPlot(t_node)
        hold off
        axis equal

        figure(2)
        semilogy(phi)

        figure(3)
        plot(t,d)
        disp(min(d(:)))
        return
    end
end

function dx=control(x,EBearings,ygBearings,funsBearings,methodField,ERanges,ygRanges,nygRanges,funsRanges)
flagUseRanges=false;
if exist('ERanges','var')
    flagUseRanges=true;
end

x=reshape(x,2,[]);
yBearings=bearingNetworkComputeBearings(x,EBearings);

if ~flagUseRanges
    dx=bearingNetworkControlDirect(EBearings,yBearings,ygBearings,funsBearings,...
        'methodBearing',methodField);
else
    [yRanges,nyRanges]=bearingNetworkComputeBearings(x,ERanges);
    alpha=[1 1];
    dx=bearingNetworkControlDirect(EBearings,yBearings,ygBearings,funsBearings,...
        'ranges',ERanges,yRanges,ygRanges,nyRanges,nygRanges,funsRanges,'alpha',alpha,...
        'methodBearing',methodField);
end
dx(:,3)=0;
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