function POClocalization_coordinateDescent
flagCollectErrors=false;
flagDisplayIt=false;
flagProgressBar=true;
if flagDisplayIt
    flagProgressBar=false;
end

t_node=testNetworkBuildTestNetwork('methodInit','noisytruth',...
    'varNoisyTruth',[0.5 0.5]);

NIt=1000;
E=testNetworkGetEdges(t_node);
A=t_node.A;
N=testNetworkGetNumberOfNodes(t_node);
NEdges=size(E,1);
[~,Qij]=testNetworkGetRelativeEssential(t_node);
% [Ri0,Ti0]=testNetworkGetRotTransl(t_node);
% GiInit=RT2G(Ri0,Ti0);
GiInit=RT2G(rot_randn([],[],N), 8*randn(3,N));

GiTruth=t_node.gitruth;

neigh=cell(1,N);
ENeigh=cell(1,N);
QNeigh=cell(1,N);
for iN=1:N
    jNeigh=find(A(iN,:));
    neigh{iN}=jNeigh;
    flagENeigh=E(:,1)==iN ...
        & any((E(:,2)*ones(1,length(jNeigh)))==(ones(NEdges,1)*jNeigh),2);
    ENeigh{iN}=E(flagENeigh,:);
    QNeigh{iN}=Qij(:,:,flagENeigh);
end

Gi=GiInit;
if flagCollectErrors
    c=zeros(NIt*N+1,1);
    t=zeros(NIt*N+1,1);
    c(1)=essentialCostNetwork(G2R(Gi),G2T(Gi),Qij,E);
    cnt=2;
end
if flagProgressBar
    w=getTextWaitBar(NIt);
end
tStart=cputime;
for it=1:NIt
    for iN=1:N
        GiNeigh=Gi(:,:,neigh{iN});
        QiNeigh=QNeigh{iN};
        Gi(:,:,iN)=singleUpdate(Gi(:,:,iN),GiNeigh,QiNeigh);
        if flagDisplayIt
            fprintf('it:%4d n:%3d\n',it,iN);
        end
        if flagCollectErrors
            cCurrent=essentialCostNetwork(G2R(Gi),G2T(Gi),Qij,E);
            tCurrent=cputime-tStart;
            c(cnt)=cCurrent;
            t(cnt)=tCurrent;
            cnt=cnt+1;
            if flagDisplayIt
                fprintf('t:%4.1f c:%.4e\n',tCurrent,cCurrent);
            end
        end
        if flagDisplayIt
            fprintf('\n');
        end
    end
    if flagProgressBar
        w(it)
    end
end
fprintf('cputime: %.2f\n',cputime-tStart);

GiEst=Gi;

figure(1)
localization_display(GiTruth,GiInit,GiEst,t_node.A)

figure(2)
semilogy(t,c)

function c=singleCost(Gi,GiNeigh,QNeigh)
NNeigh=size(GiNeigh,3);
[Ri,Ti]=G2RT(Gi);
c=0;
for iNeigh=1:NNeigh
    [Rj,Tj]=G2RT(GiNeigh(:,:,iNeigh));
    Qij=QNeigh(:,:,iNeigh);
    c=c+essentialCost(Ri,Ti,Rj,Tj,Qij);
end

function Gi=singleUpdate(Gi,GiNeigh,QNeigh)
NNeigh=size(GiNeigh,3);
[Ri,Ti]=G2RT(Gi);
gradc=0;
for iNeigh=1:NNeigh
    [Rj,Tj]=G2RT(GiNeigh(:,:,iNeigh));
    Qij=QNeigh(:,:,iNeigh);
    [~,gradcij]=essentialCost(Ri,Ti,Rj,Tj,Qij);
    gradc=gradc+gradcij(1:6);
end
d=-gradc/NNeigh;
Gi=rot3r3_exp(Gi,rot3r3_hat(Gi,d));
