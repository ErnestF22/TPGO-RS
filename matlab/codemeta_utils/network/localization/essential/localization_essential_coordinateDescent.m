function [t_node,output]=localization_essential_coordinateDescent(t_node,varargin)
flagCollectErrors=false;
flagDisplayIt=false;
flagProgressBar=false;
maxIt=1000;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'displayit'
            flagDisplayIt=true;
        case 'progressbar'
            flagProgressBar=true;
        case 'maxit'
            ivarargin=ivarargin+1;
            maxIt=varargin{ivarargin};
        case 'geterrors'
            flagCollectErrors=true;
        case 'optslieminimize'
            ivarargin=ivarargin+1;
            optsLieMinimize=[optsLieMinimize varargin{ivarargin}];
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if flagDisplayIt
    flagProgressBar=false;
end

E=testNetworkGetEdges(t_node);
A=t_node.A;
N=testNetworkGetNumberOfNodes(t_node);
NEdges=size(E,1);
[~,Qij]=testNetworkGetRelativeEssential(t_node);
[Ri0,Ti0]=testNetworkGetRotTransl(t_node);
GiInit=RT2G(Ri0,Ti0);

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
    c=zeros(maxIt*N+1,1);
    t=zeros(maxIt*N+1,1);
    c(1)=essentialCostNetwork(G2R(Gi),G2T(Gi),Qij,E);
    cnt=2;
end
if flagProgressBar
    w=getTextWaitBar(maxIt);
end
tStart=cputime;
for it=1:maxIt
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
tFinal=cputime-tStart;
fprintf('cputime: %.2f\n',tFinal);

output.cFinal=essentialCostNetwork(G2R(Gi),G2T(Gi),Qij,E);
output.tFinal=tFinal;
if flagCollectErrors
    output.c=c;
    output.t=t;
end
t_node.gi=Gi;

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
