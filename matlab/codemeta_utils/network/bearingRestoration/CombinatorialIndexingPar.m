function [Erigid,ReqEdges,lambda2,p]=CombinatorialIndexingPar(E,x,membership,ReqEdges,varargin)
sigmaNoise=0.0;
%optional parameters
ivarargin=1;
flagOperationMode='BestL2';
%flagOperationMode='FirstOption';
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'flagoperationmode'
            ivarargin=ivarargin+1;
            flagOperationMode=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if length(unique(membership))==1
    Erigid=E;
    u=bearingCluster_getBearingsScalesFromE(x,E,'noisy',sigmaNoise);
    eigval=getEigVal(E,u);
    lambda2=eigval(2);
    p.ItStart=0;
    p.ItStop=0;
    return
end
if ReqEdges<=0
    if length(unique(membership))~=1
        disp('HOLY MOLY!');
    end
    Erigid=E;
    u=bearingCluster_getBearingsScalesFromE(x,E,'noisy',sigmaNoise);
    eigval=getEigVal(E,u);
    lambda2=eigval(2);
    Extra=[toc,lambda2];
    return
end

Nedges=size(E,1);
Nnodes=max(E(:));
Nalledges=Nnodes*(Nnodes-1)/2;
Naddedges=Nalledges-Nedges;

E=sort(E,2);
Ek=combnk(1:Nnodes,2);
Er=setdiff(Ek,E,'rows');

c=unique(membership);
Ecomp=[];
for idxComp=1:length(c)
    address=find(membership==c(idxComp));
    Comps=sort(unique([E(address,1);E(address,2)]));
    Ecomp=[Ecomp;combnk(Comps,2)];
end
Er=setdiff(Er,Ecomp,'rows');
Naddedges=size(Er,1);
%V=combnk(1:Naddedges,ReqEdges);
%v=[];
%Iterations=nchoosek(Naddedges,ReqEdges);
flagTerminate=false;
result=false;
lambda2=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(flagOperationMode,'BestL2')
    v=[];
    while ~result
        [v,flagTerminate]=combnkIndexGen(Naddedges,ReqEdges,v);
        while ~flagTerminate
            Eloop=[E;Er(v,:)];
            u=bearingCluster_getBearingsScalesFromE(x,Eloop,'noisy',sigmaNoise);
            membership=bearingCluster_clustering(Eloop,u,'flagfirstcall',true);
            c=unique(membership);
            if length(c)==1
                result=true;
                eigval=getEigVal(Eloop,u);
                l=eigval(2);
                if l>lambda2
                    lambda2=l;
                    Vgood=v;
                end
            end
            [v,flagTerminate]=combnkIndexGen(Naddedges,ReqEdges,v);
        end
        v=[];
       ReqEdges=ReqEdges+1;
    end
    Erigid=Er(Vgood,:);
    ReqEdges=ReqEdges-1;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
elseif strcmp(flagOperationMode,'FirstOption')
    v=[];
    while ~result
        [v,flagTerminate]=combnkIndexGen(Naddedges,ReqEdges,v);
        Eloop=[E;Er(v,:)];
        u=bearingCluster_getBearingsScalesFromE(x,Eloop,'noisy',sigmaNoise);
        membership=bearingCluster_clustering(Eloop,u,'flagfirstcall',true);
        c=unique(membership);
        if length(c)==1
            result=true;
            eigval=getEigVal(Eloop,u);
            lambda2=eigval(2);
            Vgood=v;
        end
    if flagTerminate
        v=[];
        ReqEdges=ReqEdges+1;
    end
    end
    Erigid=Er(Vgood,:);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(flagOperationMode,'Both')
    lambda2=zeros(nchoosek(Naddedges,ReqEdges),1);
    p = Par(nchoosek(Naddedges,ReqEdges));
    parfor Idx=1:nchoosek(Naddedges,ReqEdges)
        Par.tic;
        v=CombIndexGenIdx(Naddedges,ReqEdges,Idx);
        Eloop=[E;Er(v,:)];
        u=bearingCluster_getBearingsScalesFromE(x,Eloop,'noisy',sigmaNoise);
        membership=bearingCluster_clustering(Eloop,u,'flagfirstcall',true);
        c=unique(membership);
        if length(c)==1
            eigval=getEigVal(Eloop,u);
            l=eigval(2);
             p(Idx) = Par.toc;
            lambda2(Idx)=l;
        end
    end
    stop(p);
    [~,I]=max(lambda2);
    v=CombIndexGenIdx(Naddedges,ReqEdges,I);
    Erigid=Er(v,:);
else
    error('WRONG OPERATION MODE!');
end