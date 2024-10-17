function R=rotLocWeiszfeldConsensus(E,RRel,varargin)
NNodes=max(E(:));
maxIt=30;
flagCollect=false;
flagShowProgress=false;
tolCost=0;
tolUpdate=0;

RInit=rot_randn(eye(3),[],NNodes);

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'rinit'
            ivarargin=ivarargin+1;
            RInit=varargin{ivarargin};
        case 'nnodes'
            ivarargin=ivarargin+1;
            NNodes=varargin{ivarargin};
            RInit=rot_randn(3,3,NNodes);
        case 'maxit'
            ivarargin=ivarargin+1;
            maxIt=varargin{ivarargin};
        case 'collect'
            flagCollect=true;
        case 'tolcost'
            ivarargin=ivarargin+1;
            tolCost=varargin{ivarargin};
        case 'tolupdate'
            ivarargin=ivarargin+1;
            tolUpdate=varargin{ivarargin};
        case 'showprogress'
            flagShowProgress=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

flagTolCost=tolCost~=0;
flagTolUpdate=tolUpdate~=0;

sz=size(RInit);
if flagCollect
    R=zeros([sz maxIt]);
else
    R=zeros(sz);
end
R(:,:,:,1)=RInit;
p=1;
if flagTolCost
    cPrev=rotLocWeiszfeldCost(E,RInit,RRel);
end
if flagTolUpdate
    normVUpdatePrev=Inf;
end
if flagShowProgress
    w=getTextWaitBar(maxIt);
    w(0)
end
for it=1:maxIt
    if flagCollect
        p=p+1;
        R(:,:,:,p)=R(:,:,:,p-1);
    end
    if flagTolUpdate
        normVUpdate=0;
    end
    
    %update of each node in sequence
    for iNode=1:NNodes
        v=tangentVectors(E,RRel,R(:,:,:,p),iNode);
        [vNorm,normV]=cnormalize(v);
        flagNormV=normV~=0;
        vUpdate=sum(vNorm(:,flagNormV),2)/sum(1./normV(flagNormV));
        R(:,:,iNode,p)=R(:,:,iNode,p)*rot(vUpdate);
        if flagTolUpdate
            normVUpdate=normVUpdate+sum(vUpdate.^2);
        end
    end
    if flagTolUpdate
        normVUpdate=sqrt(normVUpdate)/NNodes;
        if (normVUpdatePrev-normVUpdate)<tolUpdate
            break
        else
            normVUpdatePrev=normVUpdate;
        end
    end
    if flagTolCost
        c=rotLocWeiszfeldCost(E,R(:,:,:,p),RRel);
        if (cPrev-c)<tolCost
            break
        else
            cPrev=c;
        end
    end
    
    if flagShowProgress
        w(it)
    end
end
R=R(:,:,:,1:p);

function v=tangentVectors(E,RRel,R,iNode)
flagij=E(:,1)==iNode;
jNodeij=E(flagij,2);
Rij=multiprod(multiprod(R(:,:,iNode)',R(:,:,jNodeij)),invR(RRel(:,:,flagij)));
flagji=E(:,2)==iNode;
jNodeji=E(flagji,1);
Rji=multiprod(multiprod(R(:,:,iNode)',R(:,:,jNodeji)),RRel(:,:,flagji));
v=[logrot(Rij) logrot(Rji)];

%rot_dist(Ri'*Rj*Rij',I)