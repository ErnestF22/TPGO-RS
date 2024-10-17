function [R,Q]=rotLocFrobConsensus(E,RRel,varargin)
NNodes=max(E(:));
epsilon=1/(2*edges2maxDegree(E));
maxIt=30;
idxLeader=1;
flagCollect=false;
flagShowProgress=false;
methodFixLeader='identity';
d=size(RRel,1);

RInit=rot_randn(eye(d),[],NNodes);

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'rinit'
            ivarargin=ivarargin+1;
            RInit=varargin{ivarargin};
        case 'methodfixleader'
            ivarargin=ivarargin+1;
            methodFixLeader=lower(varargin{ivarargin});
        case 'nnodes'
            ivarargin=ivarargin+1;
            NNodes=varargin{ivarargin};
            RInit=rot_randn(3,3,NNodes);
        case 'maxit'
            ivarargin=ivarargin+1;
            maxIt=varargin{ivarargin};
        case 'collect'
            flagCollect=true;
        case 'idxleader'
            ivarargin=ivarargin+1;
            idxLeader=varargin{ivarargin};
        case 'epsilon'
            ivarargin=ivarargin+1;
            epsilon=varargin{ivarargin};
        case 'epsilonfactor'
            ivarargin=ivarargin+1;
            epsilon=epsilon*varargin{ivarargin};
        case 'showprogress'
            flagShowProgress=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

sz=size(RInit);
if flagCollect
    Q=zeros([sz maxIt]);
else
    Q=zeros(sz);
end
Q(:,:,:,1)=multiprod(RInit(:,:,1),RInit);
Q(:,:,idxLeader,1)=eye(d);
R=Q;
p=1;
if flagShowProgress
    w=getTextWaitBar(maxIt);
    w(0)
end
for it=1:maxIt
    gradc=rotLocFrobCost_grad(E,Q(:,:,:,p),RRel);
    switch methodFixLeader
        case 'identity'
            gradc(:,:,idxLeader)=0;
        case 'triangular'
            for iRow=1:d
                for iCol=iRow+1:d
                    gradc(iRow,iCol,idxLeader)=0;
                end
                gradc(iRow,iRow)=1;
            end
        otherwise
            error('Method to fix leader not recognized')
    end
    QNew=Q(:,:,:,p)-epsilon*gradc;
    if flagCollect
        p=p+1;
    end
    Q(:,:,:,p)=QNew;
    R(:,:,:,p)=rot_proj(Q(:,:,:,p));
    if flagShowProgress
        w(it)
    end
end
