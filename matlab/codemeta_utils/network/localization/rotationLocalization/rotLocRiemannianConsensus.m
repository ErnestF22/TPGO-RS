function R=rotLocRiemannianConsensus(E,RRel,varargin)
NNodes=max(E(:));
epsilon=1/(2*edges2maxDegree(E));
maxIt=30;
funs=rot3_almostGlobal_functions('type','tron','b',2);
flagCollect=false;
flagShowProgress=false;
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
        case 'epsilon'
            ivarargin=ivarargin+1;
            epsilon=varargin{ivarargin};
        case 'epsilonfactor'
            ivarargin=ivarargin+1;
            epsilon=epsilon*varargin{ivarargin};
        case 'funs'
            ivarargin=ivarargin+1;
            funs=varargin{ivarargin};
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

flagTolUpdate=tolUpdate~=0;

sz=size(RInit);
if flagCollect
    R=zeros([sz maxIt]);
else
    R=zeros(sz);
end
R(:,:,:,1)=RInit;
p=1;
if flagTolUpdate
    normVUpdatePrev=Inf;
end
if flagShowProgress
    w=getTextWaitBar(maxIt);
    w(0)
end
for it=1:maxIt
    gradc=rotLocRiemannianCost_grad(E,R(:,:,:,p),RRel,funs);
    RNew=rot_exp(R(:,:,:,p),-epsilon*gradc);
    if flagCollect
        p=p+1;
    end
    R(:,:,:,p)=RNew;
    if flagTolUpdate
        normVUpdate=sqrt(sum(rot_metric(R(:,:,:,p),gradc,gradc)))/NNodes;
        if (normVUpdatePrev-normVUpdate)<tolUpdate
            break
        else
            normVUpdatePrev=normVUpdate;
        end
   end
   if flagShowProgress
        w(it)
    end
end
R=R(:,:,:,1:p);
