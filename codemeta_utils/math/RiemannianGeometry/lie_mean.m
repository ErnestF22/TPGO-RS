function yMean=lie_mean(y,lie_funs,varargin)
maxIt=100;      %maximum number of iterations
minIt=1;        %minimum number of iterations before stopping
flagshow=false;
flagShowIt=false;   %show iterations
flagNormalizeWeights=true;
flagLineSearch=true;

%assume data is stacked in the third dimension of the array y
N=size(y,3);
yMean=y(:,:,1);
w=ones(N,1);    %default weights

tolMean=1e-14;
tolCost=1e-14;

epsilon=1;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'maxit'
            ivarargin=ivarargin+1;
            maxIt=varargin{ivarargin};
        case 'init'
            ivarargin=ivarargin+1;
            yMean=varargin{ivarargin};
        case 'weights'
            ivarargin=ivarargin+1;
            w=varargin{ivarargin};
        case 'normalizeweights'
            ivarargin=ivarargin+1;
            flagNormalizeWeights=varargin{ivarargin};
        case 'minit'
            ivarargin=ivarargin+1;
            minIt=varargin{ivarargin};
        case 'showcost'
            flagshow=true;
        case 'stepsize'
            ivarargin=ivarargin+1;
            epsilon=varargin{ivarargin};
        case 'showit'
            flagShowIt=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end
if(minIt>maxIt-1)
    maxIt=minIt+1;
end

lie_exp=lie_funs.exp;
lie_log=lie_funs.log;
lie_dist=lie_funs.dist;

tolMean=tolMean*epsilon;
tolCost=tolCost*epsilon;

W=repmat(shiftdim(shiftdim(w),-2), size(y,1), size(y,2));
if(flagNormalizeWeights)
    W=W/sum(w);
    w=w/sum(w);
end

fPrev=Inf;
if flagshow
    allf(1)=cost(yMean,y,w);
end

for it=1:maxIt
    if flagShowIt
        fprintf('It: %d\n', it)
    end
    yPrev=yMean;
    hy=lie_log(yMean,y);
    
%     hycent=hy;
%     for(iN=1:N)
%         hycent(:,:,iN)=yPrev'*hy(:,:,iN);
%     end
    
    geodCost=@(t) cost(lie_exp(yMean,t*sum(hy.*W,3)),y,w);
    if(flagLineSearch)
        topt=fminbnd(geodCost,0,1);
    else
        topt=epsilon;
    end
    yMean=lie_exp(yMean,topt*sum(hy.*W,3));
    
    fCurrent=cost(yMean,y,w);
%    fPrev-fCurrent
    if(fCurrent-fPrev>1e-7)
        warning(['Cost went up at iteration ' num2str(it)...
            ' with gradient magnitude ' num2str(sqrt(lie_funs.metric(yPrev,topt*sum(hy.*W,3),topt*sum(hy.*W,3)))) '! Stopping here'])
        break;
    end
    if it>minIt && (fPrev-fCurrent)<tolCost && lie_dist(yPrev,yMean)<tolMean
        %disp(['Terminated after ' num2str(it) ' iteration with cost ' num2str(fCurrent,'%.12e')])
        break;
    end
    fPrev=fCurrent;
    if flagshow
        allf(it+1)=cost(yMean,y,w);
    end
    
end

if flagshow
    semilogy(allf)
end

    function f=cost(ymean,y,w)
    f=sum((shiftdim(lie_dist(ymean,y)).^2).*w);
    end
    
    function d=costDer(ymean)
    hy=lie_log(yMean,y);
    d=sum(hy.*W,3);
    end
end
