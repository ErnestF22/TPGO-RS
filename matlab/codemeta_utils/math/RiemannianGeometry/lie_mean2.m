function [yMean,errors]=lie_mean2(y,lie_funs,varargin)
flagUseWeights=false;
flagNormalizeWeights=true;

itersMax=3000;

%assume data is stacked in the third dimension of the array y
N=size(y,3);
yMean=y(:,:,1);

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'init'
            ivarargin=ivarargin+1;
            yMean=varargin{ivarargin};
        case 'weights'
            ivarargin=ivarargin+1;
            w=varargin{ivarargin};
            flagUseWeights=true;
        case 'normalizeweights'
            ivarargin=ivarargin+1;
            flagNormalizeWeights=varargin{ivarargin};
        case 'init'
            ivarargin=ivarargin+1;
            yMean=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

lie_log=lie_funs.log;
lie_dist=lie_funs.dist;

if flagUseWeights
    W=repmat(shiftdim(shiftdim(w),-2), size(y,1), size(y,2));
    if(flagNormalizeWeights)
        W=W/sum(w);
        w=w/sum(w);
    end
    fun=@costWeight;
    der=@costDerWeight;
else
    fun=@costNoWeight;
    der=@costDerNoWeight;
end

[yMean,errors]=lie_minimizeArmijo(lie_funs,fun,der,yMean,'itersMax',itersMax,'method','conjgrad','armijoOpts',{'alpha',1,'beta',0.8,'sigma',1e-2},'retractionsOpts',{'polar'});
    
    function f=costWeight(ymean)
    f=sum((lie_dist(ymean,y).^2).*shiftdim(w))/2;
    end
    
    function d=costDerWeight(yMean)
    hy=lie_log(yMean,y);
    d=-sum(hy.*W,3);
    end
    
    function f=costNoWeight(ymean)
    f=sum(lie_dist(ymean,y).^2)/(2*N);
    end
    
    function d=costDerNoWeight(yMean)
    d=-sum(lie_log(yMean,y),3)/N;
    end
end
