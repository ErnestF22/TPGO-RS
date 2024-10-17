%Implement the Weiszfeld algorithm
%function [y,output]=lie_minimize_tangentAverage(lf,y,yi,varargin)
%Optional inputs
%   'maxit'     maximum number of iterations
%   'displayIt' display info on each iteration
%   'norm'      norm used for the average in the tangent space
function [y,output]=lie_minimize_tangentAverage(lf,y,yi,varargin)
flagGetErrors=false;
flagDisplayIt=false;
LNorm=2;
thresholdDist=1e-12;
maxIt=10;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'displayit'
            flagDisplayIt=true;
        case 'maxit'
            ivarargin=ivarargin+1;
            maxIt=varargin{ivarargin};
        case 'norm'
            ivarargin=ivarargin+1;
            LNorm=varargin{ivarargin};
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if flagDisplayIt
    flagGetErrors=flagDisplayIt;
end

if flagGetErrors
    c=NaN(maxIt+1,1);
    t=NaN(maxIt+1,1);
    fCost=@(y) lie_minimize_tangentAverage_cost(lf,y,yi,'norm',LNorm);
        
    c(1)=fCost(y);
    t(1)=0;
    tStart=cputime;
end

Nyi=size(yi,3);
if Nyi==1
    y=yi;
else
    if LNorm~=2
        w=zeros(1,1,Nyi);
        dimY=size(y);
    end
    if flagDisplayIt
        fprintf('it: %4d t:%5.2f c:%.4e\n',0,0,c(1));
    end
    for it=1:maxIt
        vi=lf.log(y,yi);
        if LNorm==2
            d=mean(vi,3);
        else
            for iv=1:Nyi
                w(iv)=sqrt(lf.metric(y,vi(:,:,iv),vi(:,:,iv)));
            end
            flagValid=w>thresholdDist;
            w=w.^(LNorm-2);
            wValid=w(flagValid);
            viValid=vi(:,:,flagValid);
            d=sum(viValid.*repmat(wValid,dimY),3)/sum(w);
        end
        y=lf.exp(y,d);
        if flagGetErrors
            c(it+1)=fCost(y);
            t(it+1)=cputime-tStart;
        end
        if flagDisplayIt
            fprintf('it: %4d t:%5.2f c:%.4e\n',it,t(it+1),c(it+1));
        end
    end
end

if flagGetErrors
    output.c=c;
    output.t=t;
    output.fCost=fCost;
end
