%function ymedian=lie_median(lf,y,varargin)
%Compute the median of the points y using the Weiszfeld algorithm. 

function ymedian=lie_median(lf,y,varargin)
Nit=100;
tolMeas=1e-12;
flagShowCost=false;

y0=y(:,:,1);

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'nit'
            ivarargin=ivarargin+1;
            Nit=varargin{ivarargin};
        case 'showcost'
            flagShowCost=true;
        case 'init'
            ivarargin=ivarargin+1;
            y0=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

N=size(y,3);
ymedian=y0;
dNum=0;
dDen=0;

if flagShowCost
    costEval=zeros(1,Nit+1);
    cost=@(ymedian) sum(lf.dist(ymedian,y));
    costEval(1)=cost(ymedian);
end

for it=1:Nit
    v=lf.log(ymedian,y);
    vnorm=sqrt(lf.metric(ymedian,v,v));
    
    for iMeasurement=find(vnorm>tolMeas)
        dNum=dNum+v(:,:,iMeasurement)/vnorm(iMeasurement);
        dDen=dDen+1/vnorm(iMeasurement);
    end
    
    ymedian=lf.exp(ymedian,dNum/dDen);
    if flagShowCost
        costEval(it+1)=cost(ymedian);
    end
end

if flagShowCost
    figure
    plot(0:Nit,costEval)
end
