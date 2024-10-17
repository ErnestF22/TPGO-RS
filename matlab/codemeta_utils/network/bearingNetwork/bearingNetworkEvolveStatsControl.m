function output=bearingNetworkEvolveStatsControl(ode,t,x,varargin)
Nt=length(t);
dx=zeros(size(x));
sz=[size(x,1) size(x,2)];
flagGetNorm=false;
flagGetTotalNorm=false;
flagGetCumulativeNormSq=false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'norm'
            flagGetNorm=true;
        case 'totalnorm'
            flagGetTotalNorm=true;
        case 'cumulativenormsq'
            flagGetTotalNorm=true;
            flagGetCumulativeNormSq=true;
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if flagGetNorm
    n=zeros(Nt,sz(2));
end
if flagGetTotalNorm
    totalNorm=zeros(Nt,1);
end
for it=1:Nt
    dx(:,:,it)=reshape(ode(t(it),reshape(x(:,:,it),[],1)),sz);
    if flagGetNorm
        n(it,:)=cnorm(dx(:,:,it));
    end
    if flagGetTotalNorm
        totalNorm(it)=norm(reshape(dx(:,:,it),[],1));
    end
end
if flagGetCumulativeNormSq
    totalNormAverage=filter([1/2 1/2],1,totalNorm);
    totalNormIntervals=[1;diff(t)];
    cumulativeNormSq=cumsum(totalNormAverage.*totalNormIntervals);
end
output.dx=dx;

if flagGetNorm
    output.norm=n;
end
if flagGetTotalNorm
    output.totalNorm=totalNorm;
end
if flagGetCumulativeNormSq
    output.cumulativeNormSq=cumulativeNormSq;
end
