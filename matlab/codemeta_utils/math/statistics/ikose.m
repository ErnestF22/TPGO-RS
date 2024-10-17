%Iterated k-th Order Scale Estimator
%function s=ikose(r,K,varargin)
%This function implements the IKOSE estimator proposed in
%Wang, Hanzi, Tat-Jun Chin, and David Suter.
%"Simultaneously fitting and segmenting multiple-structure data with outliers." 
%IEEE Transactions on Pattern Analysis and Machine Intelligence, 2012.

function s=ikose(r,K,varargin)
maxIt=1000;
KQuantile=0.1;
flagKIsQuantile=false;
flagShowStats=false;
factorSigmaInliers=2.5;
methodInlierCount='standard';

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'maxit'
            ivarargin=ivarargin+1;
            maxIt=varargin{ivarargin};
        case 'quantile'
            flagKIsQuantile=true;
        case 'showstats'
            flagShowStats=true;
        case 'factorsigmainliers'
            ivarargin=ivarargin+1;
            factorSigmaInliers=varargin{ivarargin};
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

n=length(r);
if flagKIsQuantile
    KQuantile=K;
end
if flagKIsQuantile || ~exist('K','var') || isempty(K)
    K=round(KQuantile*n);
end
if flagShowStats
    fprintf('K=%d, n=%d\n',K,n);
end
rSorted=sort(abs(r));
rK=rSorted(K);

nInliers=n;
for it=1:maxIt
    s=scaleEstimate(rK,K,nInliers);
    nInliersPrev=nInliers;
    switch methodInlierCount
        case 'standard'
            nInliers=sum(r<2.5*s);
        case 'cumulative'
            nInliers=round(sum(r<factorSigmaInliers*s)/absNormalCumulative(factorSigmaInliers));
    end
    if flagShowStats
        fprintf('%d: s=%.4f, nInliers=%d (was %d)\n',it,s,nInliers,nInliersPrev)
    end
    if nInliers==nInliersPrev
        break
    end
end

function s=scaleEstimate(rK,K,N)
s=rK/absNormalQuantile(K/N);

%Quantile distribution (inverse of the cumulative distribution) for the
%absolute of normally distributed variables
function x=absNormalQuantile(F)
x=erfcinv(1-F)*sqrt(2);

function F=absNormalCumulative(x)
F=1-erfc(x/sqrt(2));
