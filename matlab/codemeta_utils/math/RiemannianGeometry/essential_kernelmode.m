%Estimate mode of a sampled distribution of QREMs by using kernel density estimation
%function [QMax,d,f]=essential_kernelmode(Q,kernelSigma,varargin)
%Input
%   Q               [6 x 3 x NPoints] array of QREMs
%   kernelSigma     sigma to use for the kernel
%Optional arguments
%   'methodDistance', method
%       'full'      compute all the distances
%       'sparse'    try to skip computation of some of the distances using
%                   the triangualr inequality
function [QMax,d,f,idxMax]=essential_kernelmode(Q,kernelSigma,varargin)
methodDistance='full';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'methoddistance'
            ivarargin=ivarargin+1;
            methodDistance=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

N=size(Q,3);
switch methodDistance
    case 'full'
        d=zeros(N);
        for iN1=2:N
            for iN2=1:iN1-1
                d(iN1,iN2)=essential_dist(Q(:,:,iN1),Q(:,:,iN2));
                d(iN2,iN1)=d(iN1,iN2);
            end
        end
        f=sum(exp(-d.^2/(2*kernelSigma)),2);
    case 'sparse'
        dThreshold=2*kernelSigma;
        d=NaN(N);
        d(1,1)=0;
        for iN1=2:N
            d(iN1,iN1)=0;
            d(iN1,1)=essential_dist(Q(:,:,iN1),Q(:,:,1));
            d(1,iN1)=d(iN1,1);
            for iN2=1:iN1-1
                flagNeedToCompute=true;
                for iN3=1:iN2
                    if ~isnan(d(iN1,iN3)) && ~isnan(d(iN2,iN3)) ...
                        && abs(d(iN1,iN3)-d(iN2,iN3))>dThreshold
                        flagNeedToCompute=false;
                        break
                    end
                end
                if flagNeedToCompute
                    d(iN1,iN2)=essential_dist(Q(:,:,iN1),Q(:,:,iN2));
                    d(iN2,iN1)=d(iN1,iN2);
                end
            end
        end
        f=zeros(N,1);
        for iN=1:N
            flagValid=~isnan(d(iN,:));
            f(iN)=sum(exp(-d(iN,flagValid).^2/(2*kernelSigma)));
        end
    otherwise
        error('Unknown methodDistance')
end

[~,idxMax]=max(f);
QMax=Q(:,:,idxMax);
