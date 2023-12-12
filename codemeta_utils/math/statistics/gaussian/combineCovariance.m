function SigmaComb=combineCovariance(Sigma,varargin)
NSigma=size(Sigma,3);
w=ones(NSigma,1);

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'weights'
            ivarargin=ivarargin+1;
            w=varargin{ivarargin};
        case 'weightsaverage'
            ivarargin=ivarargin+1;
            w=varargin{ivarargin};
            w=w/sum(w);
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

for iSigma=1:NSigma
    Sigma(:,:,iSigma)=w(iSigma)*inv(Sigma(:,:,iSigma));
end
SigmaComb=inv(sum(Sigma,3));
