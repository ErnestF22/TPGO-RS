function [muComb,SigmaComb]=combineMeanAndCovariance(mu,Sigma,varargin)
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

SigmaComb=combineSigma(Sigma,'weights',w);
NMu=size(mu,2);
for iMu=1:NMu
    mu(:,iMu)=w(iMu)*Sigma(:,:,iMu)\mu(:,iMu);
end
muComb=SigmaComb*sum(mu,2);
