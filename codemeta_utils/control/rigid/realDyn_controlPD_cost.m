function [phi,gradphi]=realDyn_controlPD_cost(T,v,varargin)
TReference=zeros(3,1);
m=1;
flagComputeGrad=nargout>1;
gainTranslationError=1;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'treference'
            ivarargin=ivarargin+1;
            TReference=varargin{ivarargin};
        case 'mass'
            ivarargin=ivarargin+1;
            m=varargin{ivarargin};
        case 'gaintranslationerror'
            ivarargin=ivarargin+1;
            gainTranslationError=varargin{ivarargin};
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NT=size(T,2);
if size(TReference,2)==1 && NT>1
    TReference=repmat(TReference,[1 NT]);
end

eT=T-TReference;
phi=gainTranslationError*cnorm(eT).^2/2+m*cnorm(v).^2/2; ...
if flagComputeGrad
    gradphi=[gainTranslationError*eT;m*v];
end
