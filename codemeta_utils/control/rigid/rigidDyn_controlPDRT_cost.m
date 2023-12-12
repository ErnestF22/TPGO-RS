function [phi,gradphi]=rigidDyn_controlPDRT_cost(R,w,T,v,varargin)
RReference=eye(3);
TReference=zeros(3,1);
J=eye(3);
m=1;
gainRotationError=1;
gainTranslationError=1;
methodGradient='packed';

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'reference'
            ivarargin=ivarargin+1;
            RReference=varargin{ivarargin};
            ivarargin=ivarargin+1;
            TReference=varargin{ivarargin};
        case 'inertiamatrix'
            ivarargin=ivarargin+1;
            J=varargin{ivarargin};
        case 'mass'
            ivarargin=ivarargin+1;
            m=varargin{ivarargin};
        case 'gainrotationerror'
            ivarargin=ivarargin+1;
            gainRotationError=varargin{ivarargin};
        case 'gaintranslationerror'
            ivarargin=ivarargin+1;
            gainTranslationError=varargin{ivarargin};
        case 'methodgradient'
            ivarargin=ivarargin+1;
            methodGradient=lower(varargin{ivarargin});
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

[phiRotation,gradRotation]=rotDyn_controlPD_cost(R,w,'inertiaMatrix',J,'RReference',RReference,'gainRotationError',gainRotationError,'methodGradient',methodGradient);
[phiTranslation,gradTranslation]=realDyn_controlPD_cost(T,v,'mass',m,'TReference',TReference,'gainTranslationError',gainTranslationError);

phi=phiRotation+phiTranslation;
gradphi=[gradRotation; gradTranslation];
