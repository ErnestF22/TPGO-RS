function dx=rigidDyn_controlPDRT(x,varargin)
RReference=eye(3);
TReference=zeros(3,1);
J=eye(3);
m=1;
gainRotationError=1;
gainRotationVelocityError=1;
gainTranslationError=1;
gainTranslationVelocityError=1;
flagExtendedVector=false;
optsRotationControl={};
optsTranslationControl={};
disturbance=zeros(18,size(x,2));

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
        case 'gain'
            ivarargin=ivarargin+1;
            gainRotationError=varargin{ivarargin};
            gainRotationVelocityError=varargin{ivarargin};
            gainTranslationError=varargin{ivarargin};
            gainTranslationVelocityError=varargin{ivarargin};
        case 'gainrotation'
            ivarargin=ivarargin+1;
            gainRotationError=varargin{ivarargin};
            gainRotationVelocityError=varargin{ivarargin};
        case 'gaintranslation'
            ivarargin=ivarargin+1;
            gainTranslationError=varargin{ivarargin};
            gainTranslationVelocityError=varargin{ivarargin};            
        case 'gainrotationerror'
            ivarargin=ivarargin+1;
            gainRotationError=varargin{ivarargin};
        case 'gainrotationvelocityerror'
            ivarargin=ivarargin+1;
            gainRotationVelocityError=varargin{ivarargin};
        case 'gaintranslationerror'
            ivarargin=ivarargin+1;
            gainTranslationError=varargin{ivarargin};
        case 'gaintranslationvelocityerror'
            ivarargin=ivarargin+1;
            gainTranslationVelocityError=varargin{ivarargin};
        case 'flagextendedvector'
            ivarargin=ivarargin+1;
            flagExtendedVector=varargin{ivarargin};
        case 'disturbance'
            ivarargin=ivarargin+1;
            disturbance=disturbance+varargin{ivarargin};
        case 'optsrotationcontrol'
            ivarargin=ivarargin+1;
            optsRotationControl=[optsRotationControl varargin{ivarargin}];
        case 'optstranslationcontrol'
            ivarargin=ivarargin+1;
            optsTranslationControl=[optsTranslationControl varargin{ivarargin}];
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

dw=rotDyn_controlPD(x,'inertiaMatrix',J,'RReference',RReference,...
    'gainRotationError',gainRotationError,'gainVelocityError',gainRotationVelocityError,...
    'disturbance',disturbance,optsRotationControl{:});
dv=realDyn_controlPD(x, 'mass',m,'TReference',TReference,...
    'gainTranslationError',gainTranslationError,'gainVelocityError',gainTranslationVelocityError,...
    'disturbance',disturbance,optsTranslationControl{:});

if ~flagExtendedVector
    dx=rigidDyn_inputPack(dw,dv);
else
    [R,w,~,v]=rigidDyn_stateUnpackRT(x);
    dx=rigidDyn_inputPack(dw,dv,'state',R,w,v);
end    