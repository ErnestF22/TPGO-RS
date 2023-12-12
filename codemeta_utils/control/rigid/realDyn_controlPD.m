function dv=realDyn_controlPD(x,varargin)
flagExtendedVector=false;
TReference=zeros(3,1);
m=1;
gainTranslationError=1;
gainVelocityError=1;
disturbance=zeros(6,size(x,2));

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
        case 'gain'
            ivarargin=ivarargin+1;
            gainTranslationError=varargin{ivarargin};
            gainVelocityError=varargin{ivarargin};
        case 'gaintranslationerror'
            ivarargin=ivarargin+1;
            gainTranslationError=varargin{ivarargin};
        case 'gainvelocityerror'
            ivarargin=ivarargin+1;
            gainVelocityError=varargin{ivarargin};
        case 'disturbance'
            ivarargin=ivarargin+1;
            disturbance=disturbance+varargin{ivarargin}(end-5:end,:);
        case 'flagextendedvector'
            ivarargin=ivarargin+1;
            flagExtendedVector=varargin{ivarargin};
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NT=size(x,2);
[T,v]=realDyn_stateUnpack(x);
if size(TReference,2)==1 && NT>1
    TReference=repmat(TReference,1,NT);
end
uT=gainTranslationError*(TReference-T);
uv=-gainVelocityError*v;
dvDisturbance=realDyn_inputUnpack(disturbance);
dv=uT+uv-dvDisturbance;
if flagExtendedVector
    dv=realDyn_inputPack(dv,'state',v);
end
