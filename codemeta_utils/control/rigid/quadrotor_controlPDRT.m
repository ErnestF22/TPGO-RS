function dx=quadrotor_controlPDRT(x,varargin)
RReference=eye(3);
TReference=zeros(3,1);
J=eye(3);
m=1;
gainRotationError=1;
gainRotationVelocityError=1;
gainTranslationError=1;
gainTranslationVelocityError=1;
optsRotationControl={};
optsTranslationControl={};
disturbance=zeros(18,size(x,2));
flagDisplayInfo=false;

tolDirection=1e-14;

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
        case 'disturbance'
            ivarargin=ivarargin+1;
            disturbance=disturbance+varargin{ivarargin};
        case 'optsrotationcontrol'
            ivarargin=ivarargin+1;
            optsRotationControl=[optsRotationControl varargin{ivarargin}];
        case 'optstranslationcontrol'
            ivarargin=ivarargin+1;
            optsTranslationControl=[optsTranslationControl varargin{ivarargin}];
        case 'displayinfo'
            flagDisplayInfo=true;
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

dv=realDyn_controlPD(x, 'mass',m,'TReference',TReference,...
    'gainTranslationError',gainTranslationError,'gainVelocityError',gainTranslationVelocityError,...
    'disturbance',disturbance,optsTranslationControl{:});
R=rotDyn_stateUnpack(x);
thrust=sum(dv.*squeeze(R(:,3,:)));
signThrust=sign(thrust);
signThrust(signThrust==0)=1;

Nx=size(x,2);
RReferenceZMin=zeros(3,3,Nx);
for ix=1:Nx
    if norm(dv(:,ix))>tolDirection
        %obtain a rotation matrix with dv on the last column
        RReferenceZ=-fliplr(orthCompleteBasis(-signThrust(ix)*dv(:,ix)));
        %find the rotation matrix with dv on the last column which is closest to
        %RReference
        RReferenceZMin(:,:,ix)=rot3_projectOnGeodesicAlgebraic(RReference(:,:,ix),RReferenceZ,rot_hat(RReferenceZ,[0;0;1]));
    else
        RReferenceZMin(:,:,ix)=RReference(:,:,ix);
    end

    if flagDisplayInfo
        dvEff=thrust(ix)*R(:,3,ix);
        disp('subspace(dv,dvEff)')
        disp(subspace(dv(:,ix),dvEff)*180/pi)
        disp('[dv dvEff dv-dvEff]')
        disp([dv(:,ix) dvEff dv(:,ix)-dvEff])
    end
end

dw=rotDyn_controlPD(x,'inertiaMatrix',J,'RReference',RReferenceZMin,...
    'gainRotationError',gainRotationError,'gainVelocityError',gainRotationVelocityError,...
    'disturbance',disturbance,optsRotationControl{:});

dx=quadrotor_inputPack(dw,thrust);
