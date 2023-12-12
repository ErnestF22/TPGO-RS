function dx=rigidDyn_model(x,varargin)
Nx=size(x,2);
Gamma=zeros(3,1);
F=zeros(3,Nx);
J=eye(3);
m=1;
disturbanceGamma=zeros(3,Nx);
disturbanceF=zeros(3,Nx);
flagDisplayInfo=false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'inertiamatrix'
            ivarargin=ivarargin+1;
            J=varargin{ivarargin};
        case 'mass'
            ivarargin=ivarargin+1;
            m=varargin{ivarargin};
        case 'torque'
            ivarargin=ivarargin+1;
            Gamma=varargin{ivarargin};
        case 'force'
            ivarargin=ivarargin+1;
            F=varargin{ivarargin};
        case 'input'
            ivarargin=ivarargin+1;
            [Gamma,F]=rigidDyn_inputUnpack(varargin{ivarargin});
        case 'displayinfo'
            flagDisplayInfo=true;
        case 'disturbance'
            ivarargin=ivarargin+1;
            [disturbanceGamma,disturbanceF]=rigidDyn_inputUnpack(varargin{ivarargin});
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NPoses=size(x,2);
dx=zeros(18,NPoses);
dx(1:12,:)=rotDyn_model(x(1:12,:),'torque',Gamma+disturbanceGamma,'inertiaMatrix',J);
dx(13:18,:)=realDyn_model(x(13:18,:),'force',F+disturbanceF,'mass',m);
if flagDisplayInfo
    disp('[F disturbanceF F+disturbanceF]')
    disp([F disturbanceF F+disturbanceF])
end

