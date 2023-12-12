function dx=quadrotor_model(x,varargin)
Gamma=zeros(3,1);
J=eye(3);
m=1;
disturbance=zeros(size(x));

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
        case 'thrust'
            ivarargin=ivarargin+1;
            thrust=varargin{ivarargin};
        case 'input'
            ivarargin=ivarargin+1;
            [Gamma,thrust]=quadrotor_inputUnpack(varargin{ivarargin});
        case 'disturbance'
            ivarargin=ivarargin+1;
            disturbance=disturbance+varargin{ivarargin};
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

R=rotDyn_stateUnpack(x);
r3=squeeze(R(:,3,:));
F=r3.*repmat(thrust,3,1);
dx=rigidDyn_model(x,'inertiaMatrix',J,'mass',m,...
    'torque',Gamma,'force',F,...
    'disturbance',disturbance);
