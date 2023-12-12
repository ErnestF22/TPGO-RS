function dx=rigidDyn_disturbance_gravity(x,varargin)
m=1;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'mass'
            ivarargin=ivarargin+1;
            m=varargin{ivarargin};
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end
R=rotDyn_stateUnpack(x);
Nx=size(x,2);
dx=rigidDyn_statePackRT(zeros(3,3,Nx),zeros(3,Nx),zeros(3,Nx),m*repmat(gravityVector(),1,Nx));
