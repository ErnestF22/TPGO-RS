function dx=realDyn_disturbance_gravity(x,varargin)
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

Nx=size(x,2);
dx=[zeros(3,Nx); m*gravityVector()*ones(1,Nx)];
