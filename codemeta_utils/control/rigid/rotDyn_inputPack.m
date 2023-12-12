%Pack the derivative of the state into a vector
%function dx=rotDyn_inputPack(dw,varargin)
%Pack the derivative of the state into a vector. By default, only the
%derivative of the velocity is packed. However, if the current state is
%provided, the packed state is extended to use it. In this way, traditional
%integrators such as ode45 can also be used.
%Optional arguments
%   'state',R,w     extend dx with the derivative of R (\dot{R}=R*hat3(w))
function dx=rotDyn_inputPack(dw,varargin)
flagExtendedVector=false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'state'
            ivarargin=ivarargin+1;
            R=varargin{ivarargin};
            ivarargin=ivarargin+1;
            w=varargin{ivarargin};
            flagExtendedVector=true;
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end
if ~flagExtendedVector
    dx=dw;
else
    dx=rotDyn_statePack(multiprod(R,hat3(w)),dw);
end
