%Second order dynamical model for a rigid body rotation
%function dx=rotDyn_model(x,varargin)
%Inputs
%   x       [6 x 1] current state x (packed)
%Outputs
%   dx      derivative of the state, \dot{x} (packed)
%Optional inputs
%   'noGyroscopic'      remove gyroscopic term (i.e., -hat(w)*J*w)
%   'inertiaMatrix',J   [3 x 3] moment of inertia tensor
%   'torque', Gamma     [3 x 1] input torque (in body's frame)
function dx=rotDyn_model(x,varargin)
flagGyroscopic=true;
Gamma=zeros(3,1);
J=eye(3);

%optional parameters
ivarargin=1;
len_varargin = length(varargin);
while ivarargin<=len_varargin
    switch lower(varargin{ivarargin})
        case 'nogyroscopic'
            flagGyroscopic=false;
        case 'torque'
            ivarargin=ivarargin+1;
            Gamma=varargin{ivarargin};
        case 'inertiamatrix'
            ivarargin=ivarargin+1;
            J=varargin{ivarargin};
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 <= len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end
    end
    ivarargin=ivarargin+1;
end

[R,w]=rotDyn_stateUnpack(x);
predw=Gamma;
if flagGyroscopic
    predw=predw-rotDyn_gyroscopicTerm(w,'inertiaMatrix',J);
end
dw=J\predw;
dx=rotDyn_inputPack(dw,'state',R,w);

% % Uncomment to limit speed
% MAX_SPEED = 0.5;
% newW = dx(10:12);
% if norm(newW) > MAX_SPEED
%     newW = MAX_SPEED*newW/norm(newW);
% end
% dx(10:12)=newW;
