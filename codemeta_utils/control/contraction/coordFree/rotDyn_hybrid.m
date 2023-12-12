function [dw, state, RRef_new, wreference] = rotDyn_hybrid(x,t,varargin)
% Hybrid controller on SO(3). There are three controllers
% 1) PD controller when inside contraction region
% 2) PD controller with reference trajectory when between 1) and 3)
% 3) Vel controller when high speed
% INPUTS:
%   x := the current state
% OUTPUTS:
%   dw := the updated dynamics
%   state := the current hybrid state where 
%       "-1" == unknown
%       "1 to 3" == the controllers/states defined above

% Setup initial parameters
persistent tRefStart % Stores when the system enters region 3)
persistent tanVecToIdentity % Store the direction to go towards the identity when entering region 3)
% Make sure to clear the above variable before calling this function for
% the first time using "clear <function_name>" in the workspace
gainRotationError=1;
gainVelocityError=1;
maxD = pi/2; % Max distance for contraction region from Identity
maxW = 2; % Max angular veloicty for contraction region
% gamma1 = 1.05; % Convergence rate of the reference position
gamma1 = 1.09;
RRef_new = eye(3);
wreference = zeros(3,1);

% Optional Parameters (from varargin)
ivarargin = 1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'gain'
            ivarargin=ivarargin+1;
            gainRotationError=varargin{ivarargin};
            gainVelocityError=varargin{ivarargin};
        case 'gainrotationerror'
            ivarargin=ivarargin+1;
            gainRotationError=varargin{ivarargin};
        case 'gainvelocityerror'
            ivarargin=ivarargin+1;
            gainVelocityError=varargin{ivarargin};
        case 'maxdistance'
            ivarargin=ivarargin+1;
            maxD = varargin{ivarargin};
        case 'maxangularspeed'
            ivarargin=ivarargin+1;
            maxW = varargin{ivarargin};
        otherwise
            % Skip data of unused entry
            ivarargin=ivarargin+1;
    end
    ivarargin=ivarargin+1;
end

% Determine which controller to use
[R,w]=rotDyn_stateUnpack(x);
V = rotBundle_Lyapunov(R,w,abs(gainRotationError),gainVelocityError);
if (norm(w) > maxW)
    % Inside 3) Vel controller region (goal is to slow down enough for PD
    % controller to take over)
    state = 3;
    dw = rotDyn_controlP_velocity(x,varargin{:});
elseif (V < min(gainRotationError*maxD^2, maxW^2))% && isempty(tRefStart))
    % Inside 1) PD controller region
    state = 1;
    dw = rotDyn_controlPD_contraction(x,varargin{:});
else
    % Inside 2) PD with reference controller region
    state = 2;
    % Update the reference trajectory (from R=Identiy, w = 0) to something
    % inside the contraction region centered at x
    if isempty(tRefStart)
        tRefStart = t; % Store the first time the system enters this region
        % Determine the direction to go towards the identity
        tanVecToIdentity = rot_vee(R,rot_log(R,eye(3)));
        if (dot(tanVecToIdentity, w) <= 0)
            % if the previous direction doesn't align with system's current
            % velocity then go along the other geodesic
            tanVecToIdentity = -(2*pi-norm(tanVecToIdentity))*tanVecToIdentity/norm(tanVecToIdentity);
        end
    end
    % RRef_new requires tangent vector from eye(3) to R, so take negative
    % of tanVecToIdentity
    RRef_new = rot_expVec(eye(3),-exp(-gamma1*(t-tRefStart))*tanVecToIdentity);
    % Loop through and modify RReference
    for i = 1:2:length(varargin)
        switch lower(varargin{i})
            case 'rreference'
                varargin{i+1} = RRef_new;
                break;
        end
    end
    % Add a reference velocity
    wreference = gamma1*exp(-gamma1*(t-tRefStart))*tanVecToIdentity;
    wOpts = {'wreference',wreference};
    varargin=[varargin, wOpts];
    dw = rotDyn_controlPD_contraction(x,varargin{:});
end

end