function [R,w,RRef,t,varargout] = rotRef_simulation(mag_R,mag_RRef,mag_W,varargin)
% Simulate the dynamics on TSO(3)xSO(3)
% INPUTS:
%   mag_R := positive scalar representing max distance from R to RRef
%   mag_RRef := positive scalar representing max distance from RRef to I
%   magW := positive scalar representing maximum velocity magnitude
% OUTPUTS:
%   R := [3x3] Matrix array of the current state 
%   w := [3x1] Vector array of the current angular velocity
%   RRef := [3x3] Matrix array of the reference trajectory
%   t := simulation time array

% Define default parameters
% resetRands();
kd = 106.6667; kv = 74.6667; kref = 0.9833;
TFinal = 10;
J=diag([5;2;1]);
w0=0.02*cnormalize(randn(3,1));
% Define reference rot as a random point at mag_RRef distance away
RRef0 = rot_exp(eye(3),hat3(cnormalize(randn(3,1))*mag_RRef));
% Define current state as mag_R distance away from R_ref0
R0 = rot_exp(RRef0,RRef0*hat3(cnormalize(randn(3,1))*mag_R));

% Load optional parameters
ivarargin=1;
len_varargin = length(varargin);
while ivarargin<=len_varargin
    switch lower(varargin{ivarargin})
        case 'init_vel'
            % Set the initial angular velocity
            ivarargin=ivarargin+1;
            w0 = varargin{ivarargin};
        case 'init_r'
            % Set the initial rotation R
            ivarargin=ivarargin+1;
            R0 = varargin{ivarargin};
        case 'init_rref'
            % Set the initial reference rotation R_ref
            ivarargin=ivarargin+1;
            RRef0 = varargin{ivarargin};
        case 'endtime'
            % Set the simulation end time
            ivarargin=ivarargin+1;
            TFinal = varargin{ivarargin};
        case 'gains'
            % Initial gains must be sent as a 3x1 vector
            ivarargin=ivarargin+1;
            gains0 = varargin{ivarargin};
            kd = gains0(1); kv = gains0(2); kref = gains0(3);
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 <= len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end
    end
    ivarargin=ivarargin+1;
end

% Define simulation parameters
optsModel={'inertiaMatrix',J,'gainreference',kref};
optsCost=[optsModel {'gainRotationError',kd}];
optsControl=[optsCost {'flagCancelDynamics',true,'flagVelocityInertiaWeight',false,'gainVelocityError',kv,'augmentedsystem'}];
x0=rotDyn_statePack(R0,w0,'augmentedsystem',RRef0);
% Simulate the system using the given gains
control=@(t,x) rotDyn_controlPD(x,optsControl{:});
closedLoop=@(t,x) rotDyn_model_ref(x,optsModel{:},'torque',control(t,x));

% Simulate
figure
optsOde=odeset('OutputFcn',@odeplot,'MaxStep',0.01);
% optsOde=odeset('MaxStep',0.01);
switch 1
    case 1
        [t,x]=ode45(closedLoop,[0 TFinal],x0,optsOde);
    case 2
        [t,x]=rotDyn_odeEuler(closedLoop,[0 TFinal],x0,optsOde);
end
x=x';
[R,w,RRef]=rotDyn_stateUnpack(x,'augmentedsystem');

% Add additional things to output
varargout{1} = x;
varargout{2} = control;
end

