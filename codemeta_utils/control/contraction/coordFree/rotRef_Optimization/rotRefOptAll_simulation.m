function [x,t,varargout] = rotRefOptAll_simulation(mag_R,mag_RRef,mag_W,M_nonnatural,varargin)
% LAST EDITED: Nov. 13, 2020 by Bee Vang
% Simulate the dynamics on TSO(3)xSO(3) where the gains, metric, and
% convergence rates are dynamic
% INPUTS:
%   mag_R := positive scalar representing max distance from R to RRef
%   mag_RRef := positive scalar representing max distance from RRef to I
%   magW := positive scalar representing maximum velocity magnitude
%   M_nonnatural := a [3x3] nonnatural metric on TSO(3)xSO(3)
% OUTPUTS:
%   x := [34x1] vector with states [R,w,RRef,kd,kv,kref,M_nn,beta]
%   t := simulation time array

%% Generate default parameters  (using result from CDC 2020)
% resetRands();
kd = 106.6667; kv = 74.6667; kref = 0.9833;
TFinal = 10;
J=diag([5;2;1]);
w0=mag_W*cnormalize(randn(3,1));
% Define reference rot as a random point at mag_RRef distance away
RRef0 = rot_exp(eye(3),hat3(cnormalize(randn(3,1))*mag_RRef));
% Define current state as mag_R distance away from R_ref0
R0 = rot_exp(RRef0,RRef0*hat3(cnormalize(randn(3,1))*mag_R));
% Set an initial beta
b0 = 0.4022;

%% Load optional parameters
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
        case 'convergence_rate'
            % Set an initial convergence rate, result should be feasible
            ivarargin=ivarargin+1;
            b0 = varargin{ivarargin};
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 <= len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end
    end
    ivarargin=ivarargin+1;
end

%% Define simulation parameters
optsModel={'inertiaMatrix',J,'gainreference',kref};
optsCost=[optsModel {'gainRotationError',kd}];
optsControl=[optsCost {'flagCancelDynamics',true,'flagVelocityInertiaWeight',false,'gainVelocityError',kv,'augmentedsystem'}];
x0=rotDynOptAll_statePack(R0,w0,RRef0,kd,kv,kref,M_nonnatural,b0);
% Simulate the system using the given gains
control=@(t,x) rotRefOpt_control(x,optsControl{:},varargin{:});
closedLoop=@(t,x) rotDynOptAll_model(x,optsModel{:},'torque',control(t,x),varargin{:});

%% Simulate
figure
optsOde=odeset('OutputFcn',@odeplot,'MaxStep',0.01);
% optsOde=odeset('MaxStep',0.01);
switch 1
    case 1
        [t,x]=ode45(closedLoop,[0 TFinal],x0,optsOde);
    case 2
        [t,x]=rotDyn_odeEuler(closedLoop,[0 TFinal],x0,optsOde);
end
% Transpose for processing and usage in other functions
x=x';

% Add additional things to output
varargout{1} = control;
end

