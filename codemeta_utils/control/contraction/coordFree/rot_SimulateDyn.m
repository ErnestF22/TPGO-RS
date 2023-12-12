function [R,w,t,u,state, RRef_new, wreference] = rot_SimulateDyn(R0, w0, TFinal, kd, kv, J, beta, m)
% Simulate the attitude system using hybrid controllers
% INPUTS:
%   R0 := Initial position on SO(3) as [3x3] matrix
%   w0 := Initial velocity as a [3x1] vector
%   TFinal := Max simulation time (scalar)
%   kd := Position error gain (scalar)
%   kv := Velocity error gain (scalar)
%   J := Inertial matrix [3x3]
% OUPUTS:
%   R := Trajectory on SO(3) as array for [3x3] matrices
%   w := Trajectory on T_{R}SO(3) as [3x1] vectors
%   t := Simulation times
%   state := indicates which controller is active

% Set the simulation parameters
resetRands();
gainVel = kv; % velocity error gain
RReference = eye(3);
optsModel={'inertiaMatrix',J};
optsReference={'RReference',RReference};
optsCost=[optsModel optsReference {'gainRotationError',-kd}];
optsControl=[optsCost {'flagCancelDynamics',true,'flagVelocityInertiaWeight',false,'gainVelocityError',gainVel}];
x0=rotDyn_statePack(R0(:),w0);

% Define simulation functions
clear rotDyn_hybrid; % clear any persistent info
control=@(t,x) rotDyn_hybrid(x,t,optsControl{:});
closedLoop=@(t,x) rotDyn_model(x,optsModel{:},'torque',control(t,x));

% Simulate
% figure(1)
optsOde=odeset('OutputFcn',@odeplot,'MaxStep',0.01);
switch 1
    case 1
        [t,x]=ode45(closedLoop,[0 TFinal],x0,optsOde);
    case 2
        [t,x]=rotDyn_odeEuler(closedLoop,[0 TFinal],x0,optsOde);
end

x=x';
[R,w]=rotDyn_stateUnpack(x);

% recompute u
for i = 1:length(t)
    [u(:,i), state(:,i), RRef_new(:,:,i), wreference(:,i)] = control(t(i),x(:,i));
end
end