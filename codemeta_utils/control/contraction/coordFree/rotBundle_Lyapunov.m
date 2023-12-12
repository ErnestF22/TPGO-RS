function [V,Vdot] = rotBundle_Lyapunov(R,w,kd,kv,varargin)
% Compute the Lyapunov function for the controller
% INPUTS:
%   R := Current rotation (array of [3x3] matrices)
%   w := Current velocity (array of [3x1] vectors)
%   kd := Position error gain (scalar)
%   kv := Velocity error gain (scalar)
% OUTPUTS:
%   V := Lyapunov value(s) (should be >= 0)
%   Vdot := Derivative of V

% Paramters
RReference = eye(3);
wReference = zeros(3,1);

% Optional Parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'rreference'
            % Target Position
            ivarargin=ivarargin+1;
            RReference=varargin{ivarargin};
        case 'wreference'
            % Target angular velocity
            ivarargin=ivarargin+1;
            wReference=varargin{ivarargin};
        otherwise
            % Skip data of unused entry
            ivarargin=ivarargin+1;
    end
    ivarargin=ivarargin+1;
end

kd = abs(kd); %just incase
V = (kd*vecnorm(rot_vee(RReference,rot_log(RReference,R))).^2) + (vecnorm(wReference-w).^2);
Vdot = -kv*vecnorm(w).^2;

% if (exist('state','var'))
%     for i = 1:length(state)
%         [V(state(i)),Vdot(state(i))] = rotBundle_LyapunovW(w(:,state(i)),kv);
%     end
% end
end