function [dw] = rotRefOpt_control(x,varargin)
% LAST EDITED: Oct. 5th 2020
% Function to solve optimal controller gains on TSO(3)xSO(3) then feed new
% gains into rotDyn_controlPD.m
% Steps:
% 1) Solve optimal change in gains where the direction vector is the
% current system dynamics
% 2) Update the gains according to the optimal change
% 3) Compute control with rotDyn_controlPD  
% INPUTS:
%   x := current state [R;\omega;RRef;kd;kv;kref] (assumes this is the
%       'augmentedsystem')
% OUTPUTS:
%   dw := the updated controller for \omega [3x1] vector

%% Setup
%optional parameters
ivarargin=1;
len_varargin = length(varargin);
while ivarargin<=len_varargin
    switch lower(varargin{ivarargin})    
        case 'gainrotationerror'
            ivarargin=ivarargin+1;
%             gainRotationError=varargin{ivarargin}; % Original gain
            gainRotationError_idx = ivarargin; % Store for update later
        case 'gainvelocityerror'
            ivarargin=ivarargin+1;
%             gainVelocityError=varargin{ivarargin}; % Original gain
            gainVelocityError_idx = ivarargin; % Store for update later
        case 'gainschedulecontroller'
            ivarargin=ivarargin+1;
            GainScheduledParams = varargin{ivarargin};
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 <= len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end
    end
    ivarargin=ivarargin+1;
end

%% Controller (including additional parameters)
% Extract the updated gains
if exist('GainScheduledParams','var')
    % Update gains base on Gain Schedule, assumes using local controller
    % once the system is no longer pi distance away
   x = rotRefOpt_GainSchedule_Checker(x,GainScheduledParams);
end

% Update sent data using new gains for controller
[~,~,~,kd,kv,~,~,~] = rotDynOptAll_stateUnpack(x);
varargin{gainRotationError_idx} = kd;
varargin{gainVelocityError_idx} = kv;

% Compute the control with the updated gains
dw = rotDyn_controlPD(x,varargin{:});

end

