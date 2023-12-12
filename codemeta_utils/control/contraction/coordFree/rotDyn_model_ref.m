function [dx] = rotDyn_model_ref(x,varargin)
% Second order dynamical model for a rigid body rotation following a
% reference trajectory. Dynamics evolve on TSO(3)xSO(3)
% INPUTS:
%   x := [21x1] current state x (packed). The first 12 elements are on
%       TSO(3), last 9 are on SO(3)
% OUTPUTS:
%   dx := derivative of the state, \dot{x} (packed)

gainReference = 1; 

% Load optional parameters
ivarargin=1;
len_varargin = length(varargin);
while ivarargin<=len_varargin
    switch lower(varargin{ivarargin})
        case 'gainreference'
            % Initial gains must be sent as a 3x1 vector
            ivarargin=ivarargin+1;
            gainReference = varargin{ivarargin};
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 <= len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end
    end
    ivarargin=ivarargin+1;
end

% Call rotDyn_model to update dynamic on TSO(3)
dx_TSO3 = rotDyn_model(x(1:12),varargin{:});

% Calculate derivative for the reference trajectory on SO(3)
[RReference,~] = rotDyn_stateUnpack([x(13:21);zeros(3,1)]);
dx_SO3 = rotDyn_statePack(gainReference*rot_log(RReference,eye(3)),zeros(3,1));

% Combine all states
dx = [dx_TSO3;dx_SO3(1:9)];
end

