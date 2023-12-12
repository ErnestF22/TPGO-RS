function [dx] = rotDynOptAll_model(x, varargin)
% LAST EDITED: Nov. 13, 2020 by Bee Vang
% Second order dynamical model for rigid body rotation. Dynamics evolve on
% TSO(3)xSO(3) with dynmical gains, metric, and convergence rate
% INPUTS:
%   x := [34x1] vector with states [R,w,RRef,kd,kv,kref,M_nn,beta]
% OUTPUTS:
%   dx := derivate of teh state, \dot{x} (packed)

%% Load optional parameters (non atm)
flagSkipOptParameters = false;
ivarargin=1;
len_varargin = length(varargin);
while ivarargin<=len_varargin
    switch lower(varargin{ivarargin})
        case 'skip_optimization_parameters'
            % Skip optimization process, used for simulating unchanged
            % parameter system
            flagSkipOptParameters = true;
        case 'gainschedulecontroller'
            flagSkipOptParameters = true;
            ivarargin=ivarargin+1; % Increment 1 b/c this parameter also includes another parameter
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 <= len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end
    end
    ivarargin=ivarargin+1;
end

%% Extracct current state
[R,w,RRef,kd,kv,kref,M_nonnatural,~] = rotDynOptAll_stateUnpack(x);

% Check that the metric tensor is still positive definite
if any(eig(M_nonnatural))<= 0
    error('Metric tensor no longer positive definite');
end

%% Update system dynamics
% Call rotDyn_model to update dynamic on TSO(3)
% dx_TSO3 = rotDyn_statePack(R*hat3(w), Gamma);
dx_TSO3 = rotDyn_model(x(1:12),varargin{:});

% Calculate derivative for the reference trajectory on SO(3)
dx_SO3 = rotDyn_statePack(kref*rot_log(RRef,eye(3)),zeros(3,1));

%% Determine optimal gain, metric, and convergence rates
if flagSkipOptParameters
    dparams = zeros(13,1); % All parameters unchanged
else
    dparams = rotRefOptAll_optProblem_bounds_boundNeg(x, varargin{:});
end
%% Combine all rates
dx = [dx_TSO3;dx_SO3(1:9);dparams];
end

