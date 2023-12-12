function [betaOut] = rotRef_betaBisection(mag_R,mag_RRef,kd,kv,kp,mag_W,varargin)
%Find the largest beta and accompanying m such that the contraction matrix
% on TSO3xSO3 <= 0 using gershgorin disc as bounds
% INPUTS:
%   L_R := vector of [min;max] eigenvalues of Dlog_RRef_R
%   L_RRef := vector of [min;max] eigenvalues of Dlog_RRef_I
%   mag_R := positive scalar representing max distance from R to RRef
%   mag_RRef := positive scalar representing max distance from RRef to I
%   kd, kv, kp := scalar positive gains
%   magW := positive scalar representing maximum velocity magnitude
% OUTPUTS:
%   betaOut := scalar for best convergence rate with a solution, positive
%       means convergence and negative means divergence

%% Default Parameters
BETA_TOL = 1e-5;
MAX_ITER = 20;
beta_UB = 10; %current beta is not feasible
beta_LB = 1e-7; %previous beta is feasible
MAXGAIN = inf; % Upper bound for the gains

% Load optional parameters
ivarargin=1;
len_varargin = length(varargin);
while ivarargin<=len_varargin
    switch lower(varargin{ivarargin})
        case 'beta_lb'
            % Set a new lower beta to test
            ivarargin=ivarargin+1;
            beta_LB = varargin{ivarargin};
        case 'beta_ub'
            % Set a new lower beta to test
            ivarargin=ivarargin+1;
            beta_UB = varargin{ivarargin};            
        case 'maxgain'
            % Set an upper bound on the gains
            ivarargin=ivarargin+1;
            MAXGAIN = varargin{ivarargin};
        case 'bisection_maxiter'
            % Set the max bisection iteration
            ivarargin=ivarargin+1;
            MAX_ITER = varargin{ivarargin};
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 <= len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end    
    end
    ivarargin=ivarargin+1;
end

% If given negative gains or gains greater than MAXGAIN, the answer is not feasible
if kd < 0 || kv < 0 || kp < 0 || kd > MAXGAIN || kv > MAXGAIN || kp > MAXGAIN
    betaOut = -1e6;
    return;
end

%% Begin bisection search
% checkAboveFlag = false; %indicate which way to test the beta limit (true UB/current, false LB/current)
betaOut = beta_LB; % Last known good beta
betaTest = beta_UB; % Beta to test

for i = 1:MAX_ITER
    [m, flagFeas] = rotRef_contractionOpt_new(mag_R,mag_RRef,kd,kv,kp,mag_W,betaTest,varargin{:});
    
    if (flagFeas == true)
        %beta is feasible, next check is between current and UB
        beta_LB = betaTest;
        % Update betaOut if better results
        if betaOut < betaTest
            betaOut = betaTest;
        end
    else
        %beta is not feasible, next check is between current and LB
        beta_UB = betaTest;
    end
    
    %stop condition
    if (abs(beta_UB - beta_LB) < BETA_TOL) %&& flagFeas == true)
        fprintf('Bisection converged..., best beta = %0.3f\n',betaOut);
        break;
    end
        
    %find the average of the UB/LB to test on next iter
    betaTest = (beta_UB+beta_LB)/2; %test beta is the average of the UB/LB
end
fprintf('Bisection complete, results (kd,kv,kref,beta): %f, %f, %f, %f\n',kd,kv,kp,betaOut);

