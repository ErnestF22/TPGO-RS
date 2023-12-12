function [ betaOut ] = rotBundle_betaBisection_new(maxD, kv, kd, magW, varargin )
% Find the largest beta (convergence rate) and accompanying m such that 
% M <= 0 using gershgorin bounds
% INPUTS:
%   minD := minimum distance using log(R)
%   maxD := max distance using log(R)
%   kv := a positive scalar for the angular velocity error gain
%   kd := a negatice scalar for the rotation error gain
%   magW := max angular speed
% OUTPUTS:
%   betaOut := the best convergence rate
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
if kd < 0 || kv < 0 || kd > MAXGAIN || kv > MAXGAIN
    betaOut = -1e6;
    return;
end
%% Begin bisection search
betaOut = beta_LB; % This is the last known good convergence rate
testBeta = beta_UB; %test beta is the average of the UB/LB

for i = 1:MAX_ITER
    % Test beta using the 3 possible bounds (just need 1 solution)
    [m, flag] = rotBundle_contractionOpt(maxD, kv, kd, testBeta, magW);
    
    % Update the best guess
    if (flag == true)
        %beta is feasible, next check is between current and UB
%         checkAboveFlag = true;
        beta_LB = testBeta;
        % Update if current best beta is less than new result
        if betaOut < testBeta
            betaOut = testBeta;
        end
    else
        %beta is not feasible, next check is between current and LB
%         checkAboveFlag = false;
        beta_UB = testBeta;
    end
    
    % Stop condition
    if (abs(beta_UB - beta_LB) < BETA_TOL)% && flag == true)
        fprintf('Bisection converged..., best beta = %0.3f\n',betaOut);
        break;
    end
        
    % Find the average of the UB/LB to test on next iter
    testBeta = (beta_UB+beta_LB)/2; %test beta is the average of the UB/LB
end
fprintf('Bisection complete, results (kd,kv,beta): %f, %f, %f\n',kd,kv,betaOut);
