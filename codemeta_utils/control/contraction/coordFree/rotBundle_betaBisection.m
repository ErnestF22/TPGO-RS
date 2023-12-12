function [ m, betaOut ] = rotBundle_betaBisection( minD, maxD, kv, kd, magW )
% Find the largest beta (convergence rate) and accompanying m such that 
% M <= 0 using gershgorin bounds
% INPUTS:
%   minD := minimum distance using log(R)
%   maxD := max distance using log(R)
%   kv := a positive scalar for the angular velocity error gain
%   kd := a negatice scalar for the rotation error gain
%   magW := max angular speed
% OUTPUTS:
%   m := the 3 values for the metric on TSO(3)
%   betaOut := the best convergence rate
%% Default Parameters
BETA_TOL = 1e-5;
MAX_ITER = 100;

%% Begin bisection search
% checkAboveFlag = false; %indicate which way to test the beta limit (true UB/current, false LB/current)
beta_UB = 10; %current beta is not feasible
beta_LB = 1e-7; %previous beta is feasible
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
        betaOut = testBeta;
    else
        %beta is not feasible, next check is between current and LB
%         checkAboveFlag = false;
        beta_UB = testBeta;
    end
    
    % Stop condition
    if (abs(beta_UB - beta_LB) < BETA_TOL)% && flag == true)
        fprintf('Bisection converged...\n');
        break;
    end
        
    % Find the average of the UB/LB to test on next iter
    testBeta = (beta_UB+beta_LB)/2; %test beta is the average of the UB/LB
end

if (i == MAX_ITER)
    fprintf('Reached max iter...\n');
end
