function [ betaOut ] = TR3xR3_betaBisection( e_A, e_B, kv, kd, kp, varargin )
%Find the largest beta and accompanying m such that G <= 0 using gershgorin
%as bounds

%% Default Parameters
BETA_TOL = 1e-5;
MAX_ITER = 100;
beta_UB = 10; %current beta is not feasible
beta_LB = 1e-7; %previous beta is feasible
maxGain = inf; % Upper bound for the gains
flagNonAugSys = false; % Flag to solve for the non-augmented system

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
            maxGain = varargin{ivarargin};
        case 'nonaugmentedsystem'
            % Solve the LMI for the standard non-augmented system
            % The path to the LMI function must be added, 
            % IE run the commend: addpath('../../pointmass_3D/');
            flagNonAugSys = true;
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 <= len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end    
    end
    ivarargin=ivarargin+1;
end

% If given negative gains or gains greater than maxGain, the answer is not feasible
if kd < 0 || kv < 0 || kp < 0 || kd > maxGain || kv > maxGain || kp > maxGain
    betaOut = -1e6;
    return;
end

%% Begin bisection search
% checkAboveFlag = false; %indicate which way to test the beta limit (true UB/current, false LB/current)
betaOut = 0; % Last known good beta
betaTest = beta_UB; % Beta to test

for i = 1:MAX_ITER
    %test beta
    if flagNonAugSys
        [m, flagFeas] = contraction3D_LMIopt(e_A(1), e_A(2), kv, kd, betaTest);
    else
        [m, flagFeas] = TR3xR3_LMIopt2(e_A, e_B, kv, kd, kp, betaTest);
    end

    if (flagFeas == true)
        %beta is feasible, next check is between current and UB
%         checkAboveFlag = true;
        beta_LB = betaTest;
        betaOut = betaTest;
    else
        %beta is not feasible, next check is between current and LB
%         checkAboveFlag = false;
        beta_UB = betaTest;
    end
    
    %stop condition
    if (abs(beta_UB - beta_LB) < BETA_TOL && flagFeas == true)
        fprintf('Stop condition met...\n');
        break;
    end
        
    %find the average of the UB/LB to test on next iter
    betaTest = (beta_UB+beta_LB)/2; %test beta is the average of the UB/LB
end

