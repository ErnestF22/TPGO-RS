function [ m, betaOut ] = contraction3D_betaBisection( eMin, eMax, kvIn, kdIn, varargin )
%Find the largest beta and accompanying m such that G <= 0 using gershgorin
%as bounds

%% Default Parameters
BETA_TOL = 1e-5;
MAX_ITER = 100;
flagLMI = false; % Flaag to us LMI method or gershgorin method
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'lmi'
            % Use LMI constraint
            flagLMI = true;
        otherwise    
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

%% Begin bisection search
% checkAboveFlag = false; %indicate which way to test the beta limit (true UB/current, false LB/current)
beta_UB = 10; %current beta is not feasible
beta_LB = 1e-7; %previous beta is feasible
betaOut = beta_LB; % Last known good beta
betaTest = beta_UB; % Beta to test

for i = 1:MAX_ITER
    %test beta
    if ~flagLMI
        [m, flagFeas] = contraction3D_gershgorin(eMin, eMax, kvIn, kdIn, betaTest);
    else
        [m, flagFeas] = contraction3D_LMIopt(eMin, eMax, kvIn, kdIn, betaTest);
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
    
%     %adjust UB/LB
%     if (checkAboveFlag == true)
%         beta_LB = beta_test;
%     else
%         beta_UB = beta_test;
%     end

    %stop condition
    if (abs(beta_UB - beta_LB) < BETA_TOL && flagFeas == true)
        fprintf('Stop condition met...\n');
        break;
    end
        
    %find the average of the UB/LB to test on next iter
    betaTest = (beta_UB+beta_LB)/2; %test beta is the average of the UB/LB
end

