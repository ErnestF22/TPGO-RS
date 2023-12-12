function [kd, kv, beta, m, gridData] = rotBundle_gainGridSearch(Dist,kd_list,kv_list, magW)
% Do a grid search of k_{d} and k_{v} gains to find best convergence rate
% INPUTS:
%   Dist := scalar argument for the R related matrix
%   kd_list := array of kd gains to try
%   kv_list := array of kv gains to try
%   magW := the maximum angular speed
% OUTPUTS:
%   beta := best feasible convergence rate
%   kd := best kd for beta 
%   kv := best kv for beta
%   m := the 3 values for the metric on TSO(3)

% Default parameters
beta = 0; % best convergematnce rate
kd = 0; % best kd for beta
kv = 0; % best kv for beta
m = zeros(2);
gridData = zeros(length(kd_list)*length(kv_list),3);

% Load InProgress
if exist('InProgressData.mat','file')
    load('InProgressData.mat');
end

% Grid Search
for kd_test = kd_list
    for kv_test = kv_list
        completedKD = find(gridData(:,1)==kd_test); 
        completedKV = find(gridData(:,2)==kv_test); 
        if ~isempty(intersect(completedKD,completedKV))
            fprintf('Skipping kd=%0.3f, kv=%0.3f \n',kd_test,kv_test);
            continue;
        end
        fprintf('Checking kd=%0.3f, kv=%0.3f \n',kd_test,kv_test);
        % run bisection search for best convergence rate
        [m_test, beta_test] = rotBundle_betaBisection(0,Dist,kv_test,kd_test,magW); 
        if (beta_test > beta)
            % store the best convergence rate
            kd = kd_test;
            kv = kv_test;
            beta = beta_test;
            m = m_test;
            fprintf('\tFound better results kd=%0.3f, kv=%0.3f, beta=%0.3f\n',...
                kd,kv,beta);
        end
        gridData = [gridData; [kd_test, kv_test, beta_test]];
    end
end
end

