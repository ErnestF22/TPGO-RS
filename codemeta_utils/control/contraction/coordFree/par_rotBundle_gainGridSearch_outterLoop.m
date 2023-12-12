function [kd, kv, beta, m, gridData] = par_rotBundle_gainGridSearch_outterLoop(Dist,kd_list,kv_list, magW)
% (PARALLEL version) Do a grid search of k_{d} and k_{v} gains to find best convergence rate
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
betaData = zeros(length(kd_list),length(kv_list)); %Rows are kd, columns are kv
mData = zeros(length(kd_list),length(kv_list),4); % Store the m_test matrix as a vector
gainData = zeros(length(kd_list),length(kv_list),2); % Store the [kd;kv] in vector 

% Grid Search
parfor ii = 1:length(kd_list)
    kd_test = kd_list(ii);
    temp_betaData = zeros(1,length(kv_list));
    temp_mData = zeros(4,length(kv_list));
    temp_gainData = zeros(2,length(kv_list));
    for jj = 1:length(kv_list)
        kv_test = kv_list(jj);
%         fprintf('Checking kd=%0.3f, kv=%0.3f \n',kd_test,kv_test);
        % run bisection search for best convergence rate
        [m_test, beta_test] = rotBundle_betaBisection(0,Dist,kv_test,kd_test,magW); 
        temp_betaData(jj) = beta_test;
        temp_mData(:,jj) = m_test(:);
        temp_gainData(:,jj) = [kd_test;kv_test];
    end
    
    betaData(ii,:) = temp_betaData;
    mData(ii,:,:) = temp_mData';
    gainData(ii,:,:) = temp_gainData';
end

% Find best result and return it
beta = max(betaData(:));
[idxMax, idyMax] = find(betaData==beta,1);
bestGains = squeeze(gainData(idxMax,idyMax,:));
kd = bestGains(1); kv = bestGains(2);
best_mMat = squeeze(mData(idxMax,idyMax,:));
m = reshape(best_mMat,2,2);

% Save data for later usage
gridData.betaData = betaData;
gridData.mData = mData;
gridData.gainData = gainData;
end