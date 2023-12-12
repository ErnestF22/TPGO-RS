function [kd,kv,kp,beta,m,gridData] = par_TR3xR3_gainGridSearch_outterLoop(e_A, e_B, kd_list, kv_list, kp_list)
% Perform a gridsearch of kd,kv,kp gains to find best convergence rate.
% INPUTS:
%   e_A := eigenvalue range for the A matrix corresponding to rho(x,xd)
%      [min;max]
%   e_B := eigenvalue range for the B matrix correspondding to rho(xd,0)
%       [min;max]
%   kd_list, kv_list, kp_list := array of gains to test
% OUTPUTS:
%   kd,kv,kp := combo of gains with best result
%   beta := best convergence rate
%   gridData := the raw data as a structure

% Default params

% Create 2D data where the columns are [kd;kv;kp;beta;m1;m2;m3;m4;m5;m6]
% Store gain list as a cell array
allGains = {kd_list, kv_list, kp_list};
% Then create a 2D matrix of all posssible gains
gridData = cell(1,numel(allGains));
[gridData{:}] = ndgrid(allGains{:});
gridData = cellfun(@(X) reshape(X,[],1),gridData,'UniformOutput',false);
gridData = horzcat(gridData{:});
% extend gridData
gridData(:,4) = 0;
% Grid Search
parfor ii = 1:size(gridData,1)
    tempData=gridData(ii,:);
    kd_test = tempData(1); kv_test = tempData(2); kp_test = tempData(3);
    fprintf([datestr(datetime,'yyyymmdd HH:MM:ss ') 'Starting loop for kd=%0.3f, kv=%0.3f, kp=%0.3f\n'],...
        kd_test, kv_test, kp_test);
    % Run the bisection search
    beta_test = TR3xR3_betaBisection( e_A, e_B, kv_test, kd_test, kp_test);
    % Store data
    tempData(4) = beta_test;
% %     tempData(5:end) = [m_test(1,1),m_test(1,2),m_test(2,2),m_test(3,3),m_test(2,3),m_test(1,3)];
    gridData(ii,:) = tempData;
    fprintf([datestr(datetime,'yyyymmdd HH:MM:ss ') 'Ending loop for kd=%0.3f, kv=%0.3f, kp=%0.3f\n'],...
        kd_test, kv_test, kp_test);
end

% Find best result and return it
[beta, idxMax_row] = max(gridData(:,4)); % The idx is the row with the best result
kd = gridData(idxMax_row,1);
kv = gridData(idxMax_row,2);
kp = gridData(idxMax_row,3);

    