function [x] = rotRefOpt_GainSchedule_Checker(x, GainScheduledParams)
% Logic to determine if the gains and reference trajectory should change
% because the state is in a "better" convergence region.
% INPUTS:
%   x := the [34x1] vector with states [R,w,RRef,kd,kv,kref,M_nn,beta]
%   GainScheduledParams := a struct with fields
%       [kd,kv,kref,beta,mag_R,mag_RRef,mag_W,M_nn]. Each line is defines a
%       new contraction region.
% OUTPUTS:
%   x := the updated [34x1] vector with states [R,w,RRef,kd,kv,kref,M_nn,beta]

% Get required state information
[R,w,~,~,~,~,~,~] = rotDynOptAll_stateUnpack(x);

% Find the better convergence region and set new gains
% NOTE: when used inside ode45, the state change here does not propagate to
% the acutal system. So it is possible that the states "jump" between the
% static controller and the gain schedule, however it should not. If this
% is a problem, the model/control functions should be rewritten
distError_R = norm(rot_log(R,eye(3)));
for i = 1:length(GainScheduledParams)
    if distError_R <= GainScheduledParams(i).mag_R
        % Extract new gains
        kd_new = GainScheduledParams(i).kd;
        kv_new = GainScheduledParams(i).kv;
        beta_new = GainScheduledParams(i).beta;
        M_nn = GainScheduledParams(i).M_nn;
        M_nn = [M_nn, [0;0];[0,0 1]]; % Fix to fit rotRef form
        % Test if current state is contracting with new gains
        % Get the contraction matrix on TSO(3) as function of R
        M_contraction = rotBundle_contractionMat(@(R) R*hat3(w),beta_new,kd_new,kv_new,[M_nn(1,1),M_nn(1,2),M_nn(2,2)],'sym');
        M_contraction_R = M_contraction(R);
        if max(eig(M_contraction_R)<=1e-5)
            kd = kd_new;
            kv = kv_new;
            RRef = eye(3);
            b = beta_new;
            x = rotDynOptAll_statePack(R,w,RRef,kd,kv,0,M_nn,b);
        end
    end
end
% NOTE: If no convergence region is met, the x variable is unchanged
end

