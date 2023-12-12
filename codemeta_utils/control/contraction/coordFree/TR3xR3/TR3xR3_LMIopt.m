function [ m, flagFeasible ] = TR3xR3_LMIopt( e_A, e_B, kv, kd, kp, beta )
% Use the LMI as the objective and constraint. NOTE this method has not
% been shown to have max and min contraction eigenvalues at the max and min
% eigenvalues of the Hessisan.
% Bound method: Split into two matrix, one depends on rho(x,xd), other
% depends on rho(xd,0)
% INPUT:
%   e_A := eigenvalue range for the A matrix corresponding to rho(x,xd)
%      [min;max]
%   e_B := eigenvalue range for the B matrix correspondding to rho(xd,0)
%       [min;max]
%   kd,kv,kp := positive scalar position/velocity/ref. position gain
%   beta := positive scalar for min convergence rate
% OUTPUT:
%   m := [m1 m2 m6;m2 m3 m5;m6 m5 m4] values
%   flagFeasible := boolean where true means the system is feasible

%% Default Parameters
flagFeasible = false; %set to false by default
%% Optimize using LMI as constraints
cvx_begin sdp quiet
    variable m(3,3) semidefinite
    % Re-assign vars
    m1 = m(1,1); m2 = m(1,2); m3 = m(2,2); m4 = m(3,3); m5 = m(2,3); m6 = m(1,3);
    % Define the discs of A
    DA_1 = @(lambda_A) -2*m2*kd*lambda_A + m1*beta ...
        + abs(-m3*kd*lambda_A + m1-m2*kv + beta*m2)...
        + abs((m2-m5)*kd*lambda_A + beta*m6);
    DA_2 = @(lambda_A) 2*(m2-m3*kv) + m3*beta ...
        + abs(-m3*kd*lambda_A + m1-m2*kv + beta*m2)...
        + abs(m3*kd*lambda_A + m6-m5*kv + beta*m5);
    DA_3 = @(lambda_A) 2*m5*kd*lambda_A + m4*beta ...
        + abs((m2-m5)*kd*lambda_A + beta*m6)...
        + abs(m3*kd*lambda_A + m6-m5*kv + beta*m5);
% %     DA_1 = @(lambda_A) -2*m2*kd*lambda_A ...
% %         + abs(-m3*kd*lambda_A+m1-m2*kv)...
% %         + abs((m2-m5)*kd*lambda_A);
% %     DA_2 = @(lambda_A) 0 ...
% %         + abs(-m3*kd*lambda_A)...
% %         + abs(m3*kd*lambda_A);
% %     DA_3 = @(lambda_A) 2*m5*kd*lambda_A ...
% %         + abs((m2-m5)*kd*lambda_A)...
% %         + abs(m3*kd*lambda_A);
    % Define the discs of B
    DB_1 = @(lambda_B) 0 + abs(-m6*kp*lambda_B);
    DB_2 = @(lambda_B) 0 + abs(-m5*kp*lambda_B);
    DB_3 = @(lambda_B) -2*m4*kp*lambda_B ...
        + abs(-m5*kp*lambda_B)...
        + abs(-m6*kp*lambda_B);
% %     DB_1 = @(lambda_B) m1*beta + abs(-m6*kp*lambda_B + m6*beta);
% %     DB_2 = @(lambda_B) m3*beta + 2*(m2-m3*kv) + abs(m1-m2*kv) + abs(-m5*kp*lambda_B + m5*beta + m6-m5*kv);
% %     DB_3 = @(lambda_B) -2*m4*kp*lambda_B + m4*beta...
% %         + abs(-m6*kp*lambda_B + m6*beta)...
% %         + abs(-m5*kp*lambda_B + m5*beta);
    % Define the two set we need max of
    set_A = [DA_1(e_A(1));DA_1(e_A(2));...
        DA_2(e_A(1));DA_2(e_A(2));...
        DA_3(e_A(1));DA_3(e_A(2))];
    set_B = [DB_1(e_B(1));DB_1(e_B(2));...
        DB_2(e_B(1));DB_2(e_B(2));...
        DB_3(e_B(1));DB_3(e_B(2))];

    % Setup feasibility SDP (no objective)
    m(1,1) == 1 % fix one entry of metric tensor for scale
    max(set_A) + max(set_B) <= 0
cvx_end

if  strcmp(lower(cvx_status),'solved')
    eM = eig(m);
    if (all(eM > 0) && max(set_A)+max(set_B) <= 0)
        flagFeasible = true;
    end
end

end

