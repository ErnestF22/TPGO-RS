function [A_all_der_func] = rotRefOpt_allBounds_Matrix_der_simple_cvx(R,w,RRef,gains,M_nn,b,gains_dot,M_nn_dot,b_dot)
% Last Edited: Jan 8, 2021 by Bee Vang
% Compute the relaxed gershgorin disc bounds' derivative wrt the gains,
% metric, and convergence rate. Intended to be used within CVX_BEGIN,
% however it can be used to get a numerical result. Not intended as a
% general function handle generator.
% NOTE: To do this we relaxed abs(a+bi) <= abs(a) + abs(b)
% INPUTS:
%   R := Current rotation matrix
%   w := Current angular speed as [3x1] vector
%   RRef := Current reference rotation matrix
%   gains := The gains given as [kd;kv;kref]
%   M_nn := Metric tensor as [3x3] matrix
%   b := Convergence rate as scalar
%   gains_dot := The gains derivative as [kd_dot;kv_dot;kref_dot] (THIS CAN
%       BE A CVX VARIABLE)
%   M_nn_dot := Metric tensor derivative as [3x3] matrix with m4 ==
%       constant (THIS CAN BE A CVX VARIABLE)
%   beta_dot := Convergence rate derivative (THIS CAN BE A CVX VARIABLE)
% OUTPUTS:
%   A_all_der_func(gains_dot,Mnn_dot,beta_dot) := (a cvx handle of) all possible 
%       derivative of bounds on the contraction matrix eigenvalue on TSO(3)xSO(3)
%       wrt the gains, metric tensor, and contraction rate

%% Default Parameters
mag_R = norm(rot_log(R,RRef)); mag_RRef = norm(rot_log(RRef,eye(3)));
mag_W = norm(w);
kd = gains(1); kv = gains(2); kref = gains(3);
m1 = M_nn(1,1); m2 = M_nn(1,2); m3 = M_nn(2,2); m4 = M_nn(3,3);
m5 = M_nn(2,3); m6 = M_nn(1,3);
kd_dot = gains_dot(1); kv_dot = gains_dot(2); kref_dot = gains_dot(3);
m1_dot = M_nn_dot(1,1); m2_dot = M_nn_dot(1,2); m3_dot = M_nn_dot(2,2); m4_dot = M_nn_dot(3,3);
m5_dot = M_nn_dot(2,3); m6_dot = M_nn_dot(1,3);

%% Init matrices
% NOTE: we are computing [1 m1 m2 m3 m5 m6]*A*[[1 m1 m2 m3 m5 m6], but only
% the submatrix A(2:end,2:end) needs to be PSD which reduces to just the
% [2x2] matrix given as [m5 m6]*Bi*[m5;m6]
% Init matrix based on symbolic or not
% Populate the centroid terms
A_11 = -m2_dot*kd*mag_R/2*cot(mag_R/2) - m2*kd_dot*mag_R/2*cot(mag_R/2);
A1_mat = {A_11};
A_22 = m2_dot + m3_dot*(b-kv) + m3*(b_dot-kv_dot);
A2_mat = {A_22};
A_33 = m4*b_dot - kref_dot*m4*mag_RRef/2*cot(mag_RRef/2);
A3_mat = {A_33};

% Define possible terms resulting from absolute values
% For disc centered at D(1,1)
% Fill entries
tempMat1 = m1_dot*b + m1*b_dot;
tempMat2 = m1_dot*sqrt(3)*b + m1*sqrt(3)*b_dot...
    + m2_dot*sqrt(3)*(-1/4*mag_W^2)...
    + m5_dot*m6*sqrt(3)*(mag_W^2/(4*m4)) + m5*m6_dot*sqrt(3)*(mag_W^2/(4*m4));
D1_abs_fun{1} = {tempMat1,tempMat2,-tempMat2};
% D1_abs_fun{1} = {sqrt(3)*b*m1,...
%     -sqrt(3)*b*m1,...
%     sqrt(3)*(-1/4*(m2-m5*m6/m4)*mag_W^2+b*m1),...
%     -sqrt(3)*(-1/4*(m2-m5*m6/m4)*mag_W^2+b*m1)};
% 
% Fill entries
tempMat1 = m1_dot*sqrt(3)*1/2 + m2_dot*sqrt(3)*(-1/2*kv + b) + m2*sqrt(3)*(-1/2*kv_dot + b_dot) ...
    - m3_dot*sqrt(3)*kd/2 - m3*sqrt(3)*kd_dot/2;
tempMat2 = m3_dot*sqrt(3)*(-kd/2*mag_R/2*cot(mag_R/2)) + m3*sqrt(3)*(-kd_dot/2*mag_R/2*cot(mag_R/2))...
    + m1_dot*sqrt(3)*1/2 + m2_dot*sqrt(3)*(-1/2*kv + b) + m2*sqrt(3)*(-1/2*kv_dot + b_dot) ...
    + 2*m5*m5_dot*sqrt(3)*kd/(4*m4)*mag_R + m5^2*sqrt(3)*kd_dot/(4*m4)*mag_R;
tempMat3 = -m3_dot*sqrt(3)*(-kd/2*mag_R/2*cot(mag_R/2)) - m3*sqrt(3)*(-kd_dot/2*mag_R/2*cot(mag_R/2))...
    - m1_dot*sqrt(3)*1/2 - m2_dot*sqrt(3)*(-1/2*kv + b) - m2*sqrt(3)*(-1/2*kv_dot + b_dot)...
    + 2*m5*m5_dot*sqrt(3)*kd/(4*m4)*mag_R + m5^2*sqrt(3)*kd_dot/(4*m4)*mag_R;
D1_abs_fun{2} = {tempMat1, -tempMat1, tempMat2, tempMat3};
% D1_abs_fun{2} = {sqrt(3)*(-kd/2*m3+1/2*(m1-m2*kv)+b*m2),...
%     -sqrt(3)*(-kd/2*m3+1/2*(m1-m2*kv)+b*m2),...
%     sqrt(3)*(-kd/2*m3*mag_R/2*cot(mag_R/2)+1/2*(m1-m2*kv)+b*m2) + sqrt(3)*kd*m5^2/(4*m4)*mag_R,... % NOTE: abs(a*1i) = a for a>=0
%     -sqrt(3)*(-kd/2*m3*mag_R/2*cot(mag_R/2)+1/2*(m1-m2*kv)+b*m2) + sqrt(3)*kd*m5^2/(4*m4)*mag_R};
%
tempMat1 = m3_dot*sqrt(3)*(-1/8*mag_W^2) + 2*m5*m5_dot*sqrt(3)*(mag_W^2/(8*m4));
tempMat2 = m2_dot*sqrt(3)*(-mag_W/4) + m5_dot*m6*sqrt(3)*(mag_W/(4*m4)) + m5*m6_dot*sqrt(3)*(mag_W/(4*m4))...
    + m3_dot*sqrt(3)*(kv*mag_W/4) + m3*sqrt(3)*(kv_dot*mag_W/4) ...
    + 2*m5*m5_dot*sqrt(3)*(-mag_W*kv/(4*m4)) + m5^2*sqrt(3)*(-mag_W*kv_dot/(4*m4));
D1_abs_fun{3} = {tempMat1 + tempMat2, -tempMat1 + tempMat2,...
    tempMat1 - tempMat2, -tempMat1 - tempMat2};
% D1_abs_fun{3} = {sqrt(3)*(-1/8*(m3-m5^2/m4)*mag_W^2)+sqrt(3)*(-1/4*(m2-m5*m6/m4-m3*kv+m5^2/m4*kv)*mag_W),...
%     -sqrt(3)*(-1/8*(m3-m5^2/m4)*mag_W^2)+sqrt(3)*(-1/4*(m2-m5*m6/m4-m3*kv+m5^2/m4*kv)*mag_W),...
%     sqrt(3)*(-1/8*(m3-m5^2/m4)*mag_W^2)-sqrt(3)*(-1/4*(m2-m5*m6/m4-m3*kv+m5^2/m4*kv)*mag_W),...
%     -sqrt(3)*(-1/8*(m3-m5^2/m4)*mag_W^2)-sqrt(3)*(-1/4*(m2-m5*m6/m4-m3*kv+m5^2/m4*kv)*mag_W)};
%
tempMat1 = m2_dot*sqrt(3)*(kd/2) + m2*sqrt(3)*(kd_dot/2) ...
    + m5_dot*sqrt(3)*(-kd/2) + m5*sqrt(3)*(-kd_dot/2);
tempMat2 = m2_dot*sqrt(3)*(kd/2*mag_R/2*cot(mag_R/2)) + m2*sqrt(3)*(kd_dot/2*mag_R/2*cot(mag_R/2))...
    + m5_dot*sqrt(3)*(-kd/2*mag_R/2*cot(mag_R/2)) + m5*sqrt(3)*(-kd_dot/2*mag_R/2*cot(mag_R/2));
tempMat3 = m2_dot*sqrt(3)*(kd/2*mag_R/2) + m2*sqrt(3)*(kd_dot/2*mag_R/2)...
    + m5_dot*sqrt(3)*(-kd/2*mag_R/2) + m5*sqrt(3)*(-kd_dot/2*mag_R/2)...
    + m5_dot*m6*sqrt(3)*(-kd*mag_R/(4*m4)) + m5*m6_dot*sqrt(3)*(-kd*mag_R/(4*m4)) + m5*m6*sqrt(3)*(-kd_dot*mag_R/(4*m4));
D1_abs_fun{4} = {tempMat1, -tempMat1, tempMat2 + tempMat3,...
    tempMat2 - tempMat3, -tempMat2 + tempMat3, -tempMat2 - tempMat3};
% D1_abs_fun{4} = {sqrt(3)*(kd/2*(m2-m5)),...
%     -sqrt(3)*(kd/2*(m2-m5)),...
%     sqrt(3)*(kd/2*(m2-m5)*mag_R/2*cot(mag_R/2)) + sqrt(3)*(kd/2*(m2-m5)*mag_R/2 - m5*m6/(4*m4)*kd*mag_R),...
%     sqrt(3)*(kd/2*(m2-m5)*mag_R/2*cot(mag_R/2)) - sqrt(3)*(kd/2*(m2-m5)*mag_R/2 - m5*m6/(4*m4)*kd*mag_R),...
%     -sqrt(3)*(kd/2*(m2-m5)*mag_R/2*cot(mag_R/2)) + sqrt(3)*(kd/2*(m2-m5)*mag_R/2 - m5*m6/(4*m4)*kd*mag_R),...
%     -sqrt(3)*(kd/2*(m2-m5)*mag_R/2*cot(mag_R/2)) - sqrt(3)*(kd/2*(m2-m5)*mag_R/2 - m5*m6/(4*m4)*kd*mag_R)};
%
tempMat1 = m6_dot*sqrt(3)*(-kref/2+b) + m6*sqrt(3)*(-kref_dot/2 + b_dot);
D1_abs_fun{5} = {tempMat1, -tempMat1};
% D1_abs_fun{5} = {sqrt(3)*(-m6*kref/2+b*m6),...
%     -sqrt(3)*(-m6*kref/2+b*m6)};
%
tempMat1 = 2*m6*m6_dot*sqrt(3)*mag_W/(4*m4) + m5_dot*m6*sqrt(3)*(-mag_W*kv/(4*m4))...
    + m5*m6_dot*sqrt(3)*(-mag_W*kv/(4*m4)) + m5*m6*sqrt(3)*(-mag_W*kv_dot/(4*m4));
D1_abs_fun{6} = {tempMat1, -tempMat1};
% D1_abs_fun{6} = {sqrt(3)*((m6^2-m5*m6*kv)/(4*m4)*mag_W),...
%     -sqrt(3)*((m6^2-m5*m6*kv)/(4*m4)*mag_W)};

% For disc centered at D(2,2)
D2_abs_fun{1} = D1_abs_fun{2};
D2_abs_fun{2} = D1_abs_fun{3};
%
tempMat1 = m3_dot*sqrt(3)*kd/2 + m3*sqrt(3)*kd_dot/2 ...
    + m6_dot*sqrt(3)*1/2 + m5_dot*sqrt(3)*(-kv/2) + m5*sqrt(3)*(-kv_dot/2);
tempMat2 = m3_dot*sqrt(3)*kd/2*mag_R/2*cot(mag_R/2) + m3*sqrt(3)*kd_dot/2*mag_R/2*cot(mag_R/2)...
    + m6_dot*sqrt(3)*1/2 ...
    + m5_dot*sqrt(3)*(-kv/2) + m5*sqrt(3)*(-kv_dot/2);
tempMat3 = m3_dot*sqrt(3)*(kd/2*mag_R/2) + m3*sqrt(3)*(kd_dot/2*mag_R/2)...
    + 2*m5*m5_dot*sqrt(3)*(-kd/(4*m4)*mag_R) + m5^2*sqrt(3)*(-kd_dot/(4*m4)*mag_R);
D2_abs_fun{3} = {tempMat1, -tempMat1, tempMat2 + tempMat3,...
    tempMat2 - tempMat3, -tempMat2 + tempMat3, -tempMat2 - tempMat3};
% D2_abs_fun{3} = {sqrt(3)*(m3*kd/2+(m6-m5*kv)/2),...
%     -sqrt(3)*(m3*kd/2+(m6-m5*kv)/2),...
%     sqrt(3)*(m3*kd/2*mag_R/2*cot(mag_R/2)+(m6-m5*kv)/2) + sqrt(3)*(m3*kd/2*mag_R/2-m5^2*kd/(4*m4)*mag_R),...
%     sqrt(3)*(m3*kd/2*mag_R/2*cot(mag_R/2)+(m6-m5*kv)/2) - sqrt(3)*(m3*kd/2*mag_R/2-m5^2*kd/(4*m4)*mag_R),...
%     -sqrt(3)*(m3*kd/2*mag_R/2*cot(mag_R/2)+(m6-m5*kv)/2) + sqrt(3)*(m3*kd/2*mag_R/2-m5^2*kd/(4*m4)*mag_R),...
%     -sqrt(3)*(m3*kd/2*mag_R/2*cot(mag_R/2)+(m6-m5*kv)/2) - sqrt(3)*(m3*kd/2*mag_R/2-m5^2*kd/(4*m4)*mag_R)};
%
tempMat1 = m5_dot*(-kref/2 + b) + m5*(-kref_dot/2 + b_dot);
D2_abs_fun{4} = {tempMat1, -tempMat1};
% D2_abs_fun{4} = {-m5*kref/2+b*m5,...
%     -(-m5*kref/2+b*m5)};
%
tempMat1 = m5_dot*m6*sqrt(3)*mag_W/(4*m4) + m5*m6_dot*sqrt(3)*mag_W/(4*m4)...
    + 2*m5*m5_dot*sqrt(3)*(-kv*mag_W/(4*m4)) + m5^2*sqrt(3)*(-kv_dot*mag_W/(4*m4));
D2_abs_fun{5} = {tempMat1, -tempMat1};
% D2_abs_fun{5} = {sqrt(3)*( (m5*m6-m5^2*kv)/(4*m4)*mag_W ),...
%     -sqrt(3)*( (m5*m6-m5^2*kv)/(4*m4)*mag_W )};

% For disc centered at D(3,3)
D3_abs_fun{1} = D1_abs_fun{4};
D3_abs_fun{2} = D1_abs_fun{5};
D3_abs_fun{3} = D1_abs_fun{6};
D3_abs_fun{4} = D2_abs_fun{3};
D3_abs_fun{5} = D2_abs_fun{4};
D3_abs_fun{6} = D2_abs_fun{5};
%
tempMat1 = m5_dot*sqrt(3)*kd + m5*sqrt(3)*kd_dot;
D3_abs_fun{7} = {tempMat1, -tempMat1};
% D3_abs_fun{7} = {sqrt(3)*m5*kd,...
%     -sqrt(3)*m5*kd};
%% Permutate all possible matrices
% For disc centered at D(1,1)
flagUseRepmat = true;
for ii = 1:length(D1_abs_fun)
    tempA = D1_abs_fun{ii}; % this should now be a vector of cells of functions
    len_A = length(tempA);
    % expand the D1_Bounds by length of A
    % NOTE: repelem repeats elements consecutivly
    if flagUseRepmat
        % For some reason, repelem does not like expanding the cell when
        % there is only on cxv function, this bypasses that
        A1_mat = repmat(A1_mat,1,len_A);
        flagUseRepmat = false;
    else
        A1_mat = repelem(A1_mat,len_A);
    end
    % Add each element of A to each element of D1_Bounds
    for jj = 1:length(A1_mat)
        A1_mat{jj} = A1_mat{jj} + tempA{mod(jj,len_A)+1};
    end
end
% For disc centered at D(2,2)
flagUseRepmat = true;
for ii = 1:length(D2_abs_fun)
    tempA = D2_abs_fun{ii}; % this should now be a vector of cells of functions
    len_A = length(tempA);
    % expand the D1_Bounds by length of A
    % NOTE: repelem repeats elements consecutivly
    if flagUseRepmat
        % For some reason, repelem does not like expanding the cell when
        % there is only on cxv function, this bypasses that
        A2_mat = repmat(A2_mat,1,len_A);
        flagUseRepmat = false;
    else
        A2_mat = repelem(A2_mat,len_A);
    end
    % Add each element of A to each element of D1_Bounds
    for jj = 1:length(A2_mat)
        A2_mat{jj} = A2_mat{jj} + tempA{mod(jj,len_A)+1};
    end
end
% For disc centered at D(3,3)
flagUseRepmat = true;
for ii = 1:length(D3_abs_fun)
    tempA = D3_abs_fun{ii}; % this should now be a vector of cells of functions
    len_A = length(tempA);
    % expand the D1_Bounds by length of A
    % NOTE: repelem repeats elements consecutivly
    if flagUseRepmat
        % For some reason, repelem does not like expanding the cell when
        % there is only on cxv function, this bypasses that
        A3_mat = repmat(A3_mat,1,len_A);
        flagUseRepmat = false;
    else
        A3_mat = repelem(A3_mat,len_A);
    end
    % Add each element of A to each element of D1_Bounds
    for jj = 1:length(A3_mat)
        A3_mat{jj} = A3_mat{jj} + tempA{mod(jj,len_A)+1};
    end
end
%% Results
A_all = [A1_mat A2_mat A3_mat];
% Convert A_all to function handles
tempBounds = cvx(zeros(length(A_all),1));
for ii = 1:length(A_all)
    % Compute the max gershgorin bound value as a scalar expression
    tempBounds(ii) = A_all{ii};
end

A_all_der_func = tempBounds;
end