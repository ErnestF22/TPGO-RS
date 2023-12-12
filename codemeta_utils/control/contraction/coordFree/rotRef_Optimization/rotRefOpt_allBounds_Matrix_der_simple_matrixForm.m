function [A,B] = rotRefOpt_allBounds_Matrix_der_simple_matrixForm()
% Last Edited: Jan 8, 2021 by Bee Vang
% Compute the relaxed gershgorin disc bounds' derivative wrt the gains,
% metric, and convergence rate. Generate a function handle of
% (R,w,RRef,gains,M_nn,beta) that creates the derivative 
% bounds as A*x+b where x = [gains_dot;M_nn_dot;bdot]. The results are
% numerics that can be used inside cvx where x is a cvx variable.
% NOTE: To do this we relaxed abs(a+bi) <= abs(a) + abs(b)
% NOTE: Also create optimized function files that are faster to use than
% function handles
% INPUTS:
% OUTPUTS:
%   [A,B](R,w,RRef,gains,M_nn,beta) := A pair of function handles of all possible 
%       derivative of bounds on the contraction matrix eigenvalue on TSO(3)xSO(3)
%       wrt the gains, metric tensor, and contraction rate in the form of
%       A*x - B where x = [gains_dot;M_nn_dot;bdot]

%% Default Parameters
syms mag_R mag_RRef mag_W kd kv kref m1 m2 m3 m4 m5 m6 b real
syms kd_dot kv_dot kref_dot m1_dot m2_dot m3_dot m4_dot m5_dot m6_dot b_dot real

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
for ii = 1:length(D1_abs_fun)
    tempA = D1_abs_fun{ii}; % this should now be a vector of cells of functions
    len_A = length(tempA);
    % expand the D1_Bounds by length of A
    % NOTE: repelem repeats elements consecutivly
    A1_mat = repelem(A1_mat,len_A);
    % Add each element of A to each element of D1_Bounds
    for jj = 1:length(A1_mat)
        A1_mat{jj} = A1_mat{jj} + tempA{mod(jj,len_A)+1};
    end
end
% For disc centered at D(2,2)
for ii = 1:length(D2_abs_fun)
    tempA = D2_abs_fun{ii}; % this should now be a vector of cells of functions
    len_A = length(tempA);
    % expand the D1_Bounds by length of A
    % NOTE: repelem repeats elements consecutivly    
    A2_mat = repelem(A2_mat,len_A);
    % Add each element of A to each element of D1_Bounds
    for jj = 1:length(A2_mat)
        A2_mat{jj} = A2_mat{jj} + tempA{mod(jj,len_A)+1};
    end
end
% For disc centered at D(3,3)
for ii = 1:length(D3_abs_fun)
    tempA = D3_abs_fun{ii}; % this should now be a vector of cells of functions
    len_A = length(tempA);
    % expand the D1_Bounds by length of A
    % NOTE: repelem repeats elements consecutivly
    A3_mat = repelem(A3_mat,len_A);
    % Add each element of A to each element of D1_Bounds
    for jj = 1:length(A3_mat)
        A3_mat{jj} = A3_mat{jj} + tempA{mod(jj,len_A)+1};
    end
end
%% Results
A_all = [A1_mat A2_mat A3_mat];
% Convert A_all to function handles
tempBounds = sym(zeros(length(A_all),1));
for ii = 1:length(A_all)
    % Compute the max gershgorin bound value as a scalar expression
    tempBounds(ii) = A_all{ii};
end

A_all_der_sym = tempBounds;
% Find A*x - b = A_all_der_sym, where x = [gains_dot;M_nn_dot;bdot]
[A,B] = equationsToMatrix(A_all_der_sym,...
    [kd_dot;kv_dot;kref_dot;...
    m1_dot;m2_dot;m3_dot;m4_dot;m5_dot;m6_dot;...
    b_dot]);

% Remap [A,b] to be function handles of (R,w,RRef,gains,M_nn,b)
temp_A_func = matlabFunction(A,'Vars',[mag_R,mag_RRef,kd,kv,kref,mag_W,m1,m2,m3,m4,m5,m6,b],'File','rotRefOpt_gen_allBounds_der_A_mat');
temp_B_func = matlabFunction(B,'Vars',[mag_R,mag_RRef,kd,kv,kref,mag_W,m1,m2,m3,m4,m5,m6,b],'File','rotRefOpt_gen_allBounds_der_B_mat');
A = @(R,w,RRef,gains,M_nn,beta) temp_A_func( norm(rot_log(R,RRef)),...
    norm(rot_log(RRef,eye(3))),... % mag_RRef
    gains(1),.... % kd
    gains(2),.... % kv
    gains(3),.... % kref
    norm(w),.... % mag_W
    M_nn(1,1),... % m1
    M_nn(1,2),... % m2
    M_nn(2,2),... % m3
    M_nn(3,3),... % m4
    M_nn(2,3),... % m5
    M_nn(1,3),... % m6
    beta); % beta
B = @(R,w,RRef,gains,M_nn,beta) temp_B_func( norm(rot_log(R,RRef)),...
    norm(rot_log(RRef,eye(3))),... % mag_RRef
    gains(1),.... % kd
    gains(2),.... % kv
    gains(3),.... % kref
    norm(w),.... % mag_W
    M_nn(1,1),... % m1
    M_nn(1,2),... % m2
    M_nn(2,2),... % m3
    M_nn(3,3),... % m4
    M_nn(2,3),... % m5
    M_nn(1,3),... % m6
    beta); % beta
end