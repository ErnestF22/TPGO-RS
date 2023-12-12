function [A_all_func_kdot] = rotRefOpt_allBounds_Matrix_kdot()
% Last Edited: Jan 6, 2021 by Bee Vang
% Compute the relaxed gershgorin disc bounds derivative wrt the gains 
% [kd;kv;kref] as real functions using
% trace(Ai*X)<=0 where X = x*x' and x = [1;m1;m2;m3;m5;m6] as function
% handles  of (R,w,RRef,gains_dot,M_nn,beta)
% NOTE: To do this we relaxed abs(a+bi) <= abs(a) + abs(b)
% INPUTS:
% OUTPUTS:
%   A_all_func_kdot(R,w,RRef,gains_dot,M_nn,beta) := a cell array of all possible 
%       bounds on the contraction matrix eigenvalue on TSO(3)xSO(3)

%% Default Parameters
flagSymbolic = true;
syms kd_dot kv_dot kref_dot b m1 m2 m3 m4 m5 m6 mag_R mag_RRef mag_W real;
%% Init matrices
% NOTE: we are computing [1 m1 m2 m3 m5 m6]*A*[[1 m1 m2 m3 m5 m6], but only
% the submatrix A(2:end,2:end) needs to be PSD which reduces to just the
% [2x2] matrix given as [m5 m6]*Bi*[m5;m6]
% Init matrix based on symbolic or not
A_11 = sym(zeros(6,6)); A_22 = A_11; A_33 = A_11;
% Populate the centroid terms
A_11(1,3) = -kd_dot*mag_R/2*cot(mag_R/2); A_11 = (A_11 + A_11')/2;
A1_mat = {A_11};
A_22(1,3) = 0; A_22(1,4) = -kv_dot; A_22 = (A_22 + A_22')/2;
A2_mat = {A_22};
A_33(1,1) = -kref_dot*m4*mag_RRef/2*cot(mag_RRef/2);
A_33 = (A_33 + A_33')/2;
A3_mat = {A_33};

% Define possible terms resulting from absolute values
% For disc centered at D(1,1)
tempMat1 = create_empty_mat(flagSymbolic,6);
tempMat2 = create_empty_mat(flagSymbolic,6);
% Fill entries
tempMat1(1,2) = 0;
tempMat1 = (tempMat1 + tempMat1')/2;
tempMat2(1,2) = 0;
tempMat2(1,3) = 0;
tempMat2(5,6) = 0;
tempMat2 = (tempMat2 + tempMat2')/2;
D1_abs_fun{1} = {tempMat1,tempMat2,-tempMat2};
% D1_abs_fun{1} = {sqrt(3)*b*m1,...
%     -sqrt(3)*b*m1,...
%     sqrt(3)*(-1/4*(m2-m5*m6/m4)*mag_W^2+b*m1),...
%     -sqrt(3)*(-1/4*(m2-m5*m6/m4)*mag_W^2+b*m1)};
% 
tempMat1 = create_empty_mat(flagSymbolic,6);
tempMat2 = create_empty_mat(flagSymbolic,6);
tempMat3 = create_empty_mat(flagSymbolic,6);
% Fill entries
tempMat1(1,2) = 0;
tempMat1(1,3) = sqrt(3)*(-1/2*kv_dot);
tempMat1(1,4) = -sqrt(3)*kd_dot/2;
tempMat1 = (tempMat1 + tempMat1')/2;
tempMat2(1,4) = sqrt(3)*(-kd_dot/2*mag_R/2*cot(mag_R/2));
tempMat2(1,2) = 0;
tempMat2(1,3) = sqrt(3)*(-1/2*kv_dot);
tempMat2(5,5) = sqrt(3)*kd_dot/(4*m4)*mag_R;
tempMat2 = (tempMat2 + tempMat2')/2;
tempMat3(1,4) = -sqrt(3)*(-kd_dot/2*mag_R/2*cot(mag_R/2));
tempMat3(1,2) = 0;
tempMat3(1,3) = -sqrt(3)*(-1/2*kv_dot);
tempMat3(5,5) = sqrt(3)*kd_dot/(4*m4)*mag_R;
tempMat3 = (tempMat3 + tempMat3')/2;
D1_abs_fun{2} = {tempMat1, -tempMat1, tempMat2, tempMat3};
% D1_abs_fun{2} = {sqrt(3)*(-kd/2*m3+1/2*(m1-m2*kv)+b*m2),...
%     -sqrt(3)*(-kd/2*m3+1/2*(m1-m2*kv)+b*m2),...
%     sqrt(3)*(-kd/2*m3*mag_R/2*cot(mag_R/2)+1/2*(m1-m2*kv)+b*m2) + sqrt(3)*kd*m5^2/(4*m4)*mag_R,... % NOTE: abs(a*1i) = a for a>=0
%     -sqrt(3)*(-kd/2*m3*mag_R/2*cot(mag_R/2)+1/2*(m1-m2*kv)+b*m2) + sqrt(3)*kd*m5^2/(4*m4)*mag_R};
%
tempMat1 = create_empty_mat(flagSymbolic,6);
tempMat2 = create_empty_mat(flagSymbolic,6);
tempMat1(1,4) = 0;
tempMat1(5,5) = 0;
tempMat1 = (tempMat1 + tempMat1')/2;
tempMat2(1,3) = 0;
tempMat2(5,6) = 0;
tempMat2(1,4) = sqrt(3)*(kv_dot*mag_W/4);
tempMat2(5,5) = sqrt(3)*(-mag_W*kv_dot/(4*m4));
tempMat2 = (tempMat2 + tempMat2')/2;
D1_abs_fun{3} = {tempMat1 + tempMat2, -tempMat1 + tempMat2,...
    tempMat1 - tempMat2, -tempMat1 - tempMat2};
% D1_abs_fun{3} = {sqrt(3)*(-1/8*(m3-m5^2/m4)*mag_W^2)+sqrt(3)*(-1/4*(m2-m5*m6/m4-m3*kv+m5^2/m4*kv)*mag_W),...
%     -sqrt(3)*(-1/8*(m3-m5^2/m4)*mag_W^2)+sqrt(3)*(-1/4*(m2-m5*m6/m4-m3*kv+m5^2/m4*kv)*mag_W),...
%     sqrt(3)*(-1/8*(m3-m5^2/m4)*mag_W^2)-sqrt(3)*(-1/4*(m2-m5*m6/m4-m3*kv+m5^2/m4*kv)*mag_W),...
%     -sqrt(3)*(-1/8*(m3-m5^2/m4)*mag_W^2)-sqrt(3)*(-1/4*(m2-m5*m6/m4-m3*kv+m5^2/m4*kv)*mag_W)};
%
tempMat1 = create_empty_mat(flagSymbolic,6);
tempMat2 = create_empty_mat(flagSymbolic,6);
tempMat3 = create_empty_mat(flagSymbolic,6);
%
tempMat1(1,3) = sqrt(3)*(kd_dot/2);
tempMat1(1,5) = sqrt(3)*(-kd_dot/2);
tempMat1 = (tempMat1 + tempMat1')/2;
tempMat2(1,3) = sqrt(3)*(kd_dot/2*mag_R/2*cot(mag_R/2));
tempMat2(1,5) = sqrt(3)*(-kd_dot/2*mag_R/2*cot(mag_R/2));
tempMat2 = (tempMat2 + tempMat2')/2;
tempMat3(1,3) = sqrt(3)*(kd_dot/2*mag_R/2);
tempMat3(1,5) = sqrt(3)*(-kd_dot/2*mag_R/2);
tempMat3(5,6) = sqrt(3)*(-kd_dot*mag_R/(4*m4));
tempMat3 = (tempMat3 + tempMat3')/2;
D1_abs_fun{4} = {tempMat1, -tempMat1, tempMat2 + tempMat3,...
    tempMat2 - tempMat3, -tempMat2 + tempMat3, -tempMat2 - tempMat3};
% D1_abs_fun{4} = {sqrt(3)*(kd/2*(m2-m5)),...
%     -sqrt(3)*(kd/2*(m2-m5)),...
%     sqrt(3)*(kd/2*(m2-m5)*mag_R/2*cot(mag_R/2)) + sqrt(3)*(kd/2*(m2-m5)*mag_R/2 - m5*m6/(4*m4)*kd*mag_R),...
%     sqrt(3)*(kd/2*(m2-m5)*mag_R/2*cot(mag_R/2)) - sqrt(3)*(kd/2*(m2-m5)*mag_R/2 - m5*m6/(4*m4)*kd*mag_R),...
%     -sqrt(3)*(kd/2*(m2-m5)*mag_R/2*cot(mag_R/2)) + sqrt(3)*(kd/2*(m2-m5)*mag_R/2 - m5*m6/(4*m4)*kd*mag_R),...
%     -sqrt(3)*(kd/2*(m2-m5)*mag_R/2*cot(mag_R/2)) - sqrt(3)*(kd/2*(m2-m5)*mag_R/2 - m5*m6/(4*m4)*kd*mag_R)};
%
tempMat1 = create_empty_mat(flagSymbolic,6);
tempMat1(1,6) = sqrt(3)*(-kref_dot/2); tempMat1 = (tempMat1 + tempMat1')/2;
D1_abs_fun{5} = {tempMat1, -tempMat1};
% D1_abs_fun{5} = {sqrt(3)*(-m6*kref/2+b*m6),...
%     -sqrt(3)*(-m6*kref/2+b*m6)};
%
tempMat1 = create_empty_mat(flagSymbolic,6);
tempMat1(6,6) = 0;
tempMat1(5,6) = sqrt(3)*(-mag_W*kv_dot/(4*m4)); tempMat1 = (tempMat1+tempMat1')/2;
D1_abs_fun{6} = {tempMat1, -tempMat1};
% D1_abs_fun{6} = {sqrt(3)*((m6^2-m5*m6*kv)/(4*m4)*mag_W),...
%     -sqrt(3)*((m6^2-m5*m6*kv)/(4*m4)*mag_W)};

% For disc centered at D(2,2)
D2_abs_fun{1} = D1_abs_fun{2};
D2_abs_fun{2} = D1_abs_fun{3};
%
tempMat1 = create_empty_mat(flagSymbolic,6);
tempMat2 = create_empty_mat(flagSymbolic,6);
tempMat3 = create_empty_mat(flagSymbolic,6);
%
tempMat1(1,4) = sqrt(3)*kd_dot/2;
tempMat1(1,6) = 0;
tempMat1(1,5) = sqrt(3)*(-kv_dot/2); tempMat1 = (tempMat1 + tempMat1')/2;
tempMat2(1,4) = sqrt(3)*kd_dot/2*mag_R/2*cot(mag_R/2);
tempMat2(1,6) = 0;
tempMat2(1,5) = sqrt(3)*(-kv_dot/2); tempMat2 = (tempMat2 + tempMat2')/2;
tempMat3(1,4) = sqrt(3)*(kd_dot/2*mag_R/2);
tempMat3(5,5) = sqrt(3)*(-kd_dot/(4*m4)*mag_R); tempMat3 = (tempMat3 + tempMat3')/2;
D2_abs_fun{3} = {tempMat1, -tempMat1, tempMat2 + tempMat3,...
    tempMat2 - tempMat3, -tempMat2 + tempMat3, -tempMat2 - tempMat3};
% D2_abs_fun{3} = {sqrt(3)*(m3*kd/2+(m6-m5*kv)/2),...
%     -sqrt(3)*(m3*kd/2+(m6-m5*kv)/2),...
%     sqrt(3)*(m3*kd/2*mag_R/2*cot(mag_R/2)+(m6-m5*kv)/2) + sqrt(3)*(m3*kd/2*mag_R/2-m5^2*kd/(4*m4)*mag_R),...
%     sqrt(3)*(m3*kd/2*mag_R/2*cot(mag_R/2)+(m6-m5*kv)/2) - sqrt(3)*(m3*kd/2*mag_R/2-m5^2*kd/(4*m4)*mag_R),...
%     -sqrt(3)*(m3*kd/2*mag_R/2*cot(mag_R/2)+(m6-m5*kv)/2) + sqrt(3)*(m3*kd/2*mag_R/2-m5^2*kd/(4*m4)*mag_R),...
%     -sqrt(3)*(m3*kd/2*mag_R/2*cot(mag_R/2)+(m6-m5*kv)/2) - sqrt(3)*(m3*kd/2*mag_R/2-m5^2*kd/(4*m4)*mag_R)};
%
tempMat1 = create_empty_mat(flagSymbolic,6);
tempMat1(1,5) = (-kref_dot/2); tempMat1 = (tempMat1 + tempMat1')/2;
D2_abs_fun{4} = {tempMat1, -tempMat1};
% D2_abs_fun{4} = {-m5*kref/2+b*m5,...
%     -(-m5*kref/2+b*m5)};
%
tempMat1 = create_empty_mat(flagSymbolic,6);
tempMat1(5,6) = 0;
tempMat1(5,5) = sqrt(3)*(-kv_dot*mag_W/(4*m4)); tempMat1 = (tempMat1 + tempMat1')/2;
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
tempMat1 = create_empty_mat(flagSymbolic,6);
tempMat1(1,5) = sqrt(3)*kd_dot; tempMat1 = (tempMat1 + tempMat1')/2;
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
x = [1;m1;m2;m3;m5;m6]; X = x*x';
for ii = 1:length(A_all)
    temp_sym = trace(A_all{ii}*X);
    tempBounds(ii) = temp_sym;
end

% Convert symbolic vector to function handle
temp_func = matlabFunction(tempBounds,'Vars',[mag_R,mag_RRef,kd_dot,kv_dot,kref_dot,mag_W,m1,m2,m3,m4,m5,m6,b]);
% For convenience, remap as function of (R,w,RRef,[gains],M_nn,b)
A_all_func_kdot = @(R,w,RRef,gains_dot,M_nn,beta) temp_func( norm(rot_log(R,RRef)),... % mag_R
    norm(rot_log(RRef,eye(3))),... % mag_RRef
    gains_dot(1),.... % kd_dot
    gains_dot(2),.... % kv_dot
    gains_dot(3),.... % kref_dot
    norm(w),.... % mag_W
    M_nn(1,1),... % m1
    M_nn(1,2),... % m2
    M_nn(2,2),... % m3
    M_nn(3,3),... % m4
    M_nn(2,3),... % m5
    M_nn(1,3),... % m6
    beta); % beta
end