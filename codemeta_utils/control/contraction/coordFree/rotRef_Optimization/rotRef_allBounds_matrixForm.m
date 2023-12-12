function [] = rotRef_allBounds_matrixForm()
% Last Edited: Jan 28, 2021 by Bee Vang
% Same as rotRefOpt_allBounds_Matrix_simple.m, but generate function files
% for use in rotRef_contractionOpt_new.m.
% NOTE: To do this we relaxed abs(a+bi) <= abs(a) + abs(b)
% INPUTS:
% OUTPUTS:
% GEN FILES:
%   

%% Default Parameters
syms kd kv kref b m1 m2 m3 m4 m5 m6 mag_R mag_RRef mag_W real;
syms m5m6 m5_squared m6_squared real
%% Init matrices
% NOTE: we are computing [1 m1 m2 m3 m5 m6]*A*[[1 m1 m2 m3 m5 m6], but only
% the submatrix A(2:end,2:end) needs to be PSD which reduces to just the
% [2x2] matrix given as [m5 m6]*Bi*[m5;m6]
% Init matrix based on symbolic or not
% Populate the centroid terms
A_11 = -m2*kd*mag_R/2*cot(mag_R/2);
A1_mat = {A_11};
A_22 = m2 + m3*(b-kv);
A2_mat = {A_22};
A_33 = m4*b-kref*m4*mag_RRef/2*cot(mag_RRef/2);
A3_mat = {A_33};

% Define possible terms resulting from absolute values
% For disc centered at D(1,1)
% Fill entries
tempMat1 = m1*b;
tempMat2 = m1*b + m2*(1/4*mag_W^2)...
    - m5m6*(mag_W^2/(4*m4));
D1_abs_fun{1} = {tempMat1};
% D1_abs_fun{1} = {sqrt(3)*b*m1,...
%     -sqrt(3)*b*m1,...
%     sqrt(3)*(-1/4*(m2-m5*m6/m4)*mag_W^2+b*m1),...
%     -sqrt(3)*(-1/4*(m2-m5*m6/m4)*mag_W^2+b*m1)};
% 
% Fill entries
tempMat1 = m1*sqrt(3)*1/2 + m2*sqrt(3)*(-1/2*kv + b) ...
    - m3*sqrt(3)*kd/2;
tempMat2 = m3*sqrt(3)*(-kd/2*mag_R/2*cot(mag_R/2))...
    + m1*sqrt(3)*1/2 + m2*sqrt(3)*(-1/2*kv + b) ...
    + m5_squared*sqrt(3)*kd/(4*m4)*mag_R;
tempMat3 = -m3*sqrt(3)*(-kd/2*mag_R/2*cot(mag_R/2))...
    - m1*sqrt(3)*1/2 - m2*sqrt(3)*(-1/2*kv + b)...
    + m5_squared*sqrt(3)*kd/(4*m4)*mag_R;
D1_abs_fun{2} = {tempMat1, -tempMat1, tempMat2, tempMat3};
% D1_abs_fun{2} = {sqrt(3)*(-kd/2*m3+1/2*(m1-m2*kv)+b*m2),...
%     -sqrt(3)*(-kd/2*m3+1/2*(m1-m2*kv)+b*m2),...
%     sqrt(3)*(-kd/2*m3*mag_R/2*cot(mag_R/2)+1/2*(m1-m2*kv)+b*m2) + sqrt(3)*kd*m5^2/(4*m4)*mag_R,... % NOTE: abs(a*1i) = a for a>=0
%     -sqrt(3)*(-kd/2*m3*mag_R/2*cot(mag_R/2)+1/2*(m1-m2*kv)+b*m2) + sqrt(3)*kd*m5^2/(4*m4)*mag_R};
%
tempMat1 = m3*sqrt(3)*(-1/8*mag_W^2) + m5_squared*sqrt(3)*(mag_W^2/(8*m4));
tempMat2 = m2*sqrt(3)*(-mag_W/4) + m5m6*sqrt(3)*(mag_W/(4*m4)) ...
    + m3*sqrt(3)*(kv*mag_W/4) + m5_squared*sqrt(3)*(-mag_W*kv/(4*m4));
D1_abs_fun{3} = {tempMat1 + tempMat2, -tempMat1 + tempMat2,...
    tempMat1 - tempMat2, -tempMat1 - tempMat2};
% D1_abs_fun{3} = {sqrt(3)*(-1/8*(m3-m5^2/m4)*mag_W^2)+sqrt(3)*(-1/4*(m2-m5*m6/m4-m3*kv+m5^2/m4*kv)*mag_W),...
%     -sqrt(3)*(-1/8*(m3-m5^2/m4)*mag_W^2)+sqrt(3)*(-1/4*(m2-m5*m6/m4-m3*kv+m5^2/m4*kv)*mag_W),...
%     sqrt(3)*(-1/8*(m3-m5^2/m4)*mag_W^2)-sqrt(3)*(-1/4*(m2-m5*m6/m4-m3*kv+m5^2/m4*kv)*mag_W),...
%     -sqrt(3)*(-1/8*(m3-m5^2/m4)*mag_W^2)-sqrt(3)*(-1/4*(m2-m5*m6/m4-m3*kv+m5^2/m4*kv)*mag_W)};
%
tempMat1 = m2*sqrt(3)*(kd/2) + m5*sqrt(3)*(-kd/2);
tempMat2 = m2*sqrt(3)*(kd/2*mag_R/2*cot(mag_R/2))...
    + m5*sqrt(3)*(-kd/2*mag_R/2*cot(mag_R/2));
tempMat3 = m2*sqrt(3)*(kd/2*mag_R/2) + m5*sqrt(3)*(-kd/2*mag_R/2)...
    + m5m6*sqrt(3)*(-kd*mag_R/(4*m4));
D1_abs_fun{4} = {tempMat1, -tempMat1, tempMat2 + tempMat3,...
    tempMat2 - tempMat3, -tempMat2 + tempMat3, -tempMat2 - tempMat3};
% D1_abs_fun{4} = {sqrt(3)*(kd/2*(m2-m5)),...
%     -sqrt(3)*(kd/2*(m2-m5)),...
%     sqrt(3)*(kd/2*(m2-m5)*mag_R/2*cot(mag_R/2)) + sqrt(3)*(kd/2*(m2-m5)*mag_R/2 - m5*m6/(4*m4)*kd*mag_R),...
%     sqrt(3)*(kd/2*(m2-m5)*mag_R/2*cot(mag_R/2)) - sqrt(3)*(kd/2*(m2-m5)*mag_R/2 - m5*m6/(4*m4)*kd*mag_R),...
%     -sqrt(3)*(kd/2*(m2-m5)*mag_R/2*cot(mag_R/2)) + sqrt(3)*(kd/2*(m2-m5)*mag_R/2 - m5*m6/(4*m4)*kd*mag_R),...
%     -sqrt(3)*(kd/2*(m2-m5)*mag_R/2*cot(mag_R/2)) - sqrt(3)*(kd/2*(m2-m5)*mag_R/2 - m5*m6/(4*m4)*kd*mag_R)};
%
tempMat1 = m6*sqrt(3)*(-kref/2+b);
D1_abs_fun{5} = {tempMat1, -tempMat1};
% D1_abs_fun{5} = {sqrt(3)*(-m6*kref/2+b*m6),...
%     -sqrt(3)*(-m6*kref/2+b*m6)};
%
tempMat1 = m6_squared*sqrt(3)*mag_W/(4*m4) + m5m6*sqrt(3)*(-mag_W*kv/(4*m4));
D1_abs_fun{6} = {tempMat1, -tempMat1};
% D1_abs_fun{6} = {sqrt(3)*((m6^2-m5*m6*kv)/(4*m4)*mag_W),...
%     -sqrt(3)*((m6^2-m5*m6*kv)/(4*m4)*mag_W)};

% For disc centered at D(2,2)
D2_abs_fun{1} = D1_abs_fun{2};
D2_abs_fun{2} = D1_abs_fun{3};
%
tempMat1 = m3*sqrt(3)*kd/2 + m6*sqrt(3)*1/2 + m5*sqrt(3)*(-kv/2);
tempMat2 = m3*sqrt(3)*kd/2*mag_R/2*cot(mag_R/2) + m6*sqrt(3)*1/2 ...
    + m5*sqrt(3)*(-kv/2);
tempMat3 = m3*sqrt(3)*(kd/2*mag_R/2) - m5_squared*sqrt(3)*(kd/(4*m4)*mag_R);
D2_abs_fun{3} = {tempMat1, -tempMat1, tempMat2 + tempMat3,...
    tempMat2 - tempMat3, -tempMat2 + tempMat3, -tempMat2 - tempMat3};
% D2_abs_fun{3} = {sqrt(3)*(m3*kd/2+(m6-m5*kv)/2),...
%     -sqrt(3)*(m3*kd/2+(m6-m5*kv)/2),...
%     sqrt(3)*(m3*kd/2*mag_R/2*cot(mag_R/2)+(m6-m5*kv)/2) + sqrt(3)*(m3*kd/2*mag_R/2-m5^2*kd/(4*m4)*mag_R),...
%     sqrt(3)*(m3*kd/2*mag_R/2*cot(mag_R/2)+(m6-m5*kv)/2) - sqrt(3)*(m3*kd/2*mag_R/2-m5^2*kd/(4*m4)*mag_R),...
%     -sqrt(3)*(m3*kd/2*mag_R/2*cot(mag_R/2)+(m6-m5*kv)/2) + sqrt(3)*(m3*kd/2*mag_R/2-m5^2*kd/(4*m4)*mag_R),...
%     -sqrt(3)*(m3*kd/2*mag_R/2*cot(mag_R/2)+(m6-m5*kv)/2) - sqrt(3)*(m3*kd/2*mag_R/2-m5^2*kd/(4*m4)*mag_R)};
%
tempMat1 = m5*(-kref/2 + b);
tempMat2 = m5*(-kref/2*mag_RRef/2*cot(mag_RRef/2) + b);
D2_abs_fun{4} = {tempMat1, -tempMat1, tempMat2, -tempMat2};
% D2_abs_fun{4} = {-m5*kref/2+b*m5,...
%     -(-m5*kref/2+b*m5)};
%
tempMat1 = m5m6*sqrt(3)*mag_W/(4*m4) + m5_squared*sqrt(3)*(-kv*mag_W/(4*m4));
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
tempMat1 = m5*sqrt(3)*kd;
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

% Extract the bounds as matrix product of
% A*[m1;m2;m3;m5;m6;m5m6;m5_squared;m6_squared]-B=tempBounds
% NOTE: m4 is considered a constant so it'll show up in B
[A,B] = equationsToMatrix(tempBounds,...
    [m1;m2;m3;m5;m6;m5m6;m5_squared;m6_squared]);

% Convert the "constant" matrices to function files
temp_A_func = matlabFunction(A,'Vars',[mag_R,mag_RRef,kd,kv,kref,mag_W,m4,b],'File','rotRef_gen_allBounds_A');
temp_B_func = matlabFunction(B,'Vars',[mag_R,mag_RRef,kd,kv,kref,mag_W,m4,b],'File','rotRef_gen_allBounds_B');
end