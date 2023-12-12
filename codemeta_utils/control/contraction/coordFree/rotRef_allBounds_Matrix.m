function [A_all] = rotRef_allBounds_Matrix(mag_R,mag_RRef,kd,kv,kp,mag_W,beta,varargin)
% LAST EDITED: Jan 29, 2021 by Bee Vang.
% NOTE: The bounds here are wrong, the centroid for the top 3 disc had
% additional negative terms, in addition some simplifications were made
% (see rotRefOpt_allBounds_Matrix_simple.m for updated and more correct
% bounds. Leaving this file sicne some other functions may depend on it.
%
% Compute the relaxed gershgorin disc bounds as matrices in the form of
% trace(Ai*X)<=0 where X = x*x' and x = [1;m1;m2;m3;m5;m6]
% NOTE: To do this we relaxed abs(a+bi) <= abs(a) + abs(b)
% INPUTS:
%   mag_R := positive scalar representing max distance from R to RRef
%   mag_RRef := positive scalar representing max distance from RRef to I
%   kd, kv, kp := scalar positive gains
%   magW := positive scalar representing maximum velocity magnitude
%   beta := positive scalar representing minimum convergence rate
% OUTPUTS:
%   A_all := matrix array of all gershgorin disc constraints

%% Default Parameters
m4 = 1; % fix the scale of the metric on TSO3xSO3
flagSymbolic = false;
%optional parameters
ivarargin=1;
len_varargin = length(varargin);
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'm4'
            ivarargin = ivarargin + 1;
            m4 = varargin{ivarargin};
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 <= len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end 
    end
    ivarargin=ivarargin+1;
end
if contains(class(mag_R),'sym','IgnoreCase',true)
    flagSymbolic = true;
    syms m4 real
end
%% Init matrices
% NOTE: we are computing [1 m1 m2 m3 m5 m6]*A*[[1 m1 m2 m3 m5 m6], but only
% the submatrix A(2:end,2:end) needs to be PSD which reduces to just the
% [2x2] matrix given as [m5 m6]*Bi*[m5;m6]
% Init matrix based on symbolic or not
if flagSymbolic
    A_11 = sym(zeros(6,6)); A_22 = A_11; A_33 = A_11;
else
    A_11 = zeros(6,6); A_22 = A_11; A_33 = A_11;
end
A_11(1,3) = -kd*mag_R/2*cot(mag_R/2);% - 1/4*mag_W^2;
% Add terms from hat3(w)^2 = w*w'-norm(w)^2*eye(3)
%A_11(5,6) = mag_W^2/(4*m4); 
A_11 = (A_11 + A_11')/2;
A1_mat = {A_11};
A_22(1,3) = 1; A_22(1,4) = (beta-kv); A_22 = (A_22 + A_22')/2;
A2_mat = {A_22};
A_33(1,1) = m4*beta-kp*m4*mag_RRef/2*cot(mag_RRef/2);
A_33 = (A_33 + A_33')/2;
A3_mat = {A_33};

% Define possible terms resulting from absolute values
% For disc centered at D(1,1)
tempMat1 = create_empty_mat(flagSymbolic,6);
tempMat2 = create_empty_mat(flagSymbolic,6);
% Fill entries
tempMat1(1,2) = beta;
tempMat1 = (tempMat1 + tempMat1')/2;
tempMat2(1,2) = sqrt(3)*beta;
tempMat2(1,3) = sqrt(3)*(-1/4*mag_W^2);
tempMat2(5,6) = sqrt(3)*(mag_W^2/(4*m4));
tempMat2 = (tempMat2 + tempMat2')/2;
D1_abs_fun{1} = {tempMat1,tempMat2,-tempMat2};
% D1_abs_fun{1} = {sqrt(3)*beta*m1,...
%     -sqrt(3)*beta*m1,...
%     sqrt(3)*(-1/4*(m2-m5*m6/m4)*mag_W^2+beta*m1),...
%     -sqrt(3)*(-1/4*(m2-m5*m6/m4)*mag_W^2+beta*m1)};
% 
tempMat1 = create_empty_mat(flagSymbolic,6);
tempMat2 = create_empty_mat(flagSymbolic,6);
tempMat3 = create_empty_mat(flagSymbolic,6);
% Fill entries
tempMat1(1,2) = sqrt(3)*1/2;
tempMat1(1,3) = sqrt(3)*(-1/2*kv + beta);
tempMat1(1,4) = -sqrt(3)*kd/2;
tempMat1 = (tempMat1 + tempMat1')/2;
tempMat2(1,4) = sqrt(3)*(-kd/2*mag_R/2*cot(mag_R/2));
tempMat2(1,2) = sqrt(3)*1/2;
tempMat2(1,3) = sqrt(3)*(-1/2*kv + beta);
tempMat2(5,5) = sqrt(3)*kd/(4*m4)*mag_R;
tempMat2 = (tempMat2 + tempMat2')/2;
tempMat3(1,4) = -sqrt(3)*(-kd/2*mag_R/2*cot(mag_R/2));
tempMat3(1,2) = -sqrt(3)*1/2;
tempMat3(1,3) = -sqrt(3)*(-1/2*kv + beta);
tempMat3(5,5) = sqrt(3)*kd/(4*m4)*mag_R;
tempMat3 = (tempMat3 + tempMat3')/2;
D1_abs_fun{2} = {tempMat1, -tempMat1, tempMat2, tempMat3};
% D1_abs_fun{2} = {sqrt(3)*(-kd/2*m3+1/2*(m1-m2*kv)+beta*m2),...
%     -sqrt(3)*(-kd/2*m3+1/2*(m1-m2*kv)+beta*m2),...
%     sqrt(3)*(-kd/2*m3*mag_R/2*cot(mag_R/2)+1/2*(m1-m2*kv)+beta*m2) + sqrt(3)*kd*m5^2/(4*m4)*mag_R,... % NOTE: abs(a*1i) = a for a>=0
%     -sqrt(3)*(-kd/2*m3*mag_R/2*cot(mag_R/2)+1/2*(m1-m2*kv)+beta*m2) + sqrt(3)*kd*m5^2/(4*m4)*mag_R};
%
tempMat1 = create_empty_mat(flagSymbolic,6);
tempMat2 = create_empty_mat(flagSymbolic,6);
tempMat1(1,4) = sqrt(3)*(-1/8*mag_W^2);
tempMat1(5,5) = sqrt(3)*(mag_W^2/(8*m4));
tempMat1 = (tempMat1 + tempMat1')/2;
tempMat2(1,3) = sqrt(3)*(-mag_W/4);
tempMat2(5,6) = sqrt(3)*(mag_W/(4*m4));
tempMat2(1,4) = sqrt(3)*(kv*mag_W/4);
tempMat2(5,5) = sqrt(3)*(-mag_W*kv/(4*m4));
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
tempMat1(1,3) = sqrt(3)*(kd/2);
tempMat1(1,5) = sqrt(3)*(-kd/2);
tempMat1 = (tempMat1 + tempMat1')/2;
tempMat2(1,3) = sqrt(3)*(kd/2*mag_R/2*cot(mag_R/2));
tempMat2(1,5) = sqrt(3)*(-kd/2*mag_R/2*cot(mag_R/2));
tempMat2 = (tempMat2 + tempMat2')/2;
tempMat3(1,3) = sqrt(3)*(kd/2*mag_R/2);
tempMat3(1,5) = sqrt(3)*(-kd/2*mag_R/2);
tempMat3(5,6) = sqrt(3)*(-kd*mag_R/(4*m4));
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
tempMat1(1,6) = sqrt(3)*(-kp/2+beta); tempMat1 = (tempMat1 + tempMat1')/2;
D1_abs_fun{5} = {tempMat1, -tempMat1};
% D1_abs_fun{5} = {sqrt(3)*(-m6*kp/2+beta*m6),...
%     -sqrt(3)*(-m6*kp/2+beta*m6)};
%
tempMat1 = create_empty_mat(flagSymbolic,6);
tempMat1(6,6) = sqrt(3)*mag_W/(4*m4);
tempMat1(5,6) = sqrt(3)*(-mag_W*kv/(4*m4)); tempMat1 = (tempMat1+tempMat1')/2;
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
tempMat1(1,4) = sqrt(3)*kd/2;
tempMat1(1,6) = sqrt(3)*1/2;
tempMat1(1,5) = sqrt(3)*(-kv/2); tempMat1 = (tempMat1 + tempMat1')/2;
tempMat2(1,4) = sqrt(3)*kd/2*mag_R/2*cot(mag_R/2);
tempMat2(1,6) = sqrt(3)*1/2;
tempMat2(1,5) = sqrt(3)*(-kv/2); tempMat2 = (tempMat2 + tempMat2')/2;
tempMat3(1,4) = sqrt(3)*(kd/2*mag_R/2);
tempMat3(5,5) = sqrt(3)*(-kd/(4*m4)*mag_R); tempMat3 = (tempMat3 + tempMat3')/2;
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
tempMat1(1,5) = (-kp/2 + beta); tempMat1 = (tempMat1 + tempMat1')/2;
D2_abs_fun{4} = {tempMat1, -tempMat1};
% D2_abs_fun{4} = {-m5*kp/2+beta*m5,...
%     -(-m5*kp/2+beta*m5)};
%
tempMat1 = create_empty_mat(flagSymbolic,6);
tempMat1(5,6) = sqrt(3)*mag_W/(4*m4);
tempMat1(5,5) = sqrt(3)*(-kv*mag_W/(4*m4)); tempMat1 = (tempMat1 + tempMat1')/2;
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
tempMat1(1,5) = sqrt(3)*kd; tempMat1 = (tempMat1 + tempMat1')/2;
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
end