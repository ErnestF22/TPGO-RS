% POC_rotRef_allBoundEqns.m

%% Setup
close all;clear all; clc;

syms kd kv kp beta m1 m2 m3 m4 m5 m6 mag_R mag_RRef mag_W real;

% init each set of functions with the centriod with no absolute values
% terms
D1_Bounds = {-m2*kd*mag_R/2*cot(mag_R/2)};
D2_Bounds = {m2-m3*kv+beta*m3};
D3_Bounds = {-m4*kp*mag_RRef/2*cot(mag_RRef/2) + beta*m4};

% Define the absolute value terms as their pos/neg parts. Each absolute
% value term will be one cell with a vector of possible values
% For disc centered at D(1,1)
D1_abs_fun{1} = {sqrt(3)*beta*m1,...
    -sqrt(3)*beta*m1,...
    sqrt(3)*(-1/4*(m2-m5*m6/m4)*mag_W^2+beta*m1),...
    -sqrt(3)*(-1/4*(m2-m5*m6/m4)*mag_W^2+beta*m1)};
D1_abs_fun{2} = {sqrt(3)*(-kd/2*m3+1/2*(m1-m2*kv)+beta*m2),...
    -sqrt(3)*(-kd/2*m3+1/2*(m1-m2*kv)+beta*m2),...
    sqrt(3)*(-kd/2*m3*mag_R/2*cot(mag_R/2)+1/2*(m1-m2*kv)+beta*m2) + sqrt(3)*kd*m5^2/(4*m4)*mag_R,... % NOTE: abs(a*1i) = a for a>=0
    -sqrt(3)*(-kd/2*m3*mag_R/2*cot(mag_R/2)+1/2*(m1-m2*kv)+beta*m2) + sqrt(3)*kd*m5^2/(4*m4)*mag_R};
D1_abs_fun{3} = {sqrt(3)*(-1/8*(m3-m5^2/m4)*mag_W^2)+sqrt(3)*(-1/4*(m2-m5*m6/m4-m3*kv+m5^2/m4*kv)*mag_W),...
    -sqrt(3)*(-1/8*(m3-m5^2/m4)*mag_W^2)+sqrt(3)*(-1/4*(m2-m5*m6/m4-m3*kv+m5^2/m4*kv)*mag_W),...
    sqrt(3)*(-1/8*(m3-m5^2/m4)*mag_W^2)-sqrt(3)*(-1/4*(m2-m5*m6/m4-m3*kv+m5^2/m4*kv)*mag_W),...
    -sqrt(3)*(-1/8*(m3-m5^2/m4)*mag_W^2)-sqrt(3)*(-1/4*(m2-m5*m6/m4-m3*kv+m5^2/m4*kv)*mag_W)};
D1_abs_fun{4} = {sqrt(3)*(kd/2*(m2-m5)),...
    -sqrt(3)*(kd/2*(m2-m5)),...
    sqrt(3)*(kd/2*(m2-m5)*mag_R/2*cot(mag_R/2)) + sqrt(3)*(kd/2*(m2-m5)*mag_R/2 - m5*m6/(4*m4)*kd*mag_R),...
    sqrt(3)*(kd/2*(m2-m5)*mag_R/2*cot(mag_R/2)) - sqrt(3)*(kd/2*(m2-m5)*mag_R/2 - m5*m6/(4*m4)*kd*mag_R),...
    -sqrt(3)*(kd/2*(m2-m5)*mag_R/2*cot(mag_R/2)) + sqrt(3)*(kd/2*(m2-m5)*mag_R/2 - m5*m6/(4*m4)*kd*mag_R),...
    -sqrt(3)*(kd/2*(m2-m5)*mag_R/2*cot(mag_R/2)) - sqrt(3)*(kd/2*(m2-m5)*mag_R/2 - m5*m6/(4*m4)*kd*mag_R)};
D1_abs_fun{5} = {sqrt(3)*(-m6*kp/2+beta*m6),...
    -sqrt(3)*(-m6*kp/2+beta*m6)};
D1_abs_fun{6} = {sqrt(3)*((m6^2-m5*m6*kv)/(4*m4)*mag_W),...
    -sqrt(3)*((m6^2-m5*m6*kv)/(4*m4)*mag_W)};
% For disc centered at D(2,2)
D2_abs_fun{1} = D1_abs_fun{2};
D2_abs_fun{2} = D1_abs_fun{3};
D2_abs_fun{3} = {sqrt(3)*(m3*kd/2+(m6-m5*kv)/2),...
    -sqrt(3)*(m3*kd/2+(m6-m5*kv)/2),...
    sqrt(3)*(m3*kd/2*mag_R/2*cot(mag_R/2)+(m6-m5*kv)/2) + sqrt(3)*(m3*kd/2*mag_R/2-m5^2*kd/(4*m4)*mag_R),...
    sqrt(3)*(m3*kd/2*mag_R/2*cot(mag_R/2)+(m6-m5*kv)/2) - sqrt(3)*(m3*kd/2*mag_R/2-m5^2*kd/(4*m4)*mag_R),...
    -sqrt(3)*(m3*kd/2*mag_R/2*cot(mag_R/2)+(m6-m5*kv)/2) + sqrt(3)*(m3*kd/2*mag_R/2-m5^2*kd/(4*m4)*mag_R),...
    -sqrt(3)*(m3*kd/2*mag_R/2*cot(mag_R/2)+(m6-m5*kv)/2) - sqrt(3)*(m3*kd/2*mag_R/2-m5^2*kd/(4*m4)*mag_R)};
D2_abs_fun{4} = {-m5*kp/2+beta*m5,...
    -(-m5*kp/2+beta*m5)};
D2_abs_fun{5} = {sqrt(3)*( (m5*m6-m5^2*kv)/(4*m4)*mag_W ),...
    -sqrt(3)*( (m5*m6-m5^2*kv)/(4*m4)*mag_W )};
% For disc centered at D(3,3)
D3_abs_fun{1} = D1_abs_fun{4};
D3_abs_fun{2} = D1_abs_fun{5};
D3_abs_fun{3} = D1_abs_fun{6};
D3_abs_fun{4} = D2_abs_fun{3};
D3_abs_fun{5} = D2_abs_fun{4};
D3_abs_fun{6} = D2_abs_fun{5};
D3_abs_fun{7} = {sqrt(3)*m5*kd,...
    -sqrt(3)*m5*kd};
    
%% Construct all possible functions
% For disc centered at D(1,1)
for ii = 1:length(D1_abs_fun)
    A = D1_abs_fun{ii}; % this should now be a vector of cells of functions
    len_A = length(A);
    % expand the D1_Bounds by length of A
    % NOTE: repelem repeats elements consecutivly
    D1_Bounds = repelem(D1_Bounds,len_A);
    % Add each element of A to each element of D1_Bounds
    for jj = 1:length(D1_Bounds)
        D1_Bounds{jj} = D1_Bounds{jj} + A{mod(jj,len_A)+1};
    end
end
% For disc centered at D(2,2)
for ii = 1:length(D2_abs_fun)
    A = D2_abs_fun{ii}; % this should now be a vector of cells of functions
    len_A = length(A);
    % expand the D1_Bounds by length of A
    % NOTE: repelem repeats elements consecutivly
    D2_Bounds = repelem(D2_Bounds,len_A);
    % Add each element of A to each element of D1_Bounds
    for jj = 1:length(D2_Bounds)
        D2_Bounds{jj} = D2_Bounds{jj} + A{mod(jj,len_A)+1};
    end
end
% For disc centered at D(3,3)
for ii = 1:length(D3_abs_fun)
    A = D3_abs_fun{ii}; % this should now be a vector of cells of functions
    len_A = length(A);
    % expand the D1_Bounds by length of A
    % NOTE: repelem repeats elements consecutivly
    D3_Bounds = repelem(D3_Bounds,len_A);
    % Add each element of A to each element of D1_Bounds
    for jj = 1:length(D3_Bounds)
        D3_Bounds{jj} = D3_Bounds{jj} + A{mod(jj,len_A)+1};
    end
end

