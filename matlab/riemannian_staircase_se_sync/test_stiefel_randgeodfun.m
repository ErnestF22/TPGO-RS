close all;
clear;
clc;


d = 3;
num_rows_stiefel = 4;
N = 5;

A_stiefel = make_rand_stiefel_mat(num_rows_stiefel, d, N);

disp("A_stiefel");
disp(A_stiefel);


% stiefel_randn(d*num_rows_stiefel,N)
% st_rn = stiefel_randn(eye(d*num_rows_stiefel),N);
[stief_t,vt,stief_0,v0,vVec,~,~] = stiefel_randGeodFun(matStack(A_stiefel));

disp("stief_t(0)");
disp(stief_t(0));

ids = reshape(1:num_rows_stiefel*N, num_rows_stiefel, []);

eps = 1e-5;
for ii = ids
    stief_0_ii = stief_0(ii, :);
    if max(abs(eye(d) - stief_0_ii'*stief_0_ii), [], "all") > eps
        fprintf("%g mat is not on Stiefel\n", ii);
    end
end

% funPlot(@(t) stief_t(t)'*stief_t(t))
t=1;
stief_t(t)'*stief_t(t)
