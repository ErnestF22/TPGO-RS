function eig_bool = eigencheck_hessian_genproc(lambda, v, hess_fun_han, thr)
% EIGENCHECK_HESSIAN Check if lambda is an Eigenvalue for
% Hessian defined by hess_fun_han

if ~exist('thr','var')
  thr = 1e-3;
end

hess_v = hess_fun_han(v);

%% R
v_R = v.R;
lambda_R = lambda.R;
hess_v_R = hess_v.R;

diff_R = norm((lambda_R)*v_R(:) - hess_v_R(:),'inf');


%The infinity norm is defined as the absolute value of the largest 
% component of the vector
if diff_R < thr
    fprintf("R -> %g is an eigenvalue\n", lambda_R);
    eig_bool.R = boolean(1);
else
    fprintf("R -> %g is NOT an eigenvalue: diff %g > thr\n", lambda_R, diff_R);
    eig_bool.R = boolean(0);
end

%% T
v_T = v.T;
lambda_T = lambda.T;
hess_v_T = hess_v.T;

diff_T = norm((lambda_T)*v_T(:) - hess_v_T(:),'inf');


%The infinity norm is defined as the absolute value of the largest 
% component of the vector
if diff_T < thr
    fprintf("T -> %g is an eigenvalue\n", lambda_T)
    eig_bool.T = boolean(1);
else
    fprintf("T -> %g is NOT an eigenvalue: diff %g > thr\n", lambda_T, diff_T);
    eig_bool = boolean(0);
end


end %end file function

