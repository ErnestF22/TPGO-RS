function eig_bool = eigencheck_hessian(lambda, v, hess_fun_han, thr)
% EIGENCHECK_HESSIAN Check if lambda is an Eigenvalue for
% Hessian defined by hess_fun_han

if ~exist('thr','var')
  thr = 1e-3;
end

hess_v = hess_fun_han(v);

diff = norm((lambda)*v(:) - hess_v(:),'inf');


%The infinity norm is defined as the absolute value of the largest 
% component of the vector
if diff < thr
    fprintf("%g is an eigenvalue\n", lambda)
    eig_bool = boolean(1);
else
    fprintf("%g is NOT an eigenvalue: diff %g > thr\n", lambda, diff);
    eig_bool = boolean(0);
end

end %end file function

