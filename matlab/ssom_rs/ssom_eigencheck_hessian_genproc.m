function eig_bool = ssom_eigencheck_hessian_genproc(lambda, v, hess_fun_han, thr)
% EIGENCHECK_HESSIAN Check if lambda is an Eigenvalue for
% Hessian defined by hess_fun_han

if ~exist('thr','var')
  thr = 1e-3;
end

hess_v = [hess_fun_han(v).R(:); hess_fun_han(v).T(:); hess_fun_han(v).lambda(:)];
v_full = [v.R(:); v.T(:); v.lambda(:)];
diff = norm((lambda)*v_full(:) - hess_v(:),'inf');

if diff < thr
    fprintf("%g is a GENPROC eigenvalue\n", lambda);
    eig_bool = boolean(1);
else
    fprintf("%g is NOT a GENPROC eigenvalue: diff %g > thr\n", ...
        lambda, diff);
    eig_bool = boolean(0);
end




end %end file function

