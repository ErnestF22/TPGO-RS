function g=rsom_egrad_transl_stiefel(x,problem)
% x = reshape(x, problem.sz(1), problem.sz(3));
g=x*(problem.LR+problem.LR')+(problem.PR)';
end