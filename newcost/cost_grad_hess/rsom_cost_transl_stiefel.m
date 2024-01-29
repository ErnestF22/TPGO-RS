function c=rsom_cost_transl_stiefel(x,problem)
% x = reshape(x, problem.sz(1), problem.sz(3));
LR = problem.LR;
PR = problem.PR;
BR = problem.BR;
c = trace(x * LR * x') + trace(x * PR) + trace(BR);
end