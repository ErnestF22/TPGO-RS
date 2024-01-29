function c=rsom_cost_transl_stiefel(x,problem)
LR = problem.LR;
PR = problem.PR;
BR = problem.BR;
c = trace(x * LR * x') + trace(x * PR) + trace(BR);
end