function c=som_cost_rot_stiefel(x,problem)
xStack=matStack(x);
c=trace(xStack'*problem.L*xStack+xStack'*problem.P) + problem.fixed_cost_term;
end