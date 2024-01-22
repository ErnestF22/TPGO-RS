function c=rsom_cost_rot_stiefel(x,problem)
xCost=matStack(multitransp(x));
c=trace(xCost*problem.P) + problem.frct;
end