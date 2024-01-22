function g=rsom_egrad_rot_stiefel(x,problem)
xStack=matStack(x);
g=matUnstack(problem.P,problem.sz(1));
end