function g=som_egrad_rot_stiefel(x,problem)
xStack=matStack(x);
g=matUnstack((problem.L+problem.L')*xStack+problem.P,problem.sz(1));
end