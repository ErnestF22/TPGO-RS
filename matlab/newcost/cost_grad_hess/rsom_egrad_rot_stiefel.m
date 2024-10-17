function g=rsom_egrad_rot_stiefel(x,problem)
% xStack=matStack(x);
d = size(x, 2);
g=matUnstackH(problem.P,d);
end