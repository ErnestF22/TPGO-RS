function h=rsom_ehess_rot_stiefel(x,u,problem)
h = matUnstack((problem.L + problem.L')*matStack(u), problem.sz(1));
end