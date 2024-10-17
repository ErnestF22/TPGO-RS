function h=rsom_ehess_transl_stiefel(x,u,problem)
% x = reshape(x, problem.sz(1), problem.sz(3));
h = u*(problem.LR' + problem.LR);
end