function h=rsom_rhess_transl_stiefel(x,u,problem)
%Note: on Euclidean space, Eucl. grad = Riem. grad
h = rsom_ehess_transl_stiefel(x,u,problem);
end
