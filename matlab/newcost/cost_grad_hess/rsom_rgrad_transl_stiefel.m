function g=rsom_rgrad_transl_stiefel(x,problem)
%Note: on Euclidean space, Eucl. grad = Riem. grad
g = rsom_egrad_transl_stiefel(x,problem);
end
