function g=som_rgrad_rot_stiefel(x,problem)
g=stiefel_tangentProj(x,egrad(x,problem));
end
