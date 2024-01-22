function g=rsom_rgrad_rot_stiefel(x,problem)
g=stiefel_tangentProj(x,som_egrad_rot_stiefel(x,problem));
end
