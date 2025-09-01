function g = ssom_rgrad_R(R, T, lambdas, problem_data)

eg = ssom_egrad_R(R, T, lambdas ,problem_data);

g = stiefel_tangentProj(R, eg);

end %file function