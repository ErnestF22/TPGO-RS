function [Yt,Ht]=stiefel_rand_geod_fun_new(Y0)
v0=stiefel_randTangentNormVector(Y0);
Yt=@(t) stiefel_exp(Y0,t*v0);
Ht=@(t) tangent(Y0,v0,t);

function Ht=tangent(Y0,v0,t)
[~,Ht]=stiefel_exp(Y0,t*v0);