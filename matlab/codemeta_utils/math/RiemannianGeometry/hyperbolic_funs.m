function lf=hyperbolic_funs()
lf.exp=@hyperbolic_exp;
lf.log=@hyperbolic_log;
lf.metric=@hyperbolic_metric;

lf.dist=@hyperbolic_dist;

lf.tangentProj=@hyperbolic_tangentProj;
lf.randTangentNormVector=@hyperbolic_randTangentNormVector;

lf.mean=@hyperbolic_mean;
