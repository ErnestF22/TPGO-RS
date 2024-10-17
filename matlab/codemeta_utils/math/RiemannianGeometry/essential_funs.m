function lf=essential_funs()
lf.dim=@(x) 5;

lf.exp=@essential_exp;
lf.log=@essential_log;

lf.dist=@essential_dist;
lf.metric=@essential_metric;

lf.eye=@essential_eye;

lf.hat=@essential_hat;
lf.vee=@essential_vee;

lf.tangentProj=@essential_tangentProj;

lf.randTangentNormVector=@essential_randTangentNormVector;
