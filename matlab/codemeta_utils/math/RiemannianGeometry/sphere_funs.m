function lie_funs=sphere_funs()
lie_funs.exp=@sphere_exp;
lie_funs.log=@sphere_log;
lie_funs.dist=@sphere_dist;
lie_funs.dim=@(y) size(y,1)-1;
lie_funs.mean=@sphere_mean;
lie_funs.metric=@sphere_metric;

lie_funs.tangentProj=@sphere_tangentProj;
lie_funs.eucl2RiemGrad=@sphere_eucl2RiemGrad;

lie_funs.randTangentNormVector=@sphere_randTangentNormVector;
lie_funs.vee=@sphere_vee;
lie_funs.tangentBasis=@sphere_tangentBasis;

lie_funs.eye=@sphere_eye;
