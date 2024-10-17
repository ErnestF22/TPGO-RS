function lie_funs=stiefel_funs()
lie_funs.dim=@stiefel_dim;

lie_funs.exp=@stiefel_exp;
lie_funs.log=@stiefel_log;
lie_funs.dist=@stiefel_dist;
lie_funs.mean=@stiefel_mean;
lie_funs.metric=@stiefel_metric;
lie_funs.eye=@stiefel_eye;

lie_funs.tangentProj=@stiefel_tangentProj;
lie_funs.tangentBasis=@stiefel_tangentBasis;
lie_funs.randTangentNormVector=@stiefel_randTangentNormVector;

lie_funs.eucl2RiemGrad=@stiefel_eucl2RiemGrad;

lie_funs.retractions=@stiefel_retractions;

