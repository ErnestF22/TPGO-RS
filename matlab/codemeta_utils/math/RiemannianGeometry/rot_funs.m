function lie_funs=rot_funs()
lie_funs.dim=@rot_dim;

lie_funs.exp=@rot_exp;
lie_funs.log=@rot_log;

lie_funs.dist=@rot_distn;
lie_funs.mean=@rot_mean;
lie_funs.metric=@rot_metric;

lie_funs.eye=@rot_eye;

lie_funs.parallel=@rot_parallel;

lie_funs.tangentBasis=@rot_tangentBasis;
lie_funs.tangentProj=@rot_tangentProj;
lie_funs.eucl2RiemGrad=@rot_eucl2RiemGrad;
lie_funs.vee=@rot_vee;
lie_funs.hat=@rot_hat;

lie_funs.randTangentNormVector=@rot_randTangentNormVector;

lie_funs.anchors=@rot_anchors;

lie_funs.retractions=@rot_retractions;
