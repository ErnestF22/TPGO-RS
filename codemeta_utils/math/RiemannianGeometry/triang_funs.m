function lie_funs=triang_funs()
lie_funs.exp=@triang_exp;
lie_funs.log=@triang_log;

lie_funs.eye=@triang_eye;

lie_funs.metric=@triang_metric;
lie_funs.tangentBasis=@triang_tangentBasis;
lie_funs.tangentProj=@triang_tangentProject;

lie_funs.randTangentNormVector=@triang_randTangentNormVector;
