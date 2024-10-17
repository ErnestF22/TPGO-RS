function lf=essential_signed_funs()
lf.dim=@(x) 5;

lf.exp=@essential_exp;
lf.log=@(Q1,Q2,varargin) essential_log(Q1,Q2,'signed',varargin{:});

lf.dist=@(Q1,Q2,varargin) essential_dist(Q1,Q2,'signed',varargin{:});
lf.metric=@essential_metric;

lf.eye=@essential_eye;

lf.hat=@essential_hat;
lf.vee=@essential_vee;

lf.tangentProj=@essential_tangentProj;

lf.randTangentNormVector=@essential_randTangentNormVector;
