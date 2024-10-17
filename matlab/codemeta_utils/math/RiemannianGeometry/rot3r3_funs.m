function lie_funs=rot3r3_funs()
lie_funs.dim=6;

lie_funs.exp=@rot3r3_exp;
lie_funs.log=@rot3r3_log;

lie_funs.metric=@rot3r3_metric;

lie_funs.hat=@rot3r3_hat;
lie_funs.vee=@rot3r3_vee;
