function rot3r3_geodFun_test
G0=RT2G(rot_randn(),randn(3,1));
v=cnormalize(randn(6,1));
[Gt,vGt]=rot3r3_geodFun(G0,v);
check_der(Gt,vGt)
