function invRT_test
R=rot_randn();
T=randn(3,1);

G=RT2G(R,T);
[invR,invT]=invRT(R,T);

[invg(G) RT2G(invR,invT)]