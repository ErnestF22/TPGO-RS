function essential_toRT_test
R=rot_randn();
T=cnormalize(randn(3,1));

Q=essential_fromRT(R,T);
[RQ,TQ]=essential_toRT(Q);

disp([RT2G(R,T) RT2G(RQ,TQ)])
