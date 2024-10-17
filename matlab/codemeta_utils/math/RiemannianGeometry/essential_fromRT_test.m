function essential_fromRT_test
R1=rot_randn();
R2=rot_randn();
T1=randn(3,1);
T2=randn(3,1);

E1=R1'*hat(cnormalize(T2-T1))*R2;

E2=essential_toE(essential_fromRT(R1,T1,R2,T2));

disp([E1 E2])

