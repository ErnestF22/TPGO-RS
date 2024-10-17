function essential_build_test
R1=eye(3);
T1=zeros(3,1);
R2=rot_randn();
T2=cnormalize(rand(3,1));

Q1=essential_build(R1,T1,R2,T2);
E1=R1'*hat(T2-T1)*R2;
E1exp=essential_getE(Q1);
disp([E1 E1exp])
Q2=essential_build(R2,T2);
E2=hat(T2)*R2;
E2exp=-essential_getE(Q2);
disp([E2 E2exp])
