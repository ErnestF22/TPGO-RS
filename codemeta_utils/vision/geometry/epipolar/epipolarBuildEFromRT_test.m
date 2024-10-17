function epipolarBuildEFromRT_test
R1pose=rot_randn();
R2pose=rot_randn();
T1pose=randn(3,1);
T2pose=randn(3,1);

[R12,T12]=computeRelativePoseFromRT(R1pose,T1pose,R2pose,T2pose,'poses');

[R1ref,T1ref]=invRT(R1pose,T1pose);
[R2ref,T2ref]=invRT(R2pose,T2pose);

E1=epipolarBuildEFromRT(R1ref,T1ref,R2ref,T2ref,'references');
E2=epipolarBuildEFromRT(R1pose,T1pose,R2pose,T2pose,'poses');
E3=epipolarBuildEFromRT(R12,T12);

disp([E1 E1-E1;E2 E2-E1;E3 E3-E1])
