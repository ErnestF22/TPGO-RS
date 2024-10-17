function essential_getRTPair_test
R1=rot_randn();
R2=rot_randn();
T1=cnormalize(randn(3,1));
T2=zeros(3,1);

methodAbsolutePoses='reference';

E=epipolarBuildEFromRT(R1,T1,R2,T2,'methodAbsolutePoses',methodAbsolutePoses);

Q=essential_fromRT(R1,T1,R2,T2,'methodAbsolutePoses',methodAbsolutePoses);
EQ=essential_getE(Q);

[R1p,T1p,R2p,T2p]=essential_getRTPair(Q,'methodAbsolutePoses',methodAbsolutePoses);
Ep=epipolarBuildEFromRT(R1p,T1p,R2p,T2p,'methodAbsolutePoses',methodAbsolutePoses);

[R12,T12]=computeRelativePoseFromRT(R1,T1,R2,T2,'methodAbsolutePoses',methodAbsolutePoses);
[R12p,T12p]=computeRelativePoseFromRT(R1p,T1p,R2p,T2p,'methodAbsolutePoses',methodAbsolutePoses);

disp([R12 T12; R12p T12p])
disp([E;EQ;Ep])
