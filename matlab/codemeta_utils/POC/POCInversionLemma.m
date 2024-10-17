function POCInversionLemma

A=randn(3);
B=randn(3);
C=randn(3);
D=randn(3);
I=eye(3);

M=[A B; C D];
Minv=inv(M);

Ap=Minv(1:3,1:3);
Bp=Minv(1:3,4:6);
Cp=Minv(4:6,1:3);
Dp=Minv(4:6,4:6);

%CpTest=inv(B-A*inv(C)*D);
%ApTest=inv(A)*(I-B*CpTest);

disp([Cp CpTest])
disp([Ap ApTest])

