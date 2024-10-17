function POCMultiviewMatrix
resetRands()

X=randn(3,1);
R1=rot_randn();
R2=rot_randn();
R3=rot_randn();
T1=randn(3,1);
T2=randn(3,1);
T3=randn(3,1);

x1=hom(projectFromRT(R1,T1,X));
x2=hom(projectFromRT(R2,T2,X));
x3=hom(projectFromRT(R3,T3,X));


z=zeros(3,1);
I=eye(3);
Z=zeros(3);

M=[R1 T1 x1 z z; R2 T2 z x2 z; R3 T3 z z x3];

size(M)

%rank(M)
S1=blkdiag(R1',I,I);
M1=S1*M;

display(M1)

S2=[I -R1'*T1 Z;z' 1 z'; Z z I];
M2=M1*S2;
display(M2)

S3=[I z -R1'*x1 zeros(3,2); z' 1 z'; z' 0 1 zeros(1,2); zeros(2,5) eye(2)]; 
M3=M2*S3;
display(M3)

disp(rank(M)-rank(M3))

N=M3(4:end,4:end);

display(N)

D=[x2' z'; hat(x2) Z; z' x3'; Z hat(x3)];
N1=D*N;

display(N1)

Mp=N1([2:4 6:8],1:2);

display(Mp)

rank(Mp)

a1=Mp(1:3,1);
b1=Mp(1:3,2);
a2=Mp(4:6,1);
b2=Mp(4:6,2);

a1*b2'-b1*a2'



function x=hom(x)
x=[x;1];