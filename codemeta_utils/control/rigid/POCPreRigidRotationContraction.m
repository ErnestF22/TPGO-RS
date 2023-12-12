function POCPreRigidRotationContraction
global Av DAv k g Dg c m11 m12 m22

resetRands()
x0=zeros(6,1);
[xt,dxt]=real_randGeodFun(x0);
vt=@(t) stateToVel(xt(t));
dvt=@(t) stateToVel(dxt(t));

k=1;
c=2;
J=diag([1 2 3]);
H=rand(3);
g=@(p) H*(p-stateToPos(x0));
Dg=@(p) H;
Av=@(v) hat3(J*v);
DAv=@(v) -hat3(v)*J+hat3(J*v);
%funCheckDer(@(t) Av(vt(t))*vt(t), @(t) DAv(vt(t))*dvt(t))


m11=1;
m12=2;
m22=4;
M=kron([m11 m12; m12 m22],eye(3));
ft=@(t) field(xt(t));
Dft=@(t) Dfield(xt(t));
%funCheckDer(ft, @(t) Dft(t)*dxt(t))

MDft=@(t) MDfield(xt(t));
%funCompare(MDft,@(t) M*Dft(t))

global symm
symm=@(A) (A+A')/2;
asym=@(A) (A-A')/2;
global Sg Ag SDAv ADAv
Sg=@(p) symm(Dg(p));
Ag=@(p) asym(Dg(p));
SDAv=@(v) symm(DAv(v));
ADAv=@(v) asym(DAv(v));

SMDft=@(t) symmMDfield(xt(t));
%funCompare(SMDft,@(t) symm(MDft(t)))

global b
b=0.5;
SMDfbetat=@(t) symmMDfieldBeta(xt(t));
%funCompare(SMDfbetat, @(t) SMDft(t)+b*M)

%SMDfDcompt=@(t) symmMDfieldBetaDecomp(xt(t));
%funCompare(SMDfDcompt,SMDfbetat)

SMDfDcompt=@(t) symmMDfieldBetaDecomp1(xt(t)) ...
    +symmMDfieldBetaDecomp2(xt(t)) ...
    +symmMDfieldBetaDecomp3(xt(t)) ...
    +symmMDfieldBetaDecomp4(xt(t));
funCompare(SMDfDcompt,SMDfbetat)

SMDfDcomp1t=@(t) symmMDfieldBetaDecomp1(xt(t));
funCompare(@(t) symmMDfieldBetaDecomp1U(xt(t)),SMDfDcomp1t)

function r=stateToPos(x)
r=x(1:3,:);

function v=stateToVel(x)
v=x(4:6,:);

function U=USg(p)
global Sg
[U,~]=eig(Sg(p));

function L=LSg(p)
global Sg
U=USg(p);
L=diag(diag(U'*Sg(p)*U));

function f=field(x)
global k g c Av
p=stateToPos(x);
v=stateToVel(x);
f=[v; -k*g(p)-c*v+Av(v)*v];

function Df=Dfield(x)
global k Dg c DAv
p=stateToPos(x);
v=stateToVel(x);
Df=[zeros(3) eye(3); -k*Dg(p) -c*eye(3)+DAv(v)];

function MDf=MDfield(x)
global m11 m12 m22
global k Dg c DAv
p=stateToPos(x);
v=stateToVel(x);
I=eye(3);
MDf=[-k*m12*Dg(p) (m11-c*m12)*I+m12*DAv(v);
     -k*m22*Dg(p) (m12-c*m22)*I+m22*DAv(v)];

function SMDf=symmMDfield(x)
global m11 m12 m22
global k Dg Sg c DAv SDAv
p=stateToPos(x);
v=stateToVel(x);
I=eye(3);
C=(-k*m22*Dg(p)'+(m11-c*m12)*I+m12*DAv(v))/2;
SMDf=[-k*m12*Sg(p) C;
      C' (m12-c*m22)*I+m22*SDAv(v)];

function SMDf=symmMDfieldBeta(x)
global m11 m12 m22
global k Dg Sg c DAv SDAv
global b
p=stateToPos(x);
v=stateToVel(x);
I=eye(3);
A=-k*m12*Sg(p)+m11*b*I;
C=(-k*m22*Dg(p)'+(m11-c*m12)*I+m12*DAv(v))/2+m12*b*I;
D=(m12-c*m22+b*m22)*I+m22*SDAv(v);
SMDf=[A C;
      C' D];

function SMDf=symmMDfieldBetaDecomp(x)
global m11 m12 m22
global k Ag Sg c SDAv ADAv
global b
p=stateToPos(x);
v=stateToVel(x);
Z=zeros(3);
I=eye(3);
A=-k*m12*Sg(p)+m11*b*I;
C1=(-k*m22*Sg(p)'+(m11-c*m12)*I)/2+m12*b*I;
D1=(m12-c*m22+b*m22)*I;
C2=-k*m22*Ag(p)'/2;
C3=m12*SDAv(v)/2;
D3=m22*SDAv(v);
C4=m12*ADAv(v)/2;
SMDf=[A C1; C1' D1]...
    +[Z C2; C2' Z]...
    +[Z C3; C3' D3]...
    +[Z C4; C4' Z];

function SMDf=symmMDfieldBetaDecomp1(x)
global m11 m12 m22
global k Sg c
global b
p=stateToPos(x);
I=eye(3);
A=-k*m12*Sg(p)+m11*b*I;
C1=(-k*m22*Sg(p)'+(m11-c*m12)*I)/2;
D1=(m12-c*m22+b*m22)*I;
SMDf=[A C1+m12*b*I; C1'+m12*b*I D1];

function SMDf=symmMDfieldBetaDecomp2(x)
global m22
global k Ag 
p=stateToPos(x);
Z=zeros(3);
C2=-k*m22*Ag(p)'/2;
SMDf=[Z C2; C2' Z];

function SMDf=symmMDfieldBetaDecomp3(x)
global m12 m22
global SDAv
v=stateToVel(x);
Z=zeros(3);
C3=m12*SDAv(v)/2;
D3=m22*SDAv(v);
SMDf=[Z C3; C3 D3];

function SMDf=symmMDfieldBetaDecomp4(x)
global m12
global ADAv
v=stateToVel(x);
Z=zeros(3);
C4=m12*ADAv(v)/2;
SMDf=[Z C4; C4' Z];

function SMDf=symmMDfieldBetaDecomp1U(x)
global m11 m12 m22
global k c
global b
p=stateToPos(x);
Ug=USg(p);
Lg=LSg(p);
lg=diag(Lg);
UgBlock=blkdiag(Ug,Ug);
I=eye(3);
a=zeros(3,1);
b1=zeros(3,1);
c1=zeros(3,1);
d1=zeros(3,1);
for idx=1:3
    a(idx)=-k*m12*lg(idx)+m11*b;
    b1(idx)=(-k*m22*lg(idx)+m11-c*m12)/2+m12*b;
    c1(idx)=(-k*m22*lg(idx)+m11-c*m12)/2+m12*b;
    d1(idx)=m12+(b-c)*m22;
end
A=diag(a);
B1=diag(b1);
C1=diag(c1);
D1=diag(d1);
SMDf=UgBlock*[A B1; C1 D1]*UgBlock';

%A=-k*m12*Lg+m11*b*I;
%C1=(-k*m22*Lg+(m11-c*m12)*I)/2;
%D1=(m12-c*m22+b*m22)*I;
%SMDf=UgBlock*[A C1+m12*b*I; C1+m12*b*I D1]*UgBlock';
