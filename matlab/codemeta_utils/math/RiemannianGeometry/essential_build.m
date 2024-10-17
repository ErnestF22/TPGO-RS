%Build essential matrix E from a QREM representation
%function Q=essential_build(R1,T1,R2,T2)
function Q=essential_build(R1,T1,R2,T2)
if nargin<3
    R2=R1;
    T2=T1;
    R1=eye(3);
    T2=zeros(3,1);
end

R0=computeAlignmentRotation(T2-T1);
Q(1:3,1:3)=R0*R1;
Q(4:6,1:3)=R0*R2;

%compute a rotation that alings x to e3
function R0=computeAlignmentRotation(x)
v=cross(x,[0;0;1]);
nvSq=v'*v;
if nvSq<1e-15
    R0=eye(3);
else
    v=v/sqrt(nvSq);
    nx=sqrt(x'*x);
    a=acos(min(1,max(-1,x(3)/nx)));
    R0=rot(a*v);
end

    