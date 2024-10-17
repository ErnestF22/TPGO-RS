function POCRotationLocalizationTriangle
T1=[0;0;0];
T2=[1;0;0];
T3=[0;1;0];

t23w=T3-T2;
t32w=-t23w;
t21w=T1-T2;
t12w=-t21w;
t31w=T1-T3;
t13w=-t31w;

R10=rot_randn();
R20=rot_randn();
R30=rot_randn();

t23=R20*t23w;
t32=R30*t32w;
t21=R20*t21w;
t12=R10*t12w;
t31=R30*t31w;
t13=R10*t13w;



% [R2,~,~,~,v2]=rot_randGeodFun(R20);
% [R3,~,~,~,v3]=rot_randGeodFun(R30);
v2=t21w;
v3=t31w;
R2=rot_geodFun(R20,rot_hat(R20,v2));
R3=rot_geodFun(R30,rot_hat(R30,v3));

disp('Constraints')
disp(constraints(R10,R20,R30,t12,t21,t13,t31,t23,t32))
% figure(1)
% check_der(@(t) funAndDer(R10,R2(t),R3(t),t12,t21,t13,t31,t23,t32,v2,v3))
disp('Check of non-trivial nullspace of constraints')
disp([nullDiffConstraints(R10,R20,R30,t12,t21,t13,t31,t23,t32)...
    nullDiffConstraints(R10,R2(1),R3(1),t12,t21,t13,t31,t23,t32)])
% figure(2)
% plotfun(@(t) nullDiffConstraints(R10,R2(t),R3(t),t12,t21,t13,t31,t23,t32))


function c=nullDiffConstraints(R1,R2,R3,t12,t21,t13,t31,t23,t32)
[~,De]=constraints(R1,R2,R3,t12,t21,t13,t31,t23,t32);
z=zeros(3,1);
t21p=R2'*t21;
t31p=R3'*t31;
n=[z;t21p;t31p];
c=De*n;

function [c,dc]=funAndDer(R1,R2,R3,t12,t21,t13,t31,t23,t32,v2,v3)
[c,Dc]=constraints(R1,R2,R3,t12,t21,t13,t31,t23,t32);
dc=Dc(:,4:9)*[v2;v3];

function [e,De]=constraints(R1,R2,R3,t12,t21,t13,t31,t23,t32)
e=zeros(9,1);
t12p=R1'*t12;
t21p=R2'*t21;
e(1:3)=t12p+t21p;
t13p=R1'*t13;
t31p=R3'*t31;
e(4:6)=t13p+t31p;
t32p=R3'*t32;
t23p=R2'*t23;
e(7:9)=t32p+t23p;
if nargout>1
    Z=zeros(3);
    De=[hat(t12p) hat(t21p) Z;
        hat(t13p) Z hat(t31p);
        Z hat(t23p) hat(t32p)];
end

