function POCderivativeThetaZForDistance
%resetRands(1)

e3=[0;0;1];
%R0=eye(3);
% [Ra1,~,~,~,vVeca1]=rot_randGeodFun([],'randspeed');
% [Ra2,~,~,~,vVeca2]=rot_randGeodFun([],'randspeed');
% [Rb1,~,~,~,vVecb1]=rot_randGeodFun([],'randspeed');
% [Rb2,~,~,~,vVecb2]=rot_randGeodFun([],'randspeed');

[Q1,~,~,~,v1Vec]=essential_randGeodFun();
[Q2,~,~,~,v2Vec]=essential_randGeodFun();

vVeca1=v1Vec(1:3);
vVeca2=v1Vec(4:6);
vVecb1=v2Vec(1:3);
vVecb2=v2Vec(4:6);

Ra1=@(s) essential_getR1(Q1(s));
Ra2=@(s) essential_getR2(Q1(s));
switch 2
    case 1
        Rb1=@(s) essential_getR1(Q2(s));
        Rb2=@(s) essential_getR2(Q2(s));
        tOpt=@(s) essential_distMinAnglePair([Ra1(s);Ra2(s)],[Rb1(s); Rb2(s)],false);
    case 2
        Q2r=@(s) essential_closestRepresentative(Q1(s),Q2(s));
        Rb1=@(s) essential_getR1(Q2r(s));
        Rb2=@(s) essential_getR2(Q2r(s));
        tOpt=@(s) 0;
end

Rzt=@(t) rot(t*e3);

Razb1=@(t,s) Ra1(s)'*Rzt(t)*Rb1(s);
Razb2=@(t,s) Ra2(s)'*Rzt(t)*Rb2(s);

Rb1r=@(t,s) Rzt(t)*Rb1(s);
Rb2r=@(t,s) Rzt(t)*Rb2(s);

Log1=@(t,s) logrot(Razb1(t,s));
Log2=@(t,s) logrot(Razb2(t,s));

v1=@(t,s) Rb1(s)'*e3;
v2=@(t,s) Rb2(s)'*e3;

f1=@(t,s) 0.5*rot_dist(Ra1(s),Rb1r(t,s))^2;
f2=@(t,s) 0.5*rot_dist(Ra2(s),Rb2r)^2;

df1=@(t,s) Log1(t,s)'*v1(t,s);
df2=@(t,s) Log2(t,s)'*v2(t,s);

f=@(t,s) f1(t,s)+f2(t,s);
df=@(t,s) df1(t,s)+df2(t,s);

t0=rand;
t1=rand;
t2=rand;
t=@(z) t0+z*t1+z^2*t2/2;
dt=@(z) t1+z*t2;

%plotfun(@(s) df(tOpt(s),s),'angle')
s=rand;
%check_der(@(z) f(t(z),s),@(z) df(t(z),s)*dt,'angle')

DLog1=@(t,s) rot3_logDiff(eye(3),Ra1(s)'*Rb1r(t,s));
dsLog1=@(t,s) -DLog1(t,s)*(Razb1(t,s)'*vVeca1-vVecb1);
DLog2=@(t,s) rot3_logDiff(eye(3),Ra2(s)'*Rb2r(t,s));
dsLog2=@(t,s) -DLog2(t,s)*(Razb2(t,s)'*vVeca2-vVecb2);

% t=rand;
% check_der(@(s) Log1(t,s),@(s) dLog1(t,s),'angle')
% check_der(@(s) Log2(t,s),@(s) dLog2(t,s),'angle')

dsv1=@(t,s) -hat(v1(t,s))'*vVecb1;
dsv2=@(t,s) hat(vVecb2)'*v2(t,s);

% check_der(@(s) v1(t,s), @(s) dsv1(t,s),'angle')
% check_der(@(s) v2(t,s), @(s) dsv2(t,s),'angle')

dsdf1=@(t,s) dsLog1(t,s)'*v1(t,s)+Log1(t,s)'*dsv1(t,s);
dsdf2=@(t,s) dsLog2(t,s)'*v2(t,s)+Log2(t,s)'*dsv2(t,s);

% check_der(@(s) df1(t,s), @(s) dsdf1(t,s), 'angle')
% check_der(@(s) df2(t,s), @(s) dsdf2(t,s), 'angle')

ddf1=@(t,s) v1(t,s)'*DLog1(t,s)*v1(t,s);
ddf2=@(t,s) v2(t,s)'*DLog2(t,s)*v2(t,s);

% check_der(@(t) df1(t,s), @(t) ddf1(t,s),'angle')
% check_der(@(t) df2(t,s), @(t) ddf2(t,s),'angle')

ddTotf1=@(s) dsdf1(t(s),s)+ddf1(t(s),s)*dt(s);
%ddTotf1b=@(s) dsLog1(t(s),s)'*v1(t(s),s)+Log1(t(s),s)'*dsv1(t(s),s)+v1(t(s),s)'*DLog1(t(s),s)*v1(t(s),s)*dt(s);
ddTotf1b=@(s) (DLog1(t(s),s)*(-Rb1r(t(s),s)'*Ra1(s)*vVeca1+vVecb1+v1(t(s),s)*dt(s)))'*v1(t(s),s)...
    +Log1(t(s),s)'*hat(v1(t(s),s))*vVecb1;
ddTotf2=@(s) dsdf2(t(s),s)+ddf2(t(s),s)*dt(s);


%check_der(@(s) df(t(s),s), @(s) ddTotf1b(s)+ddTotf2(s))
%plotfun(tOpt,'angle')
dtOptDen=@(s) (ddf1(tOpt(s),s)+ddf2(tOpt(s),s));
%dtOpt=@(s) -(dsdf1(tOpt(s),s)+dsdf2(tOpt(s),s))/dtOptDen(s);
dtOptb=@(s) 1/dtOptDen(s)*e3'*...
    (Rb1(s)*(DLog1(tOpt(s),s)*Razb1(tOpt(s),s)'*vVeca1+(-DLog1(tOpt(s),s)+hat(Log1(tOpt(s),s)))*vVecb1)...
    +Rb2(s)*(DLog2(tOpt(s),s)*Razb2(tOpt(s),s)'*vVeca2+(-DLog2(tOpt(s),s)+hat(Log2(tOpt(s),s)))*vVecb2));
%check_der(tOpt,dtOptb,'angle')

dLog1=@(s) dsLog1(t(s),s)+DLog1(t(s),s)*v1(t(s),s)*dt(s);
%check_der(@(s) Log1(t(s),s), dLog1)

Log1Opt=@(s) Log1(tOpt(s),s);
Log2Opt=@(s) Log2(tOpt(s),s);

dLog1Opt=@(s) dsLog1(tOpt(s),s)+DLog1(tOpt(s),s)*v1(tOpt(s),s)*dtOpt(s);
dLog2Opt=@(s) dsLog2(tOpt(s),s)+DLog2(tOpt(s),s)*v2(tOpt(s),s)*dtOpt(s);

D1=@(s) DLog1(tOpt(s),s);
D2=@(s) DLog2(tOpt(s),s);
A11=@(s) 1/dtOptDen(s)*e3*e3'*Rb1(s)*D1(s);
A12=@(s) 1/dtOptDen(s)*e3*e3'*Rb1(s)*(-D1(s)+hat(Log1(tOpt(s),s)));
A21=@(s) 1/dtOptDen(s)*e3*e3'*Rb2(s)*D2(s);
A22=@(s) 1/dtOptDen(s)*e3*e3'*Rb2(s)*(-D2(s)+hat(Log2(tOpt(s),s)));

% dLog1Opbt=@(s) D1(s)*(...
%     +(-eye(3)+Rb1(s)'*A11(s))*Razb1(tOpt(s),s)'*vVeca1...
%     +( eye(3)+Rb1(s)'*A12(s))*vVecb1...
%     +Rb1(s)'*A21(s)*Razb2(tOpt(s),s)'*vVeca2...
%     +Rb1(s)'*A22(s)*vVecb2...
%     );
% 
% check_der(Log1Opt,dLog1Opbt)

% dLog2Opbt=@(s) D2(s)*(...
%     +Rb2(s)'*A11(s)*Razb1(tOpt(s),s)'*vVeca1...
%     +Rb2(s)'*A12(s)*vVecb1...
%     +(-eye(3)+Rb2(s)'*A21(s))*Razb2(tOpt(s),s)'*vVeca2...
%     +( eye(3)+Rb2(s)'*A22(s))*vVecb2...
%     );
% 
% check_der(Log2Opt,dLog2Opbt)

DLog=@(s) blkdiag(D1(s),D2(s))*...
    [(-eye(3)+Rb1(s)'*A11(s))*Razb1(tOpt(s),s)' (eye(3)+Rb1(s)'*A12(s))...
        Rb1(s)'*A21(s)*Razb2(tOpt(s),s)' Rb1(s)'*A22(s);
     Rb2(s)'*A11(s)*Razb1(tOpt(s),s)' Rb2(s)'*A12(s)...
        (-eye(3)+Rb2(s)'*A21(s))*Razb2(tOpt(s),s)' ( eye(3)+Rb2(s)'*A22(s))];

Log=@(s) [Log1Opt(s); Log2Opt(s)];
dLog=@(s) DLog(s)*[vVeca1; vVecb1; vVeca2; vVecb2];

check_der(Log, dLog,linspace(0,1,10))
