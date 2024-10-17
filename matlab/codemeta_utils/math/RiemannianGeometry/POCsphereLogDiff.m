function POCsphereLogDiff

I=eye(3);
[e1,de1]=sphere_randGeodFun();
[e2,de2]=sphere_randGeodFun();

Re1pi=@(t) -I+2*e1(t)*e1(t)';
Re2pi=@(t) -I+2*e2(t)*e2(t)';
%funCompare(@(t) det(Re1pi(t)), @(t) 1)
dRe2pi=@(t) 2*(de2(t)*e2(t)'+e2(t)*de2(t)');
%funCheckDer(Re2pi,dRe2pi,'angle')
dRe2piVec=@(t) 2*hat(de2(t))*e2(t);
%checkDerRot(Re2pi,dRe2piVec)

de=@(t) [de1(t);de2(t)];

H=@(t) householderRotation(e1(t),e2(t));
DH=@(t) householderRotation_DiffMat(e1(t),e2(t));
dHVec=@(t) DH(t)*de(t);
%checkDerRot(H,dHVec)

Rmin=@(t) H(t)*Re2pi(t);
%dRminVec=@(t) Re2pi(t)'*dHVec(t)+dRe2piVec(t);
DRmin=@(t) Re2pi(t)'*DH(t)-2*[zeros(3) hat(e2(t))];
dRminVec=@(t) DRmin(t)*de(t);
%checkDerRot(Rmin,dRminVec)

l=@(t) logrot(Rmin(t));
Dl=@(t) rot3_logDiff(eye(3),l(t),'tangentVec');
%funCheckDer(l, @(t) Dl(t)*dRminVec(t))

f=@(t) sphere_log(e1(t),e2(t));
fb=@(t) hat(e1(t))*l(t);
%funCompare(f,fb)

%df=@(t) -hat(l(t))*de1(t)+hat(e1(t))*Dl(t)*dRminVec(t);
%df=@(t) -hat(l(t))*de1(t)+hat(e1(t))*Dl(t)*DRmin(t)*de(t);
df=@(t) ([-hat(l(t)) zeros(3)]+hat(e1(t))*Dl(t)*DRmin(t))*de(t);
funCheckDer(f,df)


%dHVec=@(t) householderRotation_Diff(e1(t),e2(t),de1(t),de2(t));
%funCheckDer(H, @(t) rot_hat(H(t), dHVec(t)),'angle')

function checkDerRot(R,dRVec)
funCheckDer(R,@(t) rot_hat(R(t),dRVec(t)))
