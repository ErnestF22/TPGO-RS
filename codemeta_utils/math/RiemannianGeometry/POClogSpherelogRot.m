function POClogSpherelogRot
u2=sphere_randn();
[u,du]=sphere_randGeodFun();

%funCheckDer(u,du)

%Ropt=@(t) householderRotation3Min(u(t),u2);
H=@(t) householderRotation(u(t),u2);
v=@(t) u(t)+u2;
vp=@(t) cnormalize(v(t));
DH=@(t) -2*hat(vp(t))/norm(v(t));%*Dvp(t);

%vdH=@(t) DH(t)*du(t);
%dH=@(t) rot_hat(H(t),vdH(t));
%funCheckDer(H,dH,'angle')
Rpiu2=2*(u2*u2')-eye(3);
Ropt=@(t) Rpiu2*H(t);
RoptRef=@(t) 2*u2*u(t)'-H(t);
%funCompare(Ropt,RoptRef,'angle')

distSphere=@(t) sphere_dist(u(t),u2)^2/2;
distRot=@(t) rot_dist(eye(3),Ropt(t))^2/2;
%funCompare(distSphere,distRot,'angle')

logSphere=@(t) sphere_log(u(t),u2);
ddistSphere=@(t) -logSphere(t)'*du(t);
%funCheckDer(distSphere,ddistSphere,'angle')

%dRopt=@(t) Rpiu2*dH(t);
vdRopt=@(t) DH(t)*du(t);
%dRopt=@(t) rot_hat(Ropt(t),vdRopt(t));
%funCheckDer(Ropt,dRopt,'angle')

logRot=@(t) logrot(Ropt(t));
ddistRot=@(t) -2*logRot(t)'*hat(vp(t))/norm(v(t))*du(t);
%funCheckDer(distRot,ddistRot,'angle')
%funCompare(ddistSphere,ddistRot,'angle')
ddistRotSphere=@(t) -logrot(Ropt(t))'*hat(u(t))*du(t);
%funCompare(ddistSphere,ddistRotSphere,'angle')

logSphereRot=@(t) -hat(u(t))*logRot(t);
funCompare(logSphere,logSphereRot,'angle')

I=eye(3);
Pu=@(t) (I-u(t)*u(t)');
M=@(t) 2*(hat(u(t))+Pu(t)'*hat(u2)*Pu(t))/norm(v(t))^2;
%funCompare(M,@(t) hat(u(t)),'angle')

%keyboard


function DH=diffH(u1,u2)
v=u1+u2;
[Dvp,vp]=cnormalizeDiffMat(v);
DH=-2*hat(vp)*Dvp;
