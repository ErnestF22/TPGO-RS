function rigidTransform_testVelocities
R=rot_randn();
T=randn(3,1);
[Rt,~,~,~,w]=rot_randGeodFun(R);
[Tt,~,~,v]=real_randGeodFun(T);
wv=[w;v];

%Check world-velocities for "pose" representation
Gt=@(t) RT2G(Rt(t),Tt(t),'compact');
dGt=@(t) rot3r3_hat(Gt(t),wv); 
funCheckDer(Gt,dGt)

%Check camera-velocities from world-velocities for "pose" representation
invGt=@(t) invg(Gt(t));
wvInvt=@(t) rigidTransformG(Gt(t),[w;v],'poses','wc','velocities');
dinvGt=@(t) rot3r3_hat(invGt(t),wvInvt(t));
funCheckDer(invGt,dinvGt)

%Check camera-velocities from world-velocities  for "reference" representation
wvInvtB=@(t) rigidTransformG(invGt(t),[w;v],'references','wc','velocities');
dinvGtB=@(t) rot3r3_hat(invGt(t),wvInvtB(t));
funCheckDer(invGt,dinvGtB)

%Check world-velocities from camera-velocities for "pose" representation
wvtB=@(t) rigidTransformG(Gt(t),wvInvt(t),'poses','cw','velocities');
dGtB=@(t) rot3r3_hat(Gt(t),wvtB(t));
funCheckDer(Gt,dGtB)

%Check world-velocities from camera-velocities for "reference" representation
wvtC=@(t) rigidTransformG(invGt(t),wvInvt(t),'references','cw','velocities');
dGtC=@(t) rot3r3_hat(Gt(t),wvtC(t));
funCheckDer(Gt,dGtC)
