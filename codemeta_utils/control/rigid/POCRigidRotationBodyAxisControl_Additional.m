function POCRigidRotationBodyAxisControl_Additional
%Check derivative used in POCRigidRotationBodyAxisControl.m
e3=[0;0;1];
[er,der,dder]=POCRigidRotationBodyAxisControl_reference(0.5);

[R,~,~,~,w]=rot_randGeodFun();

ev=@(R,er) R'*er;
dev=@(R,w,er,der) hat(ev(R,er))*w+R'*der;

RMin=@(R,er) householderRotation3Min(ev(R,er),3);
DRMin=@(R,er) rot(pi*e3)'*householderRotation_DiffMat(ev(R,er),3);
dRMinVec=@(R,w,er,der) DRMin(R,er)*dev(R,w,er,der);

LogRMin=@(R,er) logrot(RMin(R,er));
DLogRMin=@(R,er) rot3_logDiff(eye(3),LogRMin(R,er),'tangentVec');
dLogRMin=@(R,w,er,der) DLogRMin(R,er)*dRMinVec(R,w,er,der);

evt=@(t) ev(R(t),er(t));
devt=@(t) dev(R(t),w,er(t),der(t));
%check_der(evt,devt)

RMint=@(t) RMin(R(t),er(t));
%dRMint=@(t) RMint(t)*hat(DRMin(R(t),er(t))*devt(t));
dRMint=@(t) RMint(t)*hat(dRMinVec(R(t),w,er(t),der(t)));
%check_der(RMint,dRMint)

LogRMint=@(t) LogRMin(R(t),er(t));
dLogRMint=@(t) dLogRMin(R(t),w,er(t),der(t));
%check_der(LogRMint,dLogRMint)

check_der(@(t) Rpihalf(er(t)),@(t) dRpihalf(er(t),der(t)))

function R=Rpihalf(v)
R=hat(v)+v*v';   

function dR=dRpihalf(v,dv)
dR=hat(dv)+dv*v'+v*dv';   
