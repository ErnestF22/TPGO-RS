function POCRelationDistanceSphereSO3
k=3;
I=eye(3);
e3=I(:,k);
ev=sphere_randn();

U=rot_geodFun(I,rot_hat(I,e3));
H=diag([-1 -1 1])*householderRotation(ev,k);
R=@(t) U(t)*H;
f=@(t) rot_dist(I,R(t));

dSphere=sphere_dist(e3,ev);
plotfun(f,'angle')
hold on
plot([-pi pi],[dSphere dSphere],'r:')
hold off

