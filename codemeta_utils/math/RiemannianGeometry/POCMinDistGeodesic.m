function POCMinDistGeodesic
v=[0;0;1];
U=rot_geodFun(eye(3),rot_hat(eye(3),v));
R=rot_randn();
f=@(t) rot_dist(R,U(t))^2/2;
df=@(t) logrot(R'*U(t))'*v;
skew=@(A) (A-A')/2;
df2=@(t) trace((R'*U(t))'*hat(v));

s=R(2,1)-R(1,2);
c=R(1,1)+R(2,2);
tOpt1=atan2(s,c);
tOpt2=modAngle(tOpt1+pi);

ROpt1=(R'*U(tOpt1));
ROpt2=ROpt1*diag([-1 -1 1]);

ctheta1=trace(ROpt1)
ctheta2=trace(ROpt2)
R(3,3)

I=eye(3);
check_der(f,df,'angle')
hold on
plotfun(df2,'angle','c')
plot([tOpt1 tOpt2],[0 0],'ko')
plot([tOpt1 tOpt2],[rot_dist(I,ROpt1)^2/2 rot_dist(I,ROpt2)^2/2],'mo')
hold off

