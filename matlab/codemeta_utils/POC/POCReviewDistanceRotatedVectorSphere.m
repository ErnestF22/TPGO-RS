function POCReviewDistanceRotatedVectorSphere
v1=real_randGeodFun(randn(3,1),'speed',0);
v2=real_randGeodFun(randn(3,1));

theta1=@(t) norm(v1(t));
theta2=@(t) norm(v2(t));

u1=@(t) v1(t)/theta1(t);
u2=@(t) v2(t)/theta2(t);

gamma=sphere_randn();

R1=@(t) rot(u1(t),theta1(t));
R2=@(t) rot(u2(t),theta2(t));

f=@(t) sphere_dist(R1(t)*gamma,R2(t)*gamma);
fB=@(t) acos(cos(theta1(t))*cos(theta2(t))+sin(theta1(t))*sin(theta2(t))*u1(t)'*u2(t));

funCompare(f,fB,'angle')
