function POCJacobianBearingConstraint

xi0=randn(2,1);
vi=randn(2,1);
xj0=randn(2,1);
vj=randn(2,1);
xi=@(t) xi0+t*vi;
xj=@(t) xj0+t*vj;
dxi=@(t) vi;
dxj=@(t) vj;
[R,v,R0,v0]=rot_randGeodFun(eye(2),'speed',rand);
[Q,w,Q0,w0]=rot_randGeodFun(eye(3),'speed',rand);
alpha0=rot_vee(R0,v0);
beta0=rot_vee(Q0,w0);

n=@(t) norm(xi(t)-xj(t));
dn=@(t) (xi(t)-xj(t))'*(dxi(t)-dxj(t))/n(t);

%check_der(n,dn);

u1=randn(2,1);
u2=randn(2,1);
p1=randn(3,1);
p2=randn(3,1);

f=@(t) (xi(t)-xj(t))/n(t);
Jf=@(t) (eye(2)*n(t)-(xi(t)-xj(t))*(xi(t)-xj(t))'/n(t))/(n(t)^2);
df=@(t) [Jf(t) -Jf(t)]*[dxi(t); dxj(t)];

%check_der(f,df)

f2=@(t) u1'*R(t)*u2;
df2=@(t) alpha0*u1'*R(t)*[0 -1; 1 0]*u2;

%check_der(f2,df2,'angle')

f3=@(t) p1'*Q(t)*p2;
df3=@(t) -p1'*Q(t)*hat(p2)*beta0;

check_der(f3, df3,'angle')