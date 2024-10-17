function POCDerivativeHouseholder

[T,dt]=real_randGeodFun(randn(3,1));
nv=@(t) norm(T(t));
v=@(t) cnormalize(T(t));
dv=@(t) (eye(3)-v(t)*v(t)')/nv(t)*dt(t);

% [v,dv]=sphere_randGeodFun(sphere_randn(eye(3,1)));
% check_der(v,dv)

e3=[0;0;1];
 
nd=@(t) norm(v(t)+e3);
d=@(t) cnormalize(v(t)+e3);
dd=@(t) (eye(3)-d(t)*d(t)')/nd(t)*dv(t);
%check_der(d,dd,'angle')

R0=@(t) 2*(d(t)*d(t)')-eye(3);
dR0=@(t) 2*(dd(t)*d(t)'+d(t)*dd(t)');
dR0b=@(t) 2*((eye(3)-d(t)*d(t)')/nd(t)*dv(t)*d(t)'+d(t)*dv(t)'*(eye(3)-d(t)*d(t)')/nd(t));

% plotfun(@(t) reshape(R0(t)'*R0(t),[],1))
% plotfun(@(t) det(R0(t)),'angle')

check_der(R0,dR0b,'angle')
