function POCDerivativeProjector
d=8;
[x,dx]=real_randGeodFun(randn(d,1));

nx=@(t) norm(x(t));

w=@(t) x(t)/nx(t);
dw=@(t) (eye(d)-w(t)*w(t)')*dx(t)/nx(t);

P=@(t) eye(d)-w(t)*w(t)';
dP=@(t) -dw(t)*w(t)'-w(t)*dw(t)';

%check_der(P,dP)

f=@(t) 1/2*dx(t)'*dP(t)*dx(t);
fb=@(t) -(w(t)'*dx(t))*(dw(t)'*dx(t));
fc=@(t) -1/nx(t)^2*(dx(t)'*(eye(d)-w(t)*w(t)')*dx(t))*(x(t)'*dx(t));

t=linspace(0,1);
plotfun(f,t,'rx')
hold on
plotfun(fc,t)
hold off
%x0=x(0);
%P0=P(0);
%dP0=dP(0);
%disp(eig(dP0))
