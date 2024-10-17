function POCDerRadialSpeed
d=3;
[x,~,~,dx]=real_randGeodFun(randn(d,1));

r=@(t) norm(x(t));
v=@(t) x(t)/r(t);
s=@(t) dx'*v(t);

%funCheckDer(r,s)

I=eye(d);
ds=@(t) 1/r(t)*dx'*(I-v(t)*v(t)')*dx;

%funCheckDer(s,ds)

funPlot(r,linspace(0,100))
