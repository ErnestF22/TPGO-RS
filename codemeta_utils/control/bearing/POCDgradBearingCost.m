function POCDgradBearingCost
funs=bearingCostFunctions('cosineSq');

%Reshaping function and its derivatives
f=funs.f;
df=funs.df;
ddf=funs.ddf;

xi=randn(2,1);
xg=randn(2,1);

dx=randn(2,1);
tx=@(t) xg+t*dx;


ri=@(x) norm(xi-x);
yi=@(x) (xi-x)/ri(x);
ygi=yi(xg);

[tygi,dtygi]=sphere_randGeodFun(ygi);

%Inner product
ci=@(x,ygi) yi(x)'*ygi;

Dri=@(x) -yi(x);

tri=@(t) ri(tx(t));
dtri=@(t) dx'*Dri(tx(t));
%check_der(tri,dtri)

Pyi=@(x) eye(2)-yi(x)*yi(x)';
Dyi=@(x) -Pyi(x)/ri(x);

tyi=@(t) yi(tx(t));
dtyi=@(t) Dyi(tx(t))*dx;
%check_der(tyi, dtyi)

Dci=@(x,ygi) Dyi(x)'*ygi;

tci=@(t) ci(tx(t),ygi);
dtci=@(t) Dci(tx(t),ygi)'*dx;
%check_der(tci, dtci)

%Function
phi=@(x,ygi) ri(x)*f(ci(x,ygi));


%Gradient w.r.t. location x
%Dphi=@(x) f(ci(x))*Dri(x)+ri(x)*df(ci(x))*Dci(x);
%Dphi=@(x,ygi) -f(ci(x,ygi))*yi(x)-df(ci(x,ygi))*Pyi(x)*ygi;
Dphi=@(x,ygi) (df(ci(x,ygi))*ci(x,ygi)-f(ci(x,ygi)))*yi(x)-df(ci(x,ygi))*ygi;

tphi=@(t) phi(tx(t),ygi);
dtphi=@(t) Dphi(tx(t),ygi)'*dx;
%check_der(tphi,dtphi)

% Dphi1=@(x) (df(ci(x))*ci(x)-f(ci(x)))*yi(x)
% DDphi1=@(x) (ddf(ci(x))*ci(x)+df(ci(x))-df(ci(x)))*yi(x)*Dci(x)'+(df(ci(x))*ci(x)-f(ci(x)))*Dyi(x);
% tDphi1=@(t) Dphi1(tx(t));
% dtDphi1=@(t) DDphi1(tx(t))*v;
% Dphi2=@(x) -df(ci(x))*ygi;
% DDphi2=@(x) -ddf(ci(x))*ygi*Dci(x)';
% tDphi2=@(t) Dphi2(tx(t));
% dtDphi2=@(t) DDphi2(tx(t))*v;
% check_der(tDphi1,dtDphi1)
% check_der(tDphi2,dtDphi2)

%Differential of the gradient w.r.t. location x
%DDphi=@(x) (ddf(ci(x))*ci(x)+df(ci(x))-df(ci(x)))*yi(x)*Dci(x)'+(df(ci(x))*ci(x)-f(ci(x)))*Dyi(x)...
%    -ddf(ci(x))*ygi*Dci(x)';
%DDphi=@(x) (ddf(ci(x))*ci(x)*yi(x)-ddf(ci(x))*ygi)*Dci(x)'+(df(ci(x))*ci(x)-f(ci(x)))*Dyi(x);
%DDphi=@(x) (ddf(ci(x))*(ci(x)*yi(x)-ygi)*ygi'+(df(ci(x))*ci(x)-f(ci(x)))*eye(2))*Dyi(x);
DDphi=@(x,ygi) -1/ri(x)*(ddf(ci(x,ygi))*(ci(x,ygi)*yi(x)-ygi)*ygi'+(df(ci(x,ygi))*ci(x,ygi)-f(ci(x,ygi)))*eye(2))*Pyi(x);

txDphi=@(t) Dphi(tx(t),ygi);
dtxDphi=@(t) DDphi(tx(t),ygi)*dx;
%check_der(txDphi,dtxDphi)

%Differential of the gradient w.r.t. measurement ygi
dfci=@(x,ygi) df(ci(x,ygi));
ddfci=@(x,ygi) ddf(ci(x,ygi));
DygDphi=@(x,ygi) -(df(ci(x,ygi))*eye(2)+ddf(ci(x,ygi))*Pyi(x)*ygi*yi(x)');
%DygDphi=@(x,ygi) bearingCostGeneral_termsDygrad(yi(x),ygi,dfci(x,ygi),ddfci(x,ygi));
tgDphi=@(t) Dphi(xg,tygi(t));
dtgDphi=@(t) DygDphi(xg,tygi(t))*dtygi(t);

%check_der(tgDphi,dtgDphi)

tgxDphi=@(t) Dphi(tx(t),tygi(t));
dtgxDphi=@(t) DDphi(tx(t),tygi(t))*dx+DygDphi(tx(t),tygi(t))*dtygi(t);
check_der(tgxDphi,dtgxDphi)
