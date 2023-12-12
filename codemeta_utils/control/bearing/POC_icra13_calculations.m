function POC_icra13_calculations
%xg=sym('xg',[2 1]);
xg=[1;0];
%xi=sym('xi',[2 1]);
xi=[0;0];
v=sym('v',[2 1]);
syms t;
tx=xg+t*v;
ygi=(xi-xg)/sqrt((xi-xg).'*(xi-xg));
tyi=(xi-tx)/sqrt((xi-tx).'*(xi-tx));
tri=sqrt((xi-tx).'*(xi-tx));
tc=simplify(ygi.'*tyi);
disp('tc=')
pretty(tc)
dtc=simplify(diff(tc,t));
disp('dtc=')
pretty(dtc)

vPyiygi=simplify(v.'*(eye(2)-tyi*tyi.')*ygi);
disp('vPyiygi=')
pretty(vPyiygi)

vtyi=v.'*tyi;
dvtyi=simplify(diff(vtyi,t));
disp('dvtyi')
pretty(dvtyi)
