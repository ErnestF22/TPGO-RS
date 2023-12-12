function bnFranchi_computeBearing_test
%%
%resetRands();
xi=randn(3,1);
xj=randn(3,1);
Ri=rot([0;0;2*pi*rand]);
Rj=rot([0;0;2*pi*rand]);


bij=bnFranchi_computeBearing(Ri,xi,xj);
bji=bnFranchi_computeBearing(Rj,xj,xi);

disp([Ri*bij -Rj*bji])
%%