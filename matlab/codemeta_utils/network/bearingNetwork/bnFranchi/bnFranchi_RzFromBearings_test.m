function bnFranchi_RzFromBearings_test
%%
%resetRands()
xi=randn(3,1);
xj=[1;0;0];
Ri=rot([0;0;2*pi*rand]);
Rj=rot([0;0;2*pi*rand]);

bij=bnFranchi_computeBearing(Ri,xi,xj);
bji=bnFranchi_computeBearing(Rj,xj,xi);

Rij=bnFranchi_RzFromBearings(bij,bji);

disp(norm(Ri'*Rj-Rij,'fro'))
%%
