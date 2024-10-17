function bnFranchi_Gamma_test
%%
xl=randn(3,1);
xm=randn(3,1);
xn=randn(3,1);
Rl=rot([0;0;2*pi*rand]);
Rm=rot([0;0;2*pi*rand]);
Rn=rot([0;0;2*pi*rand]);

bmn=bnFranchi_computeBearing(Rm,xm,xn);
dmn=bnFranchi_computeDistance(xm,xn);
bml=bnFranchi_computeBearing(Rm,xm,xl);
dml=bnFranchi_computeDistance(xm,xl);
bln=bnFranchi_computeBearing(Rl,xl,xn);
dln=bnFranchi_computeDistance(xl,xn);
blm=bnFranchi_computeBearing(Rl,xl,xm);
dlm=bnFranchi_computeDistance(xl,xm);

bnm=bnFranchi_computeBearing(Rn,xn,xm);
dnm=bnFranchi_computeDistance(xn,xm);
bnl=bnFranchi_computeBearing(Rn,xn,xl);
dnl=bnFranchi_computeDistance(xn,xl);

Rln=Rl'*Rn;
Rmn=Rm'*Rn;
Rlm=Rl'*Rm;

% sln=sqrt(1-bearingComputeCosine(bmn,bml).^2);
% %sml=sqrt(1-bearingComputeCosine(bnl,bnm).^2);
% %sml=sqrt(1-bearingComputeCosine(Rln*bnl,Rln*bnm).^2);
% %sml=sqrt(1-bearingComputeCosine(-bln,Rln*bnm).^2);
% %sml=sqrt(1-bearingComputeCosine(bln,-Rlm*Rmn*bnm).^2);
% %sml=sqrt(1-bearingComputeCosine(bln,Rlm*bmn).^2);
% sml=sqrt(1-bearingComputeCosine(bln,bnFranchi_RzFromBearings(blm,bml)*bmn).^2);
% disp([dln/dml sln/sml])

glmn=bnFranchi_Gamma(bmn,bml,blm,bln);
disp([dln/dml glmn])
