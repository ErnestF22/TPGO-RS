function POCrotGeodesicProjectionAlgebraic
R=rot_randn();
Rz0=rot_randn();
uz=sphere_randn();
vRz=rot_hat(Rz0,uz);
Rz=@(t) rot_exp(Rz0,t*vRz);

%verify that Householder rotation gives the rotation in Img(R'*Rz) which is
%pi from the identity
[ROpt,tOpt]=rot3_projectOnGeodesicTrigonometric(R,Rz0,vRz);
tPi=modAngle(tOpt+pi);

H=householderRotation(uz,R'*Rz0*uz);

disp('pi-rot_dist(R,Rz(tPi))')
disp(pi-rot_dist(R,Rz(tPi)))

disp('[R''*Rz(tPi) H]')
disp([R'*Rz(tPi) H])

%use Householder rotation to get ROpt
Ruzpi=2*(uz*uz')-eye(3);
disp([ROpt R*H*Ruzpi])

