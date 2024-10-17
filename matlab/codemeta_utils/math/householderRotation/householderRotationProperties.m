function householderRotationProperties
x1=sphere_randn();
x2=sphere_randn();
Re1pi=rot(pi*x1);
Re2pi=rot(pi*x2);
Re2pihalf=rot(pi/2*x2);

H=householderRotation(x1,x2);

disp('[x1 H*ek]')
disp([x1 H*x2])
disp('[x2 H''*x1]')
disp([x2 H'*x1])

RMin1=H*Re1pi;
RMin2=H*Re2pi;
disp('[sphere_dist(x1,ek) rot_dist(eye(3),RMin2)]')
disp([sphere_dist(x1,x2) rot_dist(eye(3),RMin2)])

disp('[logrot(RMin2) -logrot(RMin1) Rkpihalf*sphere_log(ek,x1)]')
disp([logrot(RMin2) -logrot(RMin1) Re2pihalf*sphere_log(x2,x1)])

disp('[hat(x1)*logrot(RMin2) hat(x1)*logrot(RMin1) sphere_log(x1,ek)]')
disp([hat(x1)*logrot(RMin2) -hat(x1)*logrot(RMin1) sphere_log(x1,x2)])

disp('[-RMin2''*sphere_log(x1,ek) H''*sphere_log(x1,ek) sphere_log(ek,x1)]')
disp([-RMin2'*sphere_log(x1,x2) H'*sphere_log(x1,x2) sphere_log(x2,x1)])
%[H12 H(R*x1,R*x2)]

