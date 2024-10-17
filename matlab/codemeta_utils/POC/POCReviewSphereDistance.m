function POCReviewSphereDistance
k1i=sphere_randn();
k12=sphere_randn();
t1i=pi*rand();
t12=pi*rand();
G1=sphere_randn();

y1i=rot_exp(eye(3),-t1i*hat3(k1i))*G1;
y12=rot_exp(eye(3),-t12*hat3(k12))*G1;

cd=cos(sphere_dist(y1i,y12));

f=cos(t1i)*cos(t12)+sin(t1i)*sin(t12)*k1i'*k12;

disp([cd f])

Gi=sphere_randn([1;0;0],0.5);
xGj=sphere_randn([1;0;0],0.5);
Gk=sphere_randn([1;0;0],0.5);
tij=acos(Gi'*Gj);
tik=acos(Gi'*Gk);
tjk=acos(Gj'*Gk);
%kij=cross(Gi)*Gj/sin(tij);
kik=cross(Gi,Gk)/sin(tik);
kjk=cross(Gj,Gk)/sin(tjk);

cd=cos(tij);
f=cos(tik)*cos(tjk)+sin(tik)*sin(tjk)*kik'*kjk;

disp([cd f])

kik=sphere_randn();
kjk=sphere_randn();
Gib=rot_exp(eye(3),-tik*hat3(kik))*Gk;
Gjb=rot_exp(eye(3),-tjk*hat3(kjk))*Gk;
tikb=sphere_dist(Gib,Gk);
tjkb=sphere_dist(Gjb,Gk);

cd=cos(sphere_dist(Gib,Gjb));
f=cos(tikb)*cos(tjkb)+sin(tikb)*sin(tjkb)*kik'*kjk;

disp([cd f])



%keyboard
