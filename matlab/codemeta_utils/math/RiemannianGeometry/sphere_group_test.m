D=4;
y1=randn(D,1); y1=y1/norm(y1);
y2=randn(D,1); y2=y2/norm(y2);
y3=randn(D,1); y3=y3/norm(y3);

e=[1;zeros(D-1,1)];

disp('check closure')
[norm(sphere_comp(y1,y2)) norm(sphere_comp(y2,y3)) norm(sphere_comp(y3,y1))]

disp('check identity')
[y1 sphere_comp(y1,e) sphere_comp(e,y1)]

disp('check associativity')
[sphere_comp(y1,sphere_comp(y2,y3)) sphere_comp(sphere_comp(y1,y2),y3)]

disp('check inverse')
sphere_comp(y1,y1,'inv')
sphere_comp(y1,sphere_comp(y1,e,'inv'))

disp('check pull-back of log')
h1=sphere_log(e,y1);
[sphere_log(y1,y2) sphere_parallel(e,h1,sphere_log(e,sphere_comp(y1,y2,'inv')))]
