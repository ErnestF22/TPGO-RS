function POCRotationZAxisAlignVectors
a=sphere_randn();
b=sphere_randn();
ap=a(1:2);
bp=b(1:2);
Bp=[b(1) b(2); -b(2) b(1)];
cf=(a(3)-b(3))^2/2;
m=Bp*ap;

R=@(t) rot([0;0;1],t);
Rp=@(t) [cos(t) -sin(t); sin(t) cos(t)];
cs=@(t) [cos(t);sin(t)];
dcs=@(t) [-sin(t); cos(t)];

%check_der(cs,dcs,'angle')

f=@(t) norm(a-R(t)*b)^2/2;
fb=@(t) norm(ap-Rp(t)*bp)^2/2+cf;
fc=@(t) -cs(t)'*m +bp'*bp/2+ap'*ap/2+cf;

csExp=Bp\ap;
%tExp=atan2(csExp(2),csExp(1));
tExp=atan2(m(2),m(1));

funCompare(f,fc,'angle')
hold on
plot(tExp,f(tExp),'bo');
