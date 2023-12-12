clear 
close all
x = 0;
y = 0;
u1 = 10;
v1 = 5;
u2 = 10;
v2 = 10;
s = 0.5;
figure(1)
plot([x,u1],[y,v1],'m')
hold on
plot(s*[x,u1],s*[y,v1],'r')
hold on
plot([x,u2],[y,v2],'b-')
hold on
plot([x,u1+u2],[y,v1+v2],'k-')
hold on
plot([x,s*u1+u2],[y,s*v1+v2],'c-')
axis equal

grid on 
