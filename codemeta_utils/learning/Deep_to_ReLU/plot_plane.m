
figure(1)
hold on

[X,Y] = meshgrid(-100:5:100,-100:5:100);
% 
% Z1 = A1(1,1)*X+A1(1,2)*Y+b1(1);
% plot3(X,Y,Z1,'r')
% hold on
% Z2 = A1(2,1)*X+A1(2,2)*Y+b1(2);
% plot3(X,Y,Z2,'g')
% hold on
% 
a = A2*A1;
b = A2*b1+b2;
% 
% Z3 = a(1,1)*X+a(1,2)*Y+b(1);
% plot3(X,Y,Z3,'b')
% hold on
Z4 = a(2,1)*X+a(2,2)*Y+b(2);
surf(X,Y,Z4)
hold on

a = A3*A2*A1;
b = A3*A2*b1+A3*b2+b3;

Z5 = a(1,1)*X+a(1,2)*Y+b(1);
surf(X,Y,Z5)
