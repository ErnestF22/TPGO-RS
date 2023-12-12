close all
clear 

t1 = -30;
t2 = 45;

A1 = [cosd(t1) -sind(t1);sind(t1) cosd(t1)];
b1 = [1;1];

A2 = [cosd(t2) -sind(t2);sind(t2) cosd(t2)];
b2 = [1;1];

A3 = [1 1];
b3 = 1;

% x is an 2D input 
x = [5;1];

% -first layer:
h1 = A1*x+b1; % where A_1 is a 2*2 matrix and b_1 is 2*1
z1 = (h1>0); % determine z
y1 = A1.*z1*x+(b1.*z1) % where y_1 is 2*1 and ReLU is an elementwise function
% -second layer:
h2 = A2*y1+b2; % where A_2 is a 2*2 matrix and b_2 is 2*1
z2 = (h2>0); % determine z
y2 = A2.*z2*y1+(b2.*z2) % where y_2 is 2*1 and ReLU is an elementwise function
% -output:
h3 = A3*y2+b3; % where A_3 is a 2*1 matrix and b_3 is scalar
z3 = (h3>0); % determine z
y3 = A3.*z3*y2+(b3.*z3) % where y_3 is scalar

[y1,y2,y3] = netForward(x,A1,b1,A2,b2,A3,b3)