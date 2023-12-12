 clear all
 close all
t1 = 60;
t2 = -20;
A1 = [cosd(t1) -sind(t1);sind(t1) cosd(t1)];
b1 = -2*[1;1];

A2 = [cosd(t2) -sind(t2);sind(t2) cosd(t2)];
b2 = -5*[1;1];


A3 = [1 1];
b3 = 1;

z1 = [0 0;0 1]; %z=[0;1]
z2 = [1 0;0 1]; %z=[1;1]
z3 = 1;

A = (z2*A2*(z1*A1));
b = z2*A2*z1*b1+z2*b2;
%y_out = Ax+b

x = -10:1:10;

nbColors=10;
colors=parula(nbColors);
cnt=0;
funImage(x,x,@(x) A(1,:)*x+b(1),...
                'method','surf','methodOpts',{'FaceColor',...
                colors(mod(cnt-1,nbColors)+1,:)});      
hold on
cnt=3;
funImage(x,x,@(x) A(2,:)*x+b(2),...
                'method','surf','methodOpts',{'FaceColor',...
                colors(mod(cnt-1,nbColors)+1,:)});
hold on
cnt=7;
funImage(x,x,@(x) 0,...
                'method','surf','methodOpts',{'FaceColor',...
                colors(mod(cnt-1,nbColors)+1,:)});
           


