function [n,Ah,bh,y,Cl,Cb,Wl,Wb,Au,bu,Ax,bx]=data
n=[-1;-1];
A=[-0.5 1;4 1;0.15 1;-2.4 1];
b=[30;-15;4.25;68];
p=[15;15];
[Ax,bx]=convexSet(A,b,p);

y=[20 20;10 25;5 5;25 8]';
Cl=1;
Cb=1;
Wl=1;
Wb=[1;1;1];
Au=[1 0;0 1;-1 0;0 -1];
bu=5*[1;1;1;1];

Ah=-Ax;
bh=bx;
Ah(4,:)=[];
bh(4,:)=[];
end