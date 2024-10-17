function spy3(m)
x=1:size(m,1);
y=1:size(m,2);
z=1:size(m,3);

[X,Y,Z]=meshgrid3(x,y,z);

idx=m~=0;

plot3(X(idx),Y(idx),Z(idx),'.');
