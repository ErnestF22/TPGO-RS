function [X,Y,Z]=meshgrid3(x,y,z);
x=shiftdim(x);
y=shiftdim(y);
z=shiftdim(shiftdim(z),-2);

d1=length(x);
d2=length(y);
d3=length(z);

X=repmat(x',[d2, 1, d3]);
Y=repmat(y,[1, d1, d3]);
Z=repmat(z,[d2,d2,1]);
