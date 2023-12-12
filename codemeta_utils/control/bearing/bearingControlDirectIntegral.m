function dz=bearingControlDirectIntegral(z,y,yg,funs)
kp=0.9;%1
ki=0.1;%.01


d=size(z,1)/2;
x=z(1:d,:);
xi=z(d+1:end,:);

g=bearingControlDirect(y,yg,funs);
dx=kp*g+ki*xi;
dxi=g;

dz=[dx;dxi];
