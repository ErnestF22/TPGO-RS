function R=quat2rot(q)

w=q(1,:);
x=q(2,:);
y=q(3,:);
z=q(4,:);

Nq = w.^2 + x.^2 + y.^2 + z.^2;
idx=Nq>0.0;
s(idx)=2./Nq;
s(~idx)=0;
X = x.*s; Y = y.*s; Z = z.*s;
wX = shiftdim(w.*X,-1); wY = shiftdim(w.*Y,-1); wZ = shiftdim(w.*Z,-1);
xX = shiftdim(x.*X,-1); xY = shiftdim(x.*Y,-1); xZ = shiftdim(x.*Z,-1);
yY = shiftdim(y.*Y,-1); yZ = shiftdim(y.*Z,-1); zZ = shiftdim(z.*Z,-1);
R=[ 1.0-(yY+zZ)   xY-wZ      xZ+wY;...
    xY+wZ      1.0-(xX+zZ)   yZ-wX;...
    xZ-wY         yZ+wX   1.0-(xX+yY)];
