%function derivationEqInfinitesimalPlanarMotion
w=sym('w',[3 1]); %angular velocity
t=sym('t',[3 1]); %translation velocity
syms x y %image coordinates
syms a b c %plane parameters
syms X Y Z %3-D coordinates

P=[X;Y;Z];
Pdot=hat(w)*P+t;
p=[x;y;1];

uvA=1/Z^2*(Pdot*P(3)-P*Pdot(3));

sr=[1 0 -x; 0 1 -y; 0 0 0];
uvB=sr*(1/Z*t+hat(w)*p);

%check that the two expressions for the flow uv are the same
disp(simplify(uvA-subs(subs(uvB,x,X/Z),y,Y/Z)))

%plane equation
n=[a;b;c];
%n'*P=0 -- divide both sides by Z --> n.'*p=1
uvC=subs(uvB,1/Z,n.'*p);
disp(uvC)

syms Ix Iy It  %image derivatives
eq=expand([Ix Iy It]*[uvC(1:2); 0]);

a1=a*t(1);
a2=a*t(2);
a3=a*t(3);
b1=b*t(1);
b2=b*t(2);
b3=b*t(3);
c1=c*t(1);
c2=c*t(2);
c3=c*t(3);
