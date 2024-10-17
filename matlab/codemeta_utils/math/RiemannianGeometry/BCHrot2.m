function [w,alpha,beta,gamma]=BCHrot2(u,v)

theta=norm(u);
phi=norm(v);

angleuv=u'*v/theta/phi;

a=sin(theta)*cos(phi/2)^2-sin(phi)*sin(theta/2)^2*cos(angleuv);
b=sin(phi)*cos(theta/2)^2-sin(theta)*sin(phi/2)^2*cos(angleuv);
c=0.5*sin(theta)*sin(phi)-2*sin(theta/2)^2*sin(phi/2)^2*cos(angleuv);

d=sqrt(a^2+b^2+2*a*b*cos(angleuv)+c^2*sin(angleuv)^2);

alpha=asin(d)/d*a/theta;
beta=asin(d)/d*b/phi;
gamma=asin(d)/d*c/theta/phi;

w=alpha*u+beta*v+gamma*cross(u,v);

z=computez(theta,phi,u,v);
hat(z)
w=sin(norm(z))/norm(z)*z;
;

function z=computez(theta,phi,u,v)
uhat=hat(u);
vhat=hat(v);
stheta=sin(theta)/theta;
sphi=sin(phi)/phi;
s2theta=sin(theta/2)^2/theta^2;
s2phi=sin(phi/2)^2/phi^2;
z=vee(stheta*uhat+sphi*vhat+0.5*stheta*sphi*[uhat*vhat-vhat*uhat]+stheta*s2phi*(uhat*vhat^2+vhat^2*uhat)+sphi*s2theta*(uhat^2*vhat+vhat*uhat^2)+2*s2theta*s2phi*(uhat^2*vhat^2-vhat^2*uhat^2));