function [w,alpha,beta,gamma]=BCHrot(u1,u2)

theta1=sqrt(sum(u1.^2));
theta2=sqrt(sum(u2.^2));

x1=u1./repmat(theta1,3,1);
x2=u2./repmat(theta2,3,1);


ctheta1=cos(theta1/2);
ctheta2=cos(theta2/2);
stheta1=sin(theta1/2);
stheta2=sin(theta2/2);

x1x2=sum(x1.*x2);
phi=2*acos(ctheta1.*ctheta2-stheta1.*stheta2.*x1x2);

sphi=sin(phi/2);

a=stheta1.*ctheta2./sphi;
b=ctheta1.*stheta2./sphi;
c=2*stheta1.*stheta2./sphi;

alpha=phi.*a./theta1;
beta=phi.*b./theta2;
gamma=phi.*c./theta1./theta2/2;

w=repmat(alpha,3,1).*u1+repmat(beta,3,1).*u2+repmat(gamma,3,1).*cross(u1,u2);

%At this point some of the w have norm > pi and are outside the cut locus
%for the tangent space. We will now fix this by inverting the rotation axis.

normw=sqrt(sum(w.^2));
fixIdx=find(normw>pi);

normw1=pi-mod(normw(fixIdx),pi);

f=normw(fixIdx)./normw1;

w(:,fixIdx)=-w(:,fixIdx)./repmat(f,3,1);

alpha(fixIdx)=alpha(fixIdx).*f;
beta(fixIdx)=beta(fixIdx).*f;
gamma(fixIdx)=gamma(fixIdx).*f;
