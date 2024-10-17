function POCLinearOptimalControl
A=[0 1;0 0];
B=[0;1];
syms theta s t
x0=[0;0];
x1=[cos(theta); sin(theta)];
Phits=expm(A*(t-s));
Phi0s=subs(Phits,t,0);
Phi01=subs(Phits,{t,s},{0,1});
Phi10=subs(Phits,{t,s},{1,0});
Wts=int((Phits*B)*(Phits*B).',s);
W01=subs(Wts,{t s},{0,1});
lambda0=W01\(x0-Phi01*x1);
us=-B.'*Phi0s.'*lambda0;
J=simplify(lambda0.'*W01*lambda0);
%x1b=Phi10*x0+int(subs(Phits,t,0)\B*us,s,0,t);
x1b=Phi10*x0+Phi10*int(Phi0s*B*us,s,0,t);
Jb=int(us.'*us,s,0,1);
subs(x1b,t,1)
simplify(J-Jb)
% [tn,xn]=ode45(@(t,x) A*x+B,[0 1],x0);
% xn=xn';
% disp(xn(:,end))
keyboard
