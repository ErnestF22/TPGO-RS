W=sym('W',[3 1]);
Wd=sym('Wd',[3 1]);
dWd=sym('dWd',[3 1]);
Re=sym('Re',[3 3]);
J=sym('J',[3 3]);
eW=W-Re'*Wd;
M=hat(Re'*Wd)*J*Re'*Wd+J*Re'*dWd;

JeW=-hat(W)*J*W+M+J*hat(eW)*Re'*Wd;
d=(2*J-trace(J)*eye(3))*Re'*Wd;
JeW2=

