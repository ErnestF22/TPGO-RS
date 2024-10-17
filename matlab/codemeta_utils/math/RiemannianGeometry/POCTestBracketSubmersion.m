function POCTestBracketSubmersion
close all

e3=[0;0;1];
pi=@(R) R*e3;
Dpi=@(R) -R*hat(e3);

%check differential Dpi
[R,dR,~,~,dRVec]=rot_randGeodFun();
piR=@(t) pi(R(t));
dpiR=@(t) Dpi(R(t))*dRVec;

%funCheckDer(piR,dpiR,'angle')

%define function on the sphere and its gradient
A=diag([1,2,3]);
b=[0.5;0;0];
f=@(x) 0.5*x'*A*x+x'*b;
gradf=@(x) (eye(3)-x*x')*(A*x+b);

%check gradient
[x,dx]=sphere_randGeodFun();

fx=@(t) f(x(t));
dfx=@(t) gradf(x(t))'*dx(t);

%funCheckDer(fx,dfx,'angle')

%generate corresponding points
R01=rot_randn();
x0=pi(R01);
R02=R01*blkdiag(rot_randn(eye(2)),1);

%disp([x0 pi(R01) pi(R02)])


%gradient on R
fb=@(R) f(pi(R));
gradfbVec=@(R) pinv(Dpi(R))*gradf(pi(R));
gradfb=@(R) rot_hat(R,gradfbVec(R));

fbR=@(t) fb(R(t));
dfbR=@(t) rot_metric(R(t),gradfb(R(t)),dR(t));

%funCheckDer(fbR,dfbR,'angle')


%generate vectors at x0 and pullbacks at two equivalent rotations R01 and
%R02
vx0a=sphere_randTangentNormVector(x0);
vx0b=sphere_randTangentNormVector(x0);

DR01=Dpi(R01);
DR02=Dpi(R02);
vR01aVec=pinv(DR01)*vx0a;
%disp([vx0a Dpi(R01)*vR01aVec])
vR01bVec=pinv(DR01)*vx0b;
vR02aVec=pinv(DR02)*vx0a;
vR02bVec=pinv(DR02)*vx0b;

vR01a=rot_hat(R01,vR01aVec);
vR01b=rot_hat(R01,vR01bVec);
vR02a=rot_hat(R02,vR02aVec);
vR02b=rot_hat(R02,vR02bVec);

%check that metrics are the same
%disp([sphere_metric(x0,vx0a,vx0b) rot_metric(R01,vR01a,vR01b) rot_metric(R02,vR02a,vR02b)])

%generate vertical and horizontal vectors
R0=eye(3);
XVec=[eye(2); zeros(1,2)]*randn(2,1);
YVec=[0;0;1];
X=rot_hat(R0,XVec);
Y=rot_hat(R0,YVec);

RY=rot_geodFun(R0,Y);
RX=rot_geodFun(R0,X);

%check if a left-invariant field is basic along a geodesic by computing its
%norm after the pushforward

XRast=@(R) Dpi(R)*XVec;
XRastNorm=@(R) sphere_metric(pi(R),XRast(R),XRast(R));
R0t=rot_randGeodFun(R0);

%funPlot(@(t) XRastNorm(R0t(t)),'angle')

%check if the pushforward of tangent vectors at equivalent rotations
%produce the same result

disp([XRast(R01) XRast(R02)])




% %apply X to f
% fRX=@(t) fb(RX(t));
% Xfb=@(R) rot_metric(R,rot_hat(R,XVec),gradfb(R));
% %funCheckDer(fRX, @(t) Xfb(RX(t)),'angle')
% 
% %plot f and Xf along RY
% subplot(2,1,1)
% funPlot(@(t) fb(RY(t)));
% subplot(2,1,2)
% funPlot(@(t) Xfb(RY(t)));


