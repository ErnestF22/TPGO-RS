function POCHouseholderRotationDerivative
D=3;
%x0=randn(D,1);
x0=[1;zeros(D-1,1)];
[x1,dx1]=real_randGeodFun(x0,'speed','quadratic');
[x2,dx2]=real_randGeodFun(-x0,'speed','quadratic');


%check_der(@(t) normDerAndDDer(x1(t),dx1(t),dx1(1)));
%check_der(@(t) funDer(x1(t),x2(t),dx1(t),dx2(t)),'function',linspace(-1,1,101))
check_der(@(t) funDerVec(x1(t),x2(t),dx1(t),dx2(t),dx1(1),dx2(1)))
check_der(@(t) derAndDDer(x1(t),x2(t),dx1(t),dx2(t),dx1(1),dx2(1)))

% t=1e-8;
% [R0,dR0]=funDer(x1(t),x2(t),dx1(t),dx2(t));
% disp(R0)


function [R0,dR0]=funDer(x1,x2,dx1,dx2)
D=size(x1,1);
I=eye(D);
nx1=norm(x1);
x1p=x1/nx1;
nx2=norm(x2);
x2p=x2/nx2;

Dx1p=(I-(x1p*x1p'))/nx1;
% dx1p=Dx1p*dx1;
Dx2p=(I-(x2p*x2p'))/nx2;
% dx2p=Dx2p*dx2;

v=x1p+x2p;
nv=norm(v);
vp=v/nv;

Dvp=(I-(vp*vp'))/nv;


R0=2*(vp*vp')-I;

dvp=Dvp*[Dx1p Dx2p]*[dx1;dx2];
%dR0=2*(vp*dvp'+dvp*vp');
% Dx(:,:,1)=Dvp*Dx1p;
% Dx(:,:,2)=Dvp*Dx2p;
% dx=[dx1 dx2];
% dR0=zeros(3);
% for i=1:2
%     dR0=dR0+2*(vp*dx(:,i)'*Dx(:,:,i)'+Dx(:,:,i)*dx(:,i)*vp');
% end
%dR0=R0*hat3(-2*hat3(vp)*dvp);
dR0=rot_hat(R0,-2*hat3(vp)*Dvp*[Dx1p Dx2p]*[dx1;dx2]);

function [R0,dR0]=funDerVec(x1,x2,dx1,dx2,ddx1,ddx2)
[dR0v,ddR0v,R0]=derAndDDer(x1,x2,dx1,dx2,ddx1,ddx2);
dR0=R0*hat3(dR0v);

function [dR0v,ddR0v,R0]=derAndDDer(x1,x2,dx1,dx2,ddx1,ddx2)
x1p=cnormalize(x1);
x2p=cnormalize(x2);
dx1p=cnormalizeDiff(x1,dx1);
ddx1p=cnormalizeDDiff(x1,dx1,ddx1);
dx2p=cnormalizeDiff(x2,dx2);
ddx2p=cnormalizeDDiff(x2,dx2,ddx2);

v=x1p+x2p;
dv=dx1p+dx2p;
ddv=ddx1p+ddx2p;

vp=cnormalize(v);
dvp=cnormalizeDiff(v,dv);
ddvp=cnormalizeDDiff(v,dv,ddv);

dR0v=-2*hat3(vp)*dvp;
ddR0v=-2*hat3(vp)*ddvp;

if nargout>2
    R0=2*vp*vp'-eye(size(x1,1));
end

function [dxp,ddxp]=normDerAndDDer(x,dx,ddx)
I=eye(size(x,1));
nx=norm(x);
xp=x/nx;

Dxp=(I-(xp*xp'))/nx;
dxp=Dxp*dx;

%ddxp=-1/nx^2*xp'*dx*(I-(xp*xp'))*dx-1/nx*(dxp*xp'+xp*dxp')*dx+Dxp*ddx;
%ddxp=-1/nx*dxp*xp'*dx-1/nx*(dxp*xp'+xp*dxp')*dx+Dxp*ddx;
ddxp=-1/nx*(2*dxp*xp'+xp*dxp')*dx+Dxp*ddx;
