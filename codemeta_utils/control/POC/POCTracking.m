function POCTracking

TFinal=5;
k=6;
x0=2;

xr=@(t) sin(2*pi*t);
%xr=@(t) 1;
tr=linspace(0,TFinal);

e=@(x,t) x-xr(t);

dx1=@(t,x) -k*e(x,t);
dx2=@(t,x) k*controlSwitching(e(x,t));
dx3=@(t,x) k*controlSuperTwisting(x,xr(t));

%plotfun(@controlSwitching,linspace(-1,1,100))
dt=0.01;
[t1,x1]=odeEuler(dx1,[0 TFinal],x0,dt);
[t2,x2]=odeEuler(dx2,[0 TFinal],x0,dt);
[t3,x3]=odeEuler(dx3,[0 TFinal],[x0;0],dt);
dxEval1=evalControl(t1,x1,dx1);
dxEval2=evalControl(t2,x2,dx2);
dxEval3=evalControl(t3,x3,dx3);


figure(1)
subplot(3,1,1)
plot(tr,xr(tr),'g',t1,x1,'b',t2,x2,'r',t3,x3(:,1),'c')
subplot(3,1,2)
plot(t1,abs(xr(t1)-x1),'b',t2,abs(xr(t2)-x2),'r',t3,abs(xr(t3)-x3(:,1)),'c')
title('Errors')
subplot(3,1,3)
plot(t1,dxEval1,'b',t2,dxEval2,'r',t3,dxEval3(:,1),'c')
title('Control')


function dx=controlSwitching(e)
delta=0.01;
ne=abs(e);
if ne>delta
    dx=-e/ne;
else
    dx=-e/delta;
end

function dx=controlSuperTwisting(x,xr)
flagLinear=false;
x=shiftdim(x);
xr=shiftdim(xr);
k1=3;
k3=1;
if flagLinear
    dx(1,:)=-k1*(x(1,:)-xr)+x(2,:);
else
    dx(1,:)=-k1*sqrtSign(x(1,:)-xr)+x(2,:);
end
dx(2,:)=-k3*sign(x(1,:)-xr);

function dxEval=evalControl(t,x,dx)
dxEval=zeros(size(x));
Nt=length(t);
for it=1:Nt
    dxEval(it,:)=dx(t(it),x(it,:));
end