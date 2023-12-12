function POCRotationFollowField
resetRands(1)
d=3;
N=2;
funsName='cosine';
alpha=1;
TFinal=2.2;
R0=rot_randn();

xLandmarks=randn(d,N);
v0=randn(d,1);
v1=randn(d,1);
v2=randn(d,1);
switch 2
    case 1
        x=@(t) v0+t*v1+0.5*t^2*v2;
        dx=@(t) v1+t*v2;
        ddx=@(t) v2;
    case 2
        x=@(t) v0+sin(t)*v1;
        dx=@(t) cos(t)*v1;
        ddx=@(t) -sin(t)*v1;
end
xg=x(0);
yg=bearingCompute(xg,xLandmarks);

funs=bearingCostFunctions(funsName);

dR=@(t,R) closedLoop(R,x(t),dx(t),ddx(t),xLandmarks,yg,funs,alpha);

%check_der(@(t) FReference(x(t),dx(t),ddx(t),xLandmarks,yg,funs,alpha),'function')

[t,R]=ode45(dR, [0 TFinal], R0(:));
R=reshape(R',3,3,[]);

 
Nt=length(t);
e=zeros(Nt,1);
de=zeros(Nt,1);
for it=1:Nt
    [e(it),de(it)]=errAndDer(R(:,:,it),x(t(it)),dx(t(it)),ddx(t(it)),xLandmarks,yg,funs,alpha);
end

figure(1)
subplot(2,1,1)
plot(t,e)
subplot(2,1,2)
plot(t,de)
hold on
plot((t(1:end-1)+t(2:end))/2,diff(e)./diff(t),'rx')
hold off
    

function dR=closedLoop(R,x,dx,ddx,xLandmarks,yg,funs,alpha)
[nRef,dnRef]=FReference(x,dx,ddx,xLandmarks,yg,funs,alpha);
dR=controlRVec(R,nRef,dnRef);

function [e,de]=errAndDer(R,x,dx,ddx,xLandmarks,yg,funs,alpha)
e3=[0;0;1];

[nRef,dnRef]=FReference(x,dx,ddx,xLandmarks,yg,funs,alpha);
uR=controlR(R,nRef,dnRef);

e=1-e3'*R*nRef;
de=-e3'*R*hat(uR)*nRef-e3'*R*dnRef;

function [y,dy,ddy]=measurements(x,dx,ddx,xLandmarks)
[y,ny]=bearingCompute(x,xLandmarks);
dy=bearingComputeDerivative(dx,y,ny);
ddy=bearingComputeSecondDerivative(dx,ddx,y,ny);

function [nRef,dnRef]=FReference(x,dx,ddx,xLandmarks,yg,funs,alpha)
% nRef=[1;0;0];
% dnRef=zeros(3,1);
[y,dy,ddy]=measurements(x,dx,ddx,xLandmarks);
Fref=bearingDynamicControlDirect(y,dy,yg,funs,alpha);
dFref=bearingDynamicControlDirectDerivative(y,dy,ddy,yg,funs,alpha);
[nRef,nnRef]=cnormalize(Fref);
dnRef=(eye(3)-nRef*nRef')*dFref/nnRef;

function n=nActual(R)
e3=[0;0;1];
n=R'*e3;

function uR=controlR(R,nRef,dnRef)
kR=0.5;
e3=[0;0;1];
nR=hat(nRef)*R'*e3;
vR=(kR+max(0,-e3'*R*dnRef))/(nR'*nR);
uR=vR*nR;


function dR=controlRVec(R,nRef,dnRef)
R=reshape(R,3,3);
uR=controlR(R,nRef,dnRef);
dR=reshape(R*hat(uR),[],1);


