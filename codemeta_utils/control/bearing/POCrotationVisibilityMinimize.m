function POCrotationVisibilityMinimize
resetRands()
y0=[0;1];
NLandmarks=5;
offset=[5;7];
xLandmarks=randn(2,NLandmarks)+offset*ones(1,NLandmarks);
%xLandmarks=zeros(2,NLandmarks)+offset*ones(1,NLandmarks);
T0=[0;2];
nyMax=20;

funs=bearingCostFunctions('barrier',cos(pi/6),0);
%funs=bearingCostFunctions('cosine');

%theta0=-pi/4;
theta0=0;
R0=rot_exp(eye(2),rot_hat(eye(2),theta0));
TFinal=30;

[y,ny]=bearingCompute(T0,xLandmarks);
phi=bearingCostRotationVisibility(R0,y,ny,y0,funs);
disp('Initial cost')
disp(phi)


dg=@(t,g) closedLoop(g,xLandmarks,y0,funs,nyMax);
g0=[R0 T0];
[t,z]=ode45(dg, [0 TFinal], g0);

Nt=length(t);
z=z';
g=reshape(z,2,3,[]);
gFinal=g(:,:,end);
[RFinal,TFinal]=G2RT(gFinal);

%compute cost
c=zeros(Nt,1);
dc=zeros(Nt,1);
for it=1:Nt
    [c(it),dc(it)]=cost(g(:,:,it),xLandmarks,y0,funs,nyMax);
end

figure(1)
plot(xLandmarks(1,:),xLandmarks(2,:),'r*')
hold on
draw2dcameraFromAxesAndCenter(R0,T0)
draw2dcameraFromAxesAndCenter(RFinal,TFinal)
hold off
axis equal

figure(2)
semilogy(t,c)

figure(3)
plot(t,dc)

function dg=closedLoop(g,xLandmarks,y0,funs,nyMax)
flagUseControl=true;
g=reshape(g,2,3);
[R,T]=G2RT(g);
[y,ny]=measurements(T,xLandmarks);
if flagUseControl
    uVec=bearingControlRotationVisibility(R,y,nyMax,y0,funs);
else
    uVec=-bearingCostRotationVisibilityGradient(R,y,ny,y0,funs);
end

dg=rotdrd_hat(g,uVec);
dg=dg(:);

function [phi,dphi]=cost(g,xLandmarks,y0,funs,nyMax)
g=reshape(g,2,3);
[R,T]=G2RT(g);
[y,ny]=measurements(T,xLandmarks);
[phi,gradPhi]=bearingCostRotationVisibility(R,y,ny,y0,funs);
uVec=bearingControlRotationVisibility(R,y,nyMax,y0,funs);
dphi=uVec'*gradPhi;

function [y,ny]=measurements(x,xLandmarks)
[y,ny]=bearingCompute(x,xLandmarks);

