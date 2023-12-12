function POCCartControlCost
seed=12;
resetRands(seed);
flagSaveScreenshot=false;


L=10;
NGrid=100;
NX=10;
TFinal=50;
k=1;

optsX={'MarkerSize',10};

X=L*rand(2,NX);
%XGoal=[L/2;L/2];
XGoal=L*rand(2,1);
uVec=ones(1,NX);

%x0=[L*rand;L*rand;2*pi*rand];
x0=[5;5;0];


[t,x]=ode45(@(t,x) fModel(x,control(x,X,XGoal,k)), [0 TFinal], x0);
x=x';

z=linspace(0,L,NGrid);
[gridX,gridY]=meshgrid(z,z);
for iX=1:NGrid
    for iY=1:NGrid
        XEval=[gridX(iX,iY);gridY(iX,iY)];
        C(iX,iY)=cartBearingCost(XEval,X,XGoal);
    end
end
figure(1)
contour(gridX,gridY,C,20)
hold on
plot(X(1,:),X(2,:),'r*',optsX{:})
plot(XGoal(1),XGoal(2),'g*',optsX{:})
plot(x0(1),x0(2),'b*')
plot(x(1,:),x(2,:))
hold off

if flagSaveScreenshot
    savefigure(['cartControlExample_' num2str(seed)],'epsc')
end

figure(2)
l=evalLyapunov(x,X,XGoal,k);
plot(l);

function dx=fModel(x,u)
ctheta=@(x) cos(x(3));
stheta=@(x) sin(x(3));
dx=[ctheta(x) 0; stheta(x) 0; 0 1]*u;

function l=evalLyapunov(x,X,XGoal,k)
nCandidate=2;
NPoints=size(x,2);
l=NaN(1,NPoints);
for ix=1:NPoints
    xEval=x(:,ix);
    [u,thetaRef]=control(xEval,X,XGoal,k);
    dx=fModel(xEval,u);
    switch nCandidate
        case 1
            dx=dx(1:2);
            gradc=cartBearingCost(xEval(1:2),X,XGoal);
            dxRef=-gradc;
            e=dx-dxRef;
        case 2
            e=cos(xEval(3)-thetaRef);
    end
    l(ix)=0.5*(e'*e);
end

function [u,thetaRef]=control(x,X,XGoal,k)
[~,gradc,Hessc]=cartBearingCost(x(1:2),X,XGoal);
fRef=-gradc;
dfRef=-Hessc*gradc;
fRefNormSq=sum(fRef.^2);
thetaRef=atan2(fRef(2),fRef(1));

u=[
    [cos(x(3)); sin(x(3))]'*fRef(1:2);...%sqrt(fRefNormSq);...
    (fRef(1)*dfRef(2)-dfRef(1)*fRef(2))/fRefNormSq-k*angleDiff(x(3),thetaRef)...
    ];
