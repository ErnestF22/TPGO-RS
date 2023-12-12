function POCCarControlCosineCost
seed=13;
flagSaveScreenshot=true;

resetRands(seed);


L=10;
NGrid=100;
NX=5;
TFinal=50;

optsX={'MarkerSize',10};

X=L*rand(2,NX);
%XGoal=[L/2;L/2];
XGoal=L*rand(2,1);

%x0=[L*rand;L*rand;2*pi*rand];
x0=[5;5];

[t,x]=ode45(@(t,x) control(x,X,XGoal), [0 TFinal], x0);
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
plot(X(1,:),X(2,:),'b*',optsX{:})
plot(XGoal(1),XGoal(2),'g*',optsX{:})
plot(x0(1),x0(2),'r*')
plot(x(1,:),x(2,:),'k')
hold off

if flagSaveScreenshot
    savefigure(['../docs/cartControlExampleCosineCost_' num2str(seed)],'epsc',0.5*[640 480])
end

function [c,dc]=cost(XEval,X,XGoal)
NX=size(X,2);
S=[0 -1; 1 0];
uVec=ones(1,NX);
YGoal=cnormalize(X-XGoal*uVec);
YEval=cnormalize(X-XEval*uVec);
c=sum((1-sum(YGoal.*YEval)).^2);
if nargout>1
    YEvalOrth=S*YEval;
    ip=sum(YEvalOrth.*YGoal);
    dc=-sum(([1;1]*ip).*YEvalOrth,2);
end



function u=control(x,X,XGoal)
[~,u]=cost(x,X,XGoal);