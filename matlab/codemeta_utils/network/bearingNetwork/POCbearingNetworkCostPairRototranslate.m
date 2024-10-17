function POCbearingNetworkCostPairRototranslate
%test function and its derivative along a bent radial line (linear+circular
%motion)

global funs
resetRands()

funs=bearingCostFunctions('cosine');

w=1;
v=[0;-1];
s=0.8;

xOffset=[0;0];
xCirc=@(t) [cos(w*t);sin(w*t)];
xTransl=@(t) s*t*v;

dxCirc=@(t) [-w*sin(w*t);w*cos(w*t)];
dxTransl=@(t) s*v;

x=@(t) xOffset+xCirc(t)+xTransl(t);
dx=@(t) dxCirc(t)+dxTransl(t);
X=[0;0];
YGoal=bearingCompute(x(0),X);


f1=@(t) costAndDerSingle(x(t),X,dx(t),YGoal);
t=linspace(0,pi,100);

%show cost along the line
figure(1)
check_der(f1,'function',t)
legend('Location','SouthEast')

%show the cost level sets, gradient and trajectory
sceneData.XGoal=x(0);
sceneData.XLandmarks=[0;0];
sceneData.XLandmarksEval=[0;0];
sceneData.L=3;
displayOpts={};%{'optsQuiver',{'lineWidth',2},'optsContour',{'lineWidth',2},...
    %'optsX',{'markerSize',10,'markerFaceColor','b'},};
figure(2)
bearingCostDisplay(sceneData,funs,displayOpts{:})
hold on
xt=funEval(x,t);
plotPoints(xt,'-')
hold off


function [c,dc]=costAndDerSingle(XEval,X,dXEval,YGoal)
global funs
[YEval,nYEval]=bearingCompute(XEval,X);
[c,gradc]=bearingCostGeneral(YEval,YGoal,nYEval,funs);
dc=gradc'*dXEval;
