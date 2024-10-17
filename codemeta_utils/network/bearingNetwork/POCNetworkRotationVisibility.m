function POCNetworkRotationVisibility
resetRands()
y0=[0;1];
TFinal=1;
funs=bearingCostFunctions('barrier',cos(pi/3),cos(pi/2));
nyMax=5;
Rtheta=@(t) rot_exp(eye(2),rot_hat(eye(2),t));

E=[ 1 2;
    2 3;
    3 1];
E=[E;fliplr(E)];


R0=cat(3,Rtheta(-pi/4),Rtheta(1/4*pi),Rtheta(pi));
x0=[0 0 1;
    0 2 3];

g0=invg(RT2G(R0,x0,'compact'));

NNodes=size(x0,2);
dg=@(t,g) closedLoop(E,g,y0,funs,nyMax);
cEval=@(g) cost(E,g,y0,funs,nyMax);
dg0=dg(0,g0);
c0=cEval(g0);

[t,z]=ode45(dg, [0 TFinal], g0);

Nt=length(t);
z=z';
g=reshape(z,2,3,NNodes,[]);
gFinal=g(:,:,:,end);
[RFinal,xFinal]=G2RT(gFinal);

bearingNetworkPlot_edges(x0,E,'b');
hold on
bearingNetworkPlot_nodes(x0,'b');
bearingNetworkPlot_cameras(R0,x0);
bearingNetworkPlot_edges(xFinal,E,'r');
bearingNetworkPlot_nodes(xFinal,'r');
bearingNetworkPlot_cameras(RFinal,xFinal);
hold off

% [RFinal,TFinal]=G2RT(gFinal);
% 
% %compute cost
% c=zeros(Nt,1);
% dc=zeros(Nt,1);
% for it=1:Nt
%     [c(it),dc(it)]=cost(g(:,:,it),xLandmarks,y0,funs,nyMax);
% end
% 
% figure(1)
% plot(xLandmarks(1,:),xLandmarks(2,:),'r*')
% hold on
% draw2dcamera(R0,T0)
% draw2dcamera(RFinal,TFinal)
% hold off
% 
% figure(2)
% semilogy(t,c)
% 
% figure(3)
% plot(t,dc)

function dg=closedLoop(E,g,y0,funs,nyMax)
flagUseControl=false;
g=reshape(g,2,3,[]);
[R,T]=G2RT(g);
[y,ny]=measurements(E,T);
if flagUseControl
    %uVec=bearingControlRotationVisibility(R,y,nyMax,y0,funs);
else
    uVec=-bearingNetworkCostRotationVisibilityGradient(E,R,y,ny,y0,funs);
end

dg=rotdrd_hat(g,uVec);
dg=dg(:);

function [phi,dphi]=cost(E,g,y0,funs,nyMax)
g=reshape(g,2,3,[]);
[R,T]=G2RT(g);
[y,ny]=measurements(E,T);
[phi,gradPhi]=bearingNetworkCostRotationVisibility(E,R,y,ny,y0,funs);
%uVec=bearingControlRotationVisibility(R,y,nyMax,y0,funs);
%dphi=uVec'*gradPhi;

function [y,ny]=measurements(E,x)
[y,ny]=bearingNetworkComputeBearings(x,E);

