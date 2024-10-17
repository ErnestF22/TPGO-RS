function groundHomographyFromRT_test

R1=rot([pi/2;0;0]);
T1=-R1*[0;-2;2];

H=groundHomographyFromRT(R1,T1,'poses');

%test point on the XY plane
NX=100;
X=[[1;1;0] [2*rand(2,NX);zeros(1,NX)]];

subplot(1,3,1)
draw3dcameraFromPose(R1,T1)
hold on
plot3(X(1,:),X(2,:),X(3,:),'b.')
hold off
axis square
axis equal
grid on

x1=projectFromRT(R1,T1,X);
l11=[-1/x1(1,1);0];
l12=[0;-1/x1(2,1)];
xl11=projectTo2DLine(x1,l11);
xl12=projectTo2DLine(x1,l12);

subplot(1,3,2)
plot(x1(1,:),x1(2,:),'b.')
hold on
plot(xl11(1,:),xl11(2,:),'g.')
plot(xl12(1,:),xl12(2,:),'g.')
hold off
axis([-1 1 -1 1])
view(0,-90)

xGround=homographyApply(H,x1,'inverse');
xGroundl1exp=homographyApply(H,xl11,'inverse');
xGroundl2exp=homographyApply(H,xl12,'inverse');

l11in0=homographyApply(H,l11,'line','inverse');
l12in0=homographyApply(H,l12,'line','inverse');
xGroundl1=projectTo2DLine(xGround,l11in0);
xGroundl2=projectTo2DLine(xGround,l12in0);

subplot(1,3,3)
plot(xGround(1,:),xGround(2,:),'b.')
hold on
plot(X(1,:),X(2,:),'ro')
plot(xGroundl1(1,:),xGroundl1(2,:),'g.')
plot(xGroundl2(1,:),xGroundl2(2,:),'g.')
plot(xGroundl1exp(1,:),xGroundl1exp(2,:),'mo')
plot(xGroundl2exp(1,:),xGroundl2exp(2,:),'mo')
hold off

function y=projectTo2DLine(x,l)
NX=size(x,2);
u=ones(1,NX);
A=[[l;1] [0;0;1]];
y=x-A(1:2,:)*((A'*A)\(A'*[x;u]-[0;1]*u));
