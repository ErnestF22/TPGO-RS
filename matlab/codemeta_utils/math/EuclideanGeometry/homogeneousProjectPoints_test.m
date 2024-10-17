function homogeneousProjectPoints_test
NX=50;
drawOpts={'side',2,'style','g'};

NRowSubplot=2;
NColSubplot=3;
l2D1=[1;1;1];
l2D2=[-1;1;0];
l3D1=[1;1;0;1];
l3D2=[1;-1;0;0];
l3D3=[0;-1;1;1];
x2D=randn(2,NX);
x3D=randn(3,NX);


subplot(NRowSubplot,NColSubplot,1)
y=homogeneousProjectPoints(x2D,l2D1);
plotxy(x2D,y)
hold on
draw2dLine(l2D1,drawOpts{:})
hold off
title('2D points on 2D line')

subplot(NRowSubplot,NColSubplot,2)
y=homogeneousProjectPoints(x2D,[l2D1 l2D2]);
plotxy(x2D,y)
title('2D points on intersection of 2D lines')

subplot(NRowSubplot,NColSubplot,4)
y=homogeneousProjectPoints(x3D,l3D1);
plot3xy(x3D,y)
hold on
draw3dPlane(l3D1,drawOpts{:})
hold off
title('3D points on 3D plane')

subplot(NRowSubplot,NColSubplot,5)
y=homogeneousProjectPoints(x3D,[l3D1 l3D2]);
plot3xy(x3D,y)
hold on
draw3dLine([l3D1 l3D2],drawOpts{:})
hold off
title('3D points on intersection of two 3D planes')

subplot(NRowSubplot,NColSubplot,6)
y=homogeneousProjectPoints(x3D,[l3D1 l3D2 l3D3]);
plot3xy(x3D,y)
title('3D points on intersection of three 3D planes')


function plotxy(x,y)
plot(x(1,:),x(2,:),'bx',y(1,:),y(2,:),'r+',0,0,'bo')
axis equal
axis square
grid on


function plot3xy(x,y)
plot3(x(1,:),x(2,:),x(3,:),'bx',y(1,:),y(2,:),y(3,:),'r+',0,0,0,'bo')
axis equal
axis square
grid on


