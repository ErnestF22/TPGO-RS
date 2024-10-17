function multiCalibration_realData_visualize
figure(1)
load calibrationDataset/april_vel_processed.mat
subplot(3,2,1)
draw3dcameraFromRT(eye(3),zeros(3,1),'scale',0.5)
hold on
plotPlanes(normalCamera,distanceCamera);
hold off
axis equal
view(0,0)
title('Camera and planes')
subplot(3,2,2)
draw3dcameraFromRT(eye(3),zeros(3,1),'scale',0.5,'shape','velodyne')
hold on
plotPlanes(normalVelodyne,distanceVelodyne);
hold off
axis equal
view(0,90)
title('Velodyne and planes')

load calibrationDataset/april_hok_processed.mat
subplot(3,2,3)
draw3dcameraFromRT(eye(3),zeros(3,1),'scale',0.5)
hold on
plotPlanes(normalCamera,distanceCamera);
hold off
axis equal
view(0,0)
title('Camera and planes')
subplot(3,2,4)
NPoints=size(hokuyo_endpts2D,2);
draw3dcameraFromRT(eye(3),zeros(3,1),'shape','hokuyo')
hold on
plotPoints([zeros(1,NPoints);  hokuyo_endpts2D])
hold off
title('Hokuyo and points')
view(90,0)
axis equal
axis(2.5*[-1 1 -1 1 -1 1])

load calibrationDataset/hok_vel_processed.mat
subplot(3,2,5)
draw3dcameraFromRT(eye(3),zeros(3,1),'scale',0.5,'shape','velodyne')
hold on
plotPlanes(normalVelodyne,distanceVelodyne);
hold off
axis equal
view(0,90)
title('Velodyne and planes')
subplot(3,2,6)
NPoints=size(hokuyo_endpts2D,2);
draw3dcameraFromRT(eye(3),zeros(3,1),'shape','hokuyo')
hold on
plotPoints([zeros(1,NPoints); hokuyo_endpts2D])
hold off
title('Hokuyo and points')
view(90,0)
axis equal
axis(2.5*[-1 1 -1 1 -1 1])



function plotPlanes(normal,distance)
NPlanes=size(normal,2);
for iPlane=1:NPlanes
    draw3dPlane([normal(:,iPlane);-distance(iPlane)],'center',zeros(3,1),'side',0.5)
end
